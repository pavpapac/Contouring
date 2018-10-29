# -*- coding: utf-8 -*-

import math
from pathlib import Path

import click
import numpy as np
import pydicom as dicom
from shapely.geometry import Polygon


def ReadContourStructures(DCM):
    """Reads a dicom (.dcm) file and returns a dictionary with contour information.

    The information including colour, name and contour data Contour data
    include all the (x,y,z) data points used to build the contour
    """
    contours = []

    ds = dicom.read_file(DCM)
    num_of_contours = len(ds.ROIContourSequence)

    # Extract contour color,name and data points
    for i in range(num_of_contours):
        contour = {}
        contour['color'] = ds.ROIContourSequence[i].ROIDisplayColor
        # contour['number']=ds.ROIContourSequence[i].RefdROINumber
        contour['name'] = ds.StructureSetROISequence[i].ROIName
        contour['contour data'] = [
            s.ContourData for s in ds.ROIContourSequence[i].ContourSequence]
        contours.append(contour)
    return contours


def FindContour(contours, str_name):
    """Finds a structure with name str_name in the contour dictionary.

    First the user asks for a contour and we verify that the name exists in
    the contour list.
    Outputs ReadContourStructures
    """

    name_agrees = [str_name == contour['name'] for contour in contours]

    # If the name is found return a positive
    # TODO Trow exception on which index to use
    if sum(name_agrees) == 1:
        contour_ind = name_agrees.index(1)
        print('The contour', str_name, ' has been selected\n')
        return contour_ind
    else:
        print('Name not found. Restart and provide a structure name from the list below\n')
        print([contour['name'] for contour in contours])
        return [contour['name'] for contour in contours]


def ExtractContourPoints(contours, contour_ind):
    """ Transform the contour matrix.

    Reads from the contour dictionary the x,y,z coordinates for each the
    structure slice and returns them in a new dictionary coordinates. The
    output will be later used to create a Polygon object
    """

    x_data = []
    y_data = []
    z_data = []
    coord = {}

    contour_data = contours[contour_ind]['contour data']

    for contour in contour_data:
        # first read all x,y,z values for each contour slice. Cordinates
        # are stored sequentially you need a step of 3 to go from x0 -> x1 etc
        x_temp = contour[0::3]
        y_temp = contour[1::3]
        z_temp = contour[2::3]

        # Convert to float
        x_temp = [float(elem) for elem in x_temp]
        y_temp = [float(elem) for elem in y_temp]
        z_temp = [float(elem) for elem in z_temp]

        # append to our data table
        x_data.append(x_temp)
        y_data.append(y_temp)
        z_data.append(z_temp)

    # Now convert the x,y,z lists to arrays and store them to a dictionary
    coord['x'] = np.asarray(x_data)
    coord['y'] = np.asarray(y_data)
    coord['z'] = np.asarray(z_data)

    return coord


def CreatePolygonList(coord):
    """For each contour of the structure creates a Polygon object.

    Using the shapely.geometry library to build the object, each contour slice
    will become a Polygon object. Creates and returns a dictionary
    (poly_contours) with key element: 'polygons' which has a list of the
    contour/Polygon objects per slice and 'z-slice' which has the CT slice z
    index (mm) that each contour exists.
    """

    poly_contourlist = []
    poly_z = []
    poly_contours = {}

    # For each structure slice extract the x,y coordinates and return them
    # in a list of tuples (x,y). This is then used to create the Polygon object
    # The z coordinates are always the same for each slice
    # so just keep the first element.

    for xsl, ysl, zsl in zip(coord['x'], coord['y'], coord['z']):
        poly_xycrds = []
        for x, y in zip(xsl, ysl):
            poly_xycrds.append((x, y))
        poly_z.append(zsl[0])
        poly_contour = Polygon(poly_xycrds)
        poly_contourlist.append(poly_contour)

    # Now create a dictionary with all the polygon objects for all slices
    # with their respective z-slice position

    poly_contours['polygons'] = poly_contourlist
    poly_contours['z-slice'] = poly_z

    return poly_contours


def CreatePolygon(dcm, roi_name):
    """ Create a Polygon from the dicom file.
    This function performs the follwoing sequentially:
    * reads a dicom file and extracts the structure
    * returns the (x,y,z) coordinates dictionary
    * returns the poly_contours dictionary which includes the list of all
    polygon contours per slice and the z-slices where the contours exist.

    This function is the main building block for the rest of the code since all
    the analysis and metrics will be done on the poly_contours
    """

    # First read all contour structures
    contours = ReadContourStructures(dcm)

    # Now extract the index for the structure requested
    contour_ind = FindContour(contours, roi_name)

    # now extract contour coordinates for the structure
    coord = ExtractContourPoints(contours, contour_ind)

    # Now return a list with all the contours slice by slice as polygon objects.
    # The polygons will be analyzed in the next section with the aid of the shapely
    # library
    poly_contours = CreatePolygonList(coord)

    return poly_contours


def MergeMultiContours(poly_contours):
    """Merges multi contours that exist in the same slice to one unified contour.

    This method is needed in the cases where multiple contours exist in the
    same slice. This also serves either as clean-up method for the existance of
    small dots left by the physician accidentally.
    """

    # First extract the z CT slice indices and polygon objects for the structure.
    z = poly_contours['z-slice']
    poly = poly_contours['polygons']
    ind = []

    # Then find the z-index that is equal to the next one. This means
    # that there are 2 contours sharing the same z index. Perform a union operation
    # on the polygons of that slice to create one unified contour.

    for i in range(len(z)-1):
        if z[i] == z[i+1]:
            poly[i+1] = poly[i+1].union(poly[i])
            ind.append(i)

    # Finally delete the extra contours and slice indices after unification.
    # We should only have unique z-slice indices from now on.

    for i in ind:
        del poly_contours['z-slice'][i], poly_contours['polygons'][i]

    return poly_contours, ind


def ExtractBoundGradChange(poly_contours):
    num_of_cnts = len(poly_contours['polygons'])
    #Initialize: boundaries
    xmin = []
    ymin = []
    xmax = []
    ymax = []

    # ...gradients (degrees)
    theta_xmin = []
    theta_ymin = []
    theta_xmax = []
    theta_ymax = []

    # ...gradient changes
    dtheta_xmin = []
    dtheta_ymin = []
    dtheta_xmax = []
    dtheta_ymax = []

    # Calculate the thickness of the CT slices (t=1.25 mm)
    t = abs(poly_contours['z-slice'][1]-poly_contours['z-slice'][0])

    # Extract x and y boundaries for each slice contour
    for polygon in poly_contours['polygons']:
        xmin.append(polygon.bounds[0])
        ymin.append(polygon.bounds[1])
        xmax.append(polygon.bounds[2])
        ymax.append(polygon.bounds[3])

    # Calculate the gradient for each x,y boundary of the contour
    for i in range(num_of_cnts-1):

        theta_xmin_temp = math.degrees(math.atan((xmin[i]-xmin[i+1])/t))
        theta_ymin_temp = math.degrees(math.atan((ymin[i]-ymin[i+1])/t))
        theta_xmax_temp = math.degrees(math.atan((xmax[i]-xmax[i+1])/t))
        theta_ymax_temp = math.degrees(math.atan((ymax[i]-ymax[i+1])/t))

        theta_xmin.append(theta_xmin_temp)
        theta_ymin.append(theta_ymin_temp)
        theta_xmax.append(theta_xmax_temp)
        theta_ymax.append(theta_ymax_temp)

    # Now calculate the the slice-to-slice gradient change (2nd gradient...)
    for i in range(num_of_cnts-2):
        dtheta_xmin.append(theta_xmin[i]-theta_xmin[i+1])
        dtheta_ymin.append(theta_ymin[i]-theta_ymin[i+1])
        dtheta_xmax.append(theta_xmax[i]-theta_xmax[i+1])
        dtheta_ymax.append(theta_ymax[i]-theta_ymax[i+1])

    dtheta_xmin = np.asarray(dtheta_xmin)
    dtheta_ymin = np.asarray(dtheta_ymin)
    dtheta_xmax = np.asarray(dtheta_xmax)
    dtheta_ymax = np.asarray(dtheta_ymax)

    # Calculate the mean absolute gradient change for each boundary.

    dtheta_xmin_av = np.mean(abs(dtheta_xmin))
    dtheta_ymin_av = np.mean(abs(dtheta_ymin))
    dtheta_xmax_av = np.mean(abs(dtheta_xmax))
    dtheta_ymax_av = np.mean(abs(dtheta_ymax))

    # Finally calculate the average of the four boundary gradient changes.
    DGrad_av = (dtheta_xmin_av+dtheta_ymin_av+dtheta_xmax_av+dtheta_ymax_av)/4

    return DGrad_av, t


def ExtractVolume(poly_contours):

    # First calculate the thickness of the CT slices
    t = abs(poly_contours['z-slice'][0]-poly_contours['z-slice'][1])
    # Now estimate the total volume
    vol = 0
    for polygon in poly_contours['polygons']:
        vol += t*polygon.area

    return vol


def ExtractIntersectionPoly(polyA, polyB):
    # Initialize
    polyABinter = {}
    polyABinter['polygons'] = []
    polyABinter['z-slice'] = []

    # First find the slices that are common (zA=zB). then calculate the intersection
    # of the contours and append to a new dictionary.

    for zA, pA in zip(polyA['z-slice'], polyA['polygons']):
        for zB, pB in zip(polyB['z-slice'], polyB['polygons']):
            if zA == zB:
                inter = pA.intersection(pB)
                polyABinter['polygons'].append(inter)
                polyABinter['z-slice'].append(zA)

    return polyABinter


def ExtractUnionPoly(polyA, polyB):
    # Initialize
    polyABunion = {}
    polyABunion['polygons'] = []
    polyABunion['z-slice'] = []

    for zA, pA in zip(polyA['z-slice'], polyA['polygons']):
        for zB, pB in zip(polyB['z-slice'], polyB['polygons']):
            if zA == zB:
                union = pA.union(pB)
                polyABunion['polygons'].append(union)
                polyABunion['z-slice'].append(zA)

    return polyABunion


def HausdorffMetrics(polyA, polyB):
    """ Calculates the Haussdorff distances for each contour slice.

    Hausdorff distance is defined as: maximum distance of the closest point
    between two contours.

    Returns the max, mean and std of all Hausdorff distances.
    """

    H = []

    # Now calculate the Hausdorff distance of the contours for each slice
    for zA, pA in zip(polyA['z-slice'], polyA['polygons']):
        for zB, pB in zip(polyB['z-slice'], polyB['polygons']):
            if zA == zB:
                htemp = pA.hausdorff_distance(pB)
                H.append(htemp)

    # Now calculate the max, mean and std of Hausdorff distances. Return it as a list.
    Hdmax = np.max(H)
    Hdmean = np.mean(H)
    Hdstd = np.std(H)

    return [Hdmax, Hdmean, Hdstd]


def CentroidDist(polyA, polyB):

    Cd = []

    for zA, pA in zip(polyA['z-slice'], polyA['polygons']):
        for zB, pB in zip(polyB['z-slice'], polyB['polygons']):
            if zA == zB:
                xA, yA = int(pA.centroid.x), int(pA.centroid.y)
                xB, yB = int(pB.centroid.x), int(pB.centroid.y)
                rAB = np.sqrt((xA-xB)**2+(yA-yB)**2)
                Cd.append(rAB)

    Cdmax = np.max(Cd)
    Cdmean = np.mean(Cd)
    Cdstd = np.std(Cd)

    return [Cdmax, Cdmean, Cdstd]


def TotalNumSlicesContoured(poly):

    minz = 10000
    maxz = -10000

    # Now calculate the thickness of the CT slices (in mm)
    t = abs(poly[0]['z-slice'][0]-poly[0]['z-slice'][1])

    # now find the min and max z-slice that has a contour from all the structures
    for p in poly:
        if minz > min(p['z-slice']):
            minz = min(p['z-slice'])
        if maxz < max(p['z-slice']):
            maxz = max(p['z-slice'])

    # Now calculate the number of total number of slices contoured (by any structure)

    num_tot_sl_contrd = int((maxz-minz)/t)+1

    return num_tot_sl_contrd


def Jaccard(poly):
    """calculates the intersection and union areas for all contours per slice.

    Given a list of poly_contours objects (output of CreatePolygon) calculates
    the intersection and union areas for all contours per slice. It extracts
    the Jaccard index defined as intersection area/union area per slice. It
    also calculates the average Jaccard index and the average Jaccard excluding
    all slices were at least one physician did not contour anything. The later
    helps us evaluate how well physicians perform in the slices they all
    contoured.
    """

    poly_inter = 0
    poly_union = 0
    jac_ind = []
    Jaccard_av = 0

    # first extract the intersection and union of each slice for the first two
    # polygons in the list. Then calculate the Jaccard index as the ratio
    # of (intersection area) / (union area).

    poly_inter = ExtractIntersectionPoly(poly[0], poly[1])
    poly_union = ExtractUnionPoly(poly[0], poly[1])

    for p in poly[2:]:
        poly_inter = ExtractIntersectionPoly(poly_inter, p)
        poly_union = ExtractUnionPoly(poly_union, p)

    common_zslices = len(poly_inter['z-slice'])

    for i, u in zip(poly_inter['polygons'], poly_union['polygons']):
        inter_area = i.area
        union_area = u.area
        jac_ind.append(inter_area/union_area)
        # weights.append(union_area)
        # sum_area+=union_area

    # Now calculate the total number of slices not contoured by at least one
    # physician and the percentage of slice contoured by all physicians.
    # Assign 0'z to allthe slcies that at least one physician did not contour anything
    # (maximum discrepancy)

    tot_zslices = TotalNumSlicesContoured(poly)
    not_common_zslices = tot_zslices-common_zslices
    perc_common_zslices = 100*(common_zslices/tot_zslices)
    jac_ind.extend([0]*not_common_zslices)

#    # Now calculate the weighted average Jaccard index. The weihgts are calculated
#    # as the ratio: contoured area per slice/total contoured area.
#
#    weights=[w/sum_area for w in weights]
#
#    for w,j in zip(weights,jac_ind):
#        Jaccard_wav+=w*j
#    Jaccard_wav=Jaccard_wav/sum(weights)

    # Now calclate the mean jaccard index for all slices. That includes also
    # the slices assigned 0 due to at least one physician not contoured it.
    # Calculate also the average including only the common contoured slices, that is
    # all the slices exclufing the zeros to see how well physicians did in the slices
    # where they all contoured

    Jaccard_av = np.mean(jac_ind)
    Jaccard_av_common = np.mean([ind for ind in jac_ind if ind != 0])

    return perc_common_zslices, jac_ind, Jaccard_av, Jaccard_av_common


def HausdorffDict(poly):
    """Returns the max, mean and std Haussdorf distance for all combinations

     Given a list of poly_contour objects returns the max, mean and std
    Hausdorff distance for all combinations (order selection does not matter
    here). For example calculates Haussdorf for structure 1 with 2,3..
    structure 2 with 3,4..etc etc INPUT: poly is a list of Polygon contour
    objects as created by the CreatePolygon() method. OUTPUT: Dictionary with
    derived hausdorff metrics (max,mean,std) for each comparison performed.
    Thus, the length of each key in the dictionary should be equal to the total
    number of combinations possible.
    """

    Haus_dict = {'max': [], 'mean': [], 'std': []}

    # We first perform a loop over all structures (poly). For each structure (pi) we will
    # calculate the Hausdorff distance to the remaining structures in the list
    # Return then the metrics (max,mean,std) in a dictionary.

    for pi in poly:
        for pj in poly[poly.index(pi)+1:]:

            Haus = HausdorffMetrics(pi, pj)
            Haus_dict['max'].append(Haus[0])
            Haus_dict['mean'].append(Haus[1])
            Haus_dict['std'].append(Haus[2])

    return Haus_dict


def CentroidDict(poly):

    Centroid_dict = {'max': [], 'mean': [], 'std': []}

    # We first perform a loop over all structures (poly). For each structure (pi) we will
    # calculate the distances between centroids to the remaining structures in the list
    # Return then the metrics (max,mean,std) in a dictionary.

    for pi in poly:
        for pj in poly[poly.index(pi)+1:]:

            C = CentroidDist(pi, pj)
            Centroid_dict['max'].append(C[0])
            Centroid_dict['mean'].append(C[1])
            Centroid_dict['std'].append(C[2])

    return Centroid_dict


@click.command()
@click.option('--folder',
              default='.',
              envvar='CONTOUR_PATH',
              help='Root path to the dicom')
@click.option('--files',
              type=click.Path(),
              multiple=True,
              help='Specific dicom files')
@click.argument('structure', type=click.STRING)
@click.argument('energy', type=click.INT)
def main(folder, files, structure, energy):
    """ The program compares the structures found in different images in a dicom file(s).
    """

    poly = []
    DGrad_av = []
    vol = []

    root = Path(folder)
    if not root.is_dir():
        raise click.ClickException('The folder provided is not accessible')

    dcm_list = []
    if files:
        dcm_list = [root / x for x in files]
    else:
        dcm_list = list(root.glob('*.dcm'))

    # First create a polygon object
    # Merge any multi polygons in the same slice
    # Calculate the average gradient change of the contour
    # Calculate the structure volume
    # Append the polygon lists, gradient changes and volumes for each structure to a poly list
    for dcm in dcm_list:
        if not dcm.is_file():
            click.echo(f'File {dcm} does not exists')
            continue
        poly_contours = CreatePolygon(dcm.name, structure)
        p, ind = MergeMultiContours(poly_contours)
        dg, t = ExtractBoundGradChange(p)
        v = ExtractVolume(p)
        poly.append(p)
        DGrad_av.append(dg)
        vol.append(v)

    if poly == []:
        raise click.ClickException('No contour(s) found in the files specified')

    # Calculate the Jaccard similarity index per slice (jac_ind), the average
    # Jaccard and the weighted average Jaccard.
    perc_common_zslices, jac_ind, Jaccard_av, Jaccard_av_common = Jaccard(poly)

    # Calculate the Haussdorf dictionary which includes max, mean and std of the
    # Hausdorff distance for all slices.
    Haus_dict = HausdorffDict(poly)

    # Calculate the centroid dictionary which includes the max,mean and std distances
    # between the contours for all slices
    Centroid_dict = CentroidDict(poly)

    click.echo('----------CONTOUR ANALYSIS----------')
    for dcm, v, g in zip(dcm_list, vol, DGrad_av):
        click.echo(f'Structure {structure} of dicom file {dcm.name} has volume of: {v} mm^3')
        click.echo(f'And an average gradient change of: {g} degrees')

    click.echo(f'The percentage of CT slices contoured by all physicians: {perc_common_zslices:0.2f} %')
    click.echo(f'The Jaccard similarity indices per slice: {jac_ind}')
    click.echo(f'The average Jaccard index: {Jaccard_av:0.2f}')
    click.echo(f'And for all the common slices contoured: {Jaccard_av_common:0.2f}')
    click.echo(f'The Hausdorff distance metrics: {Haus_dict}')
    click.echo(f'The centroid misalignment metrics: {Centroid_dict}')

    # Now prepare the data to be written in a textfile by converting numbers to
    # strings. Caclualte also here any specific metrics you prefer saving in the
    # file.
    Vav = str(np.mean(vol))
    Vstd = str(np.std(vol))
    Gav = str(np.mean(DGrad_av))
    Gstd = str(np.std(DGrad_av))
    Jav = str(Jaccard_av)
    Jav_cmn = str(Jaccard_av_common)

    Hmax = str(Haus_dict['max'][0])
    Hmean = str(Haus_dict['mean'][0])
    Hstd = str(Centroid_dict['std'][0])
    Cmax = str(Centroid_dict['max'][0])
    Cmean = str(Centroid_dict['mean'][0])
    Cstd = str(Centroid_dict['std'][0])

    # Now write the output values in the file line-by-line
    f = open('HNSCC.out', 'a+')
    f.write(str(energy)+','+Vav+','+Vstd+','+Gav+','+Gstd+','+Jav+','+Jav_cmn +
            ','+Hmax+','+Hmean+','+Hstd+','+Cmax+','+Cmean+','+Cstd+'\n')
    f.close()


if __name__ == '__main__':
    main()
