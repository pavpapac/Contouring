#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 19 09:45:22 2018

@author: pavlos
"""

from __future__ import absolute_import, division, print_function
import pydicom as dicom
import math
import numpy as np
from shapely.geometry import Polygon

def ReadContourStructures(DCM):
    # FUNCTION READCONTOURSTRUCTURES: Reads a dicom (.dcm) file and returns 
    #a dictionary with contour information including colour, name and contour data
    #Contour data include all the (x,y,z) data points used to build the contour
    contours=[]
    
    #First read the dicom file    
    ds=dicom.read_file(DCM)
    num_of_contours=len(ds.ROIContourSequence)
    
    #Extract contour color,name and data points
    for i in range(num_of_contours):
        contour={}
        contour['color']=ds.ROIContourSequence[i].ROIDisplayColor
        #contour['number']=ds.ROIContourSequence[i].RefdROINumber
        contour['name']=ds.StructureSetROISequence[i].ROIName
        contour['contour data']=[s.ContourData for s in ds.ROIContourSequence[i].ContourSequence]
        contours.append(contour)
    return contours

    
def FindContour(contours,str_name):
   #FUNCTION FINDCONTOUR: Finds a structure with name str_name in the contour 
   #dictionary. Output of function ReadContourStructures
   
      
   # First the user asks for a contour and we verify that the name exists in 
   # the contour list
   
   name_agrees=[str_name==contour['name'] for contour in contours]
   
   #If the name is found return a positive 
   if sum(name_agrees)==1:
        contour_ind=name_agrees.index(1)
        print('The contour', str_name,' has been selected\n')
        return contour_ind
   else:
        print('Name not found. Restart and provide a structure name from the list below\n') 
        print([contour['name'] for contour in contours])
        return [contour['name'] for contour in contours]
    

def ExtractContourPoints(contours,contour_ind):
       #FUNCTION EXTRACTCONTOURPOINTS: Reads from the contour dictionary the 
       #x,y,z coordinates for th structure selected and returns them in a new dictionary coord
       
       x_data=[]
       y_data=[]
       z_data=[]
       coord={}
       
       contour_data=contours[contour_ind]['contour data']
       
       for i in range(len(contour_data)):
           #first read all x,y,z values for each contour slice. Cordinates 
           #are stored sequentially you need a step of 3 to go from x0 -> x1 etc
           x_temp=contour_data[i][0::3]
           y_temp=contour_data[i][1::3]
           z_temp=contour_data[i][2::3]
           
           #Convert to float
           x_temp=[float(elem) for elem in x_temp]
           y_temp=[float(elem) for elem in y_temp]
           z_temp=[float(elem) for elem in z_temp]
           
           # append to our data table
           x_data.append(x_temp)
           y_data.append(y_temp)
           z_data.append(z_temp)
       
       #Now convert the x,y,z lists to arrays and store them to a dictionary
       coord['x']=np.asarray(x_data)
       coord['y']=np.asarray(y_data)
       coord['z']=np.asarray(z_data)
       
       return coord
        
def CreatePolygonList(coord):
        #FUNCTION CREATEPOLYGONLIST: For each contour of the structure 
        #creates a Polygon object using the shapely.geometry library. Thus each contour
        #will becomea Polygon object. Creates and returns a dictionary (poly_contours)
        #with key element: 'polygons' which has a list of the contour/Polygon objects per slice
        #and 'z-slice' which has the CT slice z index (mm) that each contour exists. 
        
        poly_list=[]
        poly_z=[]
        poly_contours={}
        
        #First get the number of contours we have in the structure
        num_contours=len(coord['x'])
        
        # for each contour slice (z) read the x,y coordinates and create a polygon 
        # oblject.
        for i in range(num_contours):
            num_crds=len(coord['z'][i])
            
            poly_crds=[]
            for j in range(num_crds):
                #Get the x,y coordinates for each z-slice
                poly_crds.append((coord['x'][i][j],coord['y'][i][j]))
            
            #Append the z-slice
            poly_z.append(coord['z'][i][j])
            #Create a Polygon object based on the coordinates
            poly_contour=Polygon(poly_crds)
            poly_list.append(poly_contour)
            
        #Create a dictionary with all the polygon objects with their 
        #respective z-position
        poly_contours['polygons']=poly_list
        poly_contours['z-slice']=poly_z
        
        return poly_contours
    
def CreatePolygon(dcm,roi_name):
    # FUNCTION CREATEPOLYGON: This function performs the follwoing sequentially:
    #first reads a dicom file and extracts the structure, then returns the (x,y,z)
    #coordinates dictionary and finally return the poly_contours dictionary which
    #includes the list of all polygon contours per slice and the z-slices where the
    #contours exist. This function is the main building block for the rest of the code
    #since all the analysis and metrics will be done on the poly_contours 
    
    #First read all contour structures
    contours=ReadContourStructures(dcm)
    
    #Now extract the index for the structure requested
    contour_ind=FindContour(contours,roi_name)
    
    #now extract contour coordinates for the structure 
    coord=ExtractContourPoints(contours,contour_ind)
    
    #Now return a list with all the contours slice by slice as polygon objects.
    #The polygons will be analyzed in the next section with the aid of the shapely
    # library
    poly_contours=CreatePolygonList(coord)
    
    return poly_contours

def MergeMultiContours(poly_contours):
    #FUNCTION MERGEMULTICONTOURS: Merges multi contours that exist in the same 
    #slice to one unified contour. This method is needed in the cases where 
    #multiple contours exist in the same slice. This also serves either as 
    #clean-up method for the existance of small dots left by the physician accidentally.
    
    #First extract the z CT slice indices and polygon objects for the structure.
    z=poly_contours['z-slice']
    poly=poly_contours['polygons']
    ind=[]
    
    #Then find the z-index that is equal to the next one. This means 
    #that there are 2 contours sharing the same z index. Perform a union operation
    #on the polygons of that slice to create one unified contour. 
    
    for i in range(len(z)-1):
      if z[i]==z[i+1]:
        poly[i+1]=poly[i+1].union(poly[i])
        ind.append(i)
        
    #Finally delete the extra contours and slice indices after unification. 
    #We should only have unique z-slice indices from now on. 
    
    for i in ind:
        del poly_contours['z-slice'][i],poly_contours['polygons'][i]
        
    return poly_contours,ind
    
        
def ExtractBoundGradChange(poly_contours):
    num_of_cnts=len(poly_contours['polygons'])
    #Initialize: boundaries
    xmin=[]
    ymin=[]
    xmax=[]
    ymax=[]
    
    #...gradients (degrees)
    theta_xmin=[]
    theta_ymin=[]
    theta_xmax=[]
    theta_ymax=[]
    
    #...gradient changes
    dtheta_xmin=[]
    dtheta_ymin=[]
    dtheta_xmax=[]
    dtheta_ymax=[]
    
    #Calculate the thickness of the CT slices (hard-coded: t=1.25 mm)
    t=1.25
    
    #Extract x and y boundaries for each slice contour
    for polygon in poly_contours['polygons']:
        xmin.append(polygon.bounds[0])
        ymin.append(polygon.bounds[1])
        xmax.append(polygon.bounds[2])
        ymax.append(polygon.bounds[3])
    
    #Calculate the gradient for each x,y boundary of the contour
    for i in range(num_of_cnts-1):
        
                
        theta_xmin_temp=math.degrees(math.atan((xmin[i]-xmin[i+1])/t))
        theta_ymin_temp=math.degrees(math.atan((ymin[i]-ymin[i+1])/t))
        theta_xmax_temp=math.degrees(math.atan((xmax[i]-xmax[i+1])/t))
        theta_ymax_temp=math.degrees(math.atan((ymax[i]-ymax[i+1])/t))
                       
        theta_xmin.append(theta_xmin_temp)
        theta_ymin.append(theta_ymin_temp)
        theta_xmax.append(theta_xmax_temp)
        theta_ymax.append(theta_ymax_temp)
    
    #Now calculate the the slice-to-slice gradient change (2nd gradient...)
    for i in range(num_of_cnts-2):
        dtheta_xmin.append(theta_xmin[i]-theta_xmin[i+1])
        dtheta_ymin.append(theta_ymin[i]-theta_ymin[i+1])
        dtheta_xmax.append(theta_xmax[i]-theta_xmax[i+1])
        dtheta_ymax.append(theta_ymax[i]-theta_ymax[i+1])
    
    dtheta_xmin=np.asarray(dtheta_xmin)
    dtheta_ymin=np.asarray(dtheta_ymin)
    dtheta_xmax=np.asarray(dtheta_xmax)
    dtheta_ymax=np.asarray(dtheta_ymax)
    
    #Calculate the mean absolute gradient change for each boundary.
    
    dtheta_xmin_av=np.mean(abs(dtheta_xmin))
    dtheta_ymin_av=np.mean(abs(dtheta_ymin))
    dtheta_xmax_av=np.mean(abs(dtheta_xmax))
    dtheta_ymax_av=np.mean(abs(dtheta_ymax))
    
    #Finally calculate the average of the four boundary gradient changes.
    DGrad_av=(dtheta_xmin_av+dtheta_ymin_av+dtheta_xmax_av+dtheta_ymax_av)/4
    
    return DGrad_av
    
def ExtractVolume(poly_contours):
    
    #First calculate the thickness of the CT slices
    t=abs(poly_contours['z-slice'][0]-poly_contours['z-slice'][1])
    #Now estimate the total volume         
    vol=0
    for polygon in poly_contours['polygons']:
        vol+=t*polygon.area
    
    return vol

def NumSlicesContoured(polyA,polyB):
  #Calculate the thickness of the CT slices (in mm)
  t=abs(polyA['z-slice'][0]-polyA['z-slice'][1])
  
  #Find the minum and maximum z-slice in both contoura  
  minz=min(polyA['z-slice'])
  maxz=max(polyA['z-slice'])
 
  if minz>min(polyB['z-slice']):
      minz=min(polyB['z-slice'])
  if maxz<max(polyB['z-slice']):
      maxz=max(polyB['z-slice'])
  
  #Now calculate the number of total number of slices contoured (by any structure)
  #This will be used to evaluate how many slices we don't have any intersection.
  
  num_tot_slices=int((maxz-minz)/t)
  
  return num_tot_slices

def FindSameSlices(polyA,polyB):
    #FUNCTION FINDSAMESLICES: Returns the z-slice index of polygon B that is the same
    #with polygon A. This is needed in order to since the two structures might 
    #not have contours always on the same slice.
    
    indA=[]
    indB=[]
    
    for j in range(len(polyA['z-slice'])):
     for i in range(len(polyB['z-slice'])): 
        if polyA['z-slice'][j]==polyB['z-slice'][i]:
            indA.append(j)
            indB.append(i)
            break       
        
    num_sl=len(indA)
    
    return indA,indB,num_sl
  
  

def ExtractIntersectionPoly(polyA,polyB):
  indB=[]
  polyABinter={}
  
  
  polyABinter['polygons']=[]
  polyABinter['z-slice']=[]
    
  #First find the slices that are common
  indA,indB,num_slices=FindSameSlices(polyA,polyB)
  
  #Now find the intersection for all the slices that are in common and create a new polygon list. 
  for i in range(num_slices):
    pA=polyA['polygons'][indA[i]]
    pB=polyB['polygons'][indB[i]]
    z=polyA['z-slice'][indA[i]]
    i=pA.intersection(pB)
    
    polyABinter['polygons'].append(i)
    polyABinter['z-slice'].append(z)
           
   
  return polyABinter

def ExtractUnionPoly(polyA,polyB):
    indB=[]
    polyABunion={}
    polyABunion['polygons']=[]
    polyABunion['z-slice']=[]
 
    #First find the slices that are common
    indA,indB,num_slices=FindSameSlices(polyA,polyB)
        
    #Now find the intersection and union area for all the slices that agree. 
    for i in range(num_slices):
         pA=polyA['polygons'][indA[i]]
         pB=polyB['polygons'][indB[i]]
         z=polyA['z-slice'][indA[i]]
         u=pA.union(pB)
        
         polyABunion['polygons'].append(u)
         polyABunion['z-slice'].append(z)
  
    return polyABunion

def HausdorffMetrics(polyA,polyB):
  #FUNCTION HAUSDORFF: Calculates the Haussdorff distances for each contour slice
  # between two polygon lists. Returns the max,mean and std of all Hausdorf distances.
  #Hausdorff distance is defined as: maximum distance of the closest point between
  #two contours.    
  
  #Initialization
    indB=[]
    H=[]  
    
    #First find the slices that are common
    indA,indB,num_slices=FindSameSlices(polyA,polyB)
    
    #Now calculate the Hausdorff distance of the contours for each slice
    for i in range(num_slices):
         pA=polyA['polygons'][indA[i]]
         pB=polyB['polygons'][indB[i]]
         htemp=pA.hausdorff_distance(pB)
         H.append(htemp)
    
    #Now calculate the max, mean and std of Hausdorff distances. Return it as a list.
    Hdmax=np.max(H)
    Hdmean=np.mean(H)
    Hdstd=np.std(H)
        
    return [Hdmax,Hdmean,Hdstd]

def Jaccard(poly):
    #function Jaccard: returns the Jaccard index of a list of polygons ('poly').
    #
    
    poly_inter=0
    poly_union=0
    jac_ind=[]
    weights=[]
    sum_area=0
    Jaccard_wav=0
    Jaccard_av=0
    
    #first extract the intersection and union of each slice for the first two 
    #polygons in the list. Then calculate the Jaccard index as the ratio 
    #of (intersection area) / (union area).

    num_of_str=len(poly)
    poly_inter=ExtractIntersectionPoly(poly[0],poly[1])
    poly_union=ExtractUnionPoly(poly[0],poly[1])
        
    if num_of_str>2:
        for i in range(num_of_str-2):
            poly_inter=ExtractIntersectionPoly(poly_inter,poly[i+2])
            poly_union=ExtractUnionPoly(poly_union,poly[i+2])
        for i in range(len(poly_inter['polygons'])):
            inter_area=poly_inter['polygons'][i].area
            union_area=poly_union['polygons'][i].area
            jac_ind.append(inter_area/union_area)
            weights.append(union_area)
            sum_area+=union_area
    else:
        for i in range(len(poly_inter['polygons'])):
            inter_area=poly_inter['polygons'][i].area
            union_area=poly_union['polygons'][i].area
            jac_ind.append(inter_area/union_area)
            weights.append(union_area)
            sum_area+=union_area
    
    # Now calculate the weighted average Jaccard index. The weihgts are calculated 
    # as the ratio: contoured area per slice/total contoured area.
      
        
    weights=[w/sum_area for w in weights]
        
    for i in range(len(weights)):
        Jaccard_wav+=weights[i]*jac_ind[i]
        
    Jaccard_wav=Jaccard_wav/sum(weights)
    Jaccard_av=np.mean(jac_ind)
    
    return jac_ind,Jaccard_av,Jaccard_wav
    
def HausdorffDict(poly):
    #FUNCTION HAUSDORFF POLY(): 
    #Given a series of structures (polygon lists) contoured by a physician 
    #returns the max, mean and std Haussdorf distance for all combinations 
    #(order selection does not matter here). For example calculates Haussdorf 
    #for structure 1 with 2,3.. structure 2 with 3,4..etc etc
    #INPUT: poly is a list of Polygon contour objects as created by the CreatePolygon() method.
    #OUTPUT: Dictionary with derived hausdorff metrics (max,mean,std) for each comparison performed. 
    #Thus, the length of each key in the dictionary should be equal to the total number 
    # of combinations possible. 
    
    
    #Initialization
    num_of_str=len(poly)
    Haus_dict={}
    Haus_dict['max']=[]
    Haus_dict['mean']=[]
    Haus_dict['std']=[]
    
    
    # We first perform a loop over all structures (num_of_str). For each structure (i) we will
    #calculate the Hausdorff distance to the remaining structures in the list
    #num_of_str - (i+1)  (i+1 because first index starts with 0 at python).
    #Return then the metrics (max,mean,std) in a dictionary. 
    
    for i in range(num_of_str):
        
        for j in range(num_of_str-(i+1)):
            
            j+=i+1
            Haus=HausdorffMetrics(poly[i],poly[j])
            Haus_dict['max'].append(Haus[0])
            Haus_dict['mean'].append(Haus[1])
            Haus_dict['std'].append(Haus[2])
            
            
    return Haus_dict
        
def main():
    
    ########## INITIALIZATION ##############
    poly=[]
    DGrad_av=[]
    vol=[]
    
        
    ############## USER INPUT ##############  
    dcm_list=input('Provide the names of all dicom files to be analyzed (format dcm1,dcm2,dcm3,...) \n')
    str_name=input('Provide the name of the structures contoured in each dicom file (example: GTV,CTV etc) \n')
    dcm_list=dcm_list.split(',')
    num_of_str=len(dcm_list)
    
    ########## CONTOUR EXTRACTION & ANALYSIS ############
    
    for i in range(num_of_str):
      
      #First create a polygon object
      poly_contours=CreatePolygon(dcm_list[i],str_name)
      
      #Now merge any multi polygons in the same slice
      p,ind=MergeMultiContours(poly_contours)
      
      #Calculate the average gradient change of the contour 
      dg=ExtractBoundGradChange(p)
      
      #Calculate the structure volume
      v=ExtractVolume(p)
      
      #Append the polygon lists, gradient changes and volumes for each structure to a poly list
      poly.append(p)
      DGrad_av.append(dg)
      vol.append(v)
    
    ############# SIMILARITY METRICS ##################
    # Calculate the Jaccard similarity index per slice (jac_ind), the average Jaccard
    # and the weighted average Jaccard.
    jac_ind,Jaccard_av,Jaccard_wav=Jaccard(poly)
    
    #Calculate the Haussdorf dictionary which includes max, mean and std of the 
    # Haussdorf distance for all slices. 
    Haus_dict=HausdorffDict(poly)
    
    ############# USER OUTPUT ################
    
    print('----------CONTOUR ANALYSIS----------\n')
    for i in range(num_of_str):     
        print('Structure',str_name,'of dicom file',dcm_list[i],' has volume of: ', vol[i],' mm^3 \n and an average gradient change of: ', DGrad_av[i],' degrees \n')
        
    print('The CT slices (z - mm) contoured were: \n',poly_contours['z-slice']) 
    print('The Jaccard similarity indices per slice were: \n',jac_ind)
    print('The average and weighted average Jaccard indices were: \n', Jaccard_av,Jaccard_wav)
    print('The Hausdorff metrics were: ',Haus_dict)
    
    ########### FILE OUTPUT ##################
    
    
     
        
if __name__ == '__main__':
  main()

 
    
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    