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
    # Function ReadContour: accepts a dicom (.dcm) file and returns a dictionary
    # with all the existing contours in the plan.
    
    ds=dicom.read_file(DCM) # first read .dcm file
    contours=[]
    num_of_contours=len(ds.ROIContourSequence)
    for i in range(num_of_contours):
        contour={}
        contour['color']=ds.ROIContourSequence[i].ROIDisplayColor
        #contour['number']=ds.ROIContourSequence[i].RefdROINumber
        contour['name']=ds.StructureSetROISequence[i].ROIName
        contour['contour data']=[s.ContourData for s in ds.ROIContourSequence[i].ContourSequence]
        contours.append(contour)
    return contours

    
def FindContour(contours,roi_name):
   #Function ReadContourPoints: accepts the dicom read structure and contour names
   # (output for ReadContourNames and returns an array with all the contour 
   #points for each slice z.
   
      
   # First the user asks for a contour and we verify that the name exists in 
   # the contour list
   
   name_agrees=[roi_name==contour['name'] for contour in contours]
   if sum(name_agrees)==1:
        contour_ind=name_agrees.index(1)
        print('The contour', roi_name,' has been selected\n')
        return contour_ind
   else:
        print('Name not found. Restart and provide a structure name from the list below\n') 
        print([contour['name'] for contour in contours])
        return [contour['name'] for contour in contours]
    

def ExtractContourPoints(contours,contour_ind):
       #initialization 
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
        #initialize lists and dictionaries
        poly_list=[]
        poly_z=[]
        poly_contours={}

        num_contours=len(coord['x'])
        
        # for each contour slice (z) read the x,y coordinates and create a polygon 
        # oblject.
        for i in range(num_contours):
            num_crds=len(coord['z'][i])
            
            poly_crds=[]
            for j in range(num_crds):
                poly_crds.append((coord['x'][i][j],coord['y'][i][j]))
            
            poly_z.append(coord['z'][i][j])
            poly_contour=Polygon(poly_crds)
            poly_list.append(poly_contour)
            
        #Create a dictionary with all the polygon objects with their 
        #respective z-position
        poly_contours['polygons']=poly_list
        poly_contours['z-slice']=poly_z
        
        return poly_contours
    
def CreatePolygon(dcm,roi_name):
    
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
    #Merges multi contours that exist in the same slice to one unified contour.
    #This serves also as a clean-up method fot he existance of small contour dots
    
    z=poly_contours['z-slice']
    poly=poly_contours['polygons']
    ind=[]
    
    #First find the index for the z slices that equal uniify the polygons. 
    #This returns a multi-polygon structure
    
    for i in range(len(z)-1):
      if z[i]==z[i+1]:
        poly[i+1]=poly[i+1].union(poly[i])
        ind.append(i)
        
    
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
    
    #Calculate the thickness of the CT slices
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

def ExtractIntersectionPoly(polyA,polyB):
  indB=[]
  num_slices_found=0
  polyABinter={}
  
  
  polyABinter['polygons']=[]
  polyABinter['z-slice']=[]
    
    
  #First find the number of total slices contoured by any structure.
  #Now find the slices of polyB that have the same z coordinate as polyA. 
  #These are the slices that we find contours of both structures.
  
  for z in polyA['z-slice']:
    found_slice=0
    for i in range(len(polyB['z-slice'])): 
        if z==polyB['z-slice'][i]:
            found_slice=1
            break
    num_slices_found+=found_slice               
    indB.append(i)
  
  
  #Now find the intersection for all the slices that agree and create a new polygon list. 
  for i in range(len(polyA['polygons'])):
    pA=polyA['polygons'][i]
    pB=polyB['polygons'][indB[i]]
    z=polyB['z-slice'][indB[i]]
    i=pA.intersection(pB)
    
    polyABinter['polygons'].append(i)
    polyABinter['z-slice'].append(z)
           
   
  return polyABinter

def ExtractUnionPoly(polyA,polyB):
    indB=[]
    polyABunion={}
    polyABunion['polygons']=[]
    polyABunion['z-slice']=[]
 
  #Now find the slices of polyB that have the same z coordinate as polyA. 
  #These are the slices that we find contours of both structures.
  
    for z in polyA['z-slice']:
     for i in range(len(polyB['z-slice'])): 
        if z==polyB['z-slice'][i]:
            break
     indB.append(i)
    
    
    #Now find the intersection and union area for all the slices that agree. 
    for i in range(len(polyA['polygons'])):
         pA=polyA['polygons'][i]
         pB=polyB['polygons'][indB[i]]
         z=polyB['z-slice'][indB[i]]
         u=pA.union(pB)
        
         polyABunion['polygons'].append(u)
         polyABunion['z-slice'].append(z)
  
    return polyABunion

def Jaccard(poly):
    #function Jaccard: returns the Jaccard index of a list of polygons ('poly').
    #
    
    poly_inter=0
    poly_union=0
    jac_ind=[]
    
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
    else:
        for i in range(len(poly_inter['polygons'])):
            inter_area=poly_inter['polygons'][i].area
            union_area=poly_union['polygons'][i].area
            jac_ind.append(inter_area/union_area)
    Jaccard_av=np.mean(jac_ind)
    
    return jac_ind,Jaccard_av

def main():
    
    poly=[]
    DGrad_av=[]
    vol=[]
    
        
    ############## User input & communication ##################  
    dcm_list=input('Provide the names of all dicom files to be analyzed (format dcm1,dcm2,dcm3,...) \n')
    #roi_list=input('Provide the names of the structures contoured in each dicom file: roi1,roi2,roi3,... \n')
    dcm_list=dcm_list.split(',')
    #roi_list=roi_list.split(',')
    roi_list=['Target','Target']
    num_of_str=len(roi_list)
    
    for i in range(num_of_str):
      
      ############ Analyze contour ###############
      #First create a polygon object
      poly_contours=CreatePolygon(dcm_list[i],roi_list[i])
      
      #Now merge any multi polygons in the same slice
      p,ind=MergeMultiContours(poly_contours)
      #Calculate the average gradient change of the contour 
      dg=ExtractBoundGradChange(p)
      #Calculate the structure volume
      v=ExtractVolume(p)
      
      poly.append(p)
      DGrad_av.append(dg)
      vol.append(v)
    
    jac_ind,Jaccard_av=Jaccard(poly)
    
    for i in range(num_of_str):     
        print('----------CONTOUR ANALYSIS----------\n')
        print('Structure',roi_list[i],'of dicom file',dcm_list[i],' has volume of: ', vol[i],' mm^3 \n and an average gradient change of: ', DGrad_av[i],' degrees \n')
        
        
    print('The average Jaccard index of all slices is equal to: \n', Jaccard_av)
    print('The Jaccard index per slice is: \n',jac_ind)
    print('The z-slice index is: \n',poly_contours['z-slice'])   
     
        
if __name__ == '__main__':
  main()

 
    
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    