#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 27 17:28:32 2021

@author: jakubpazio
"""

import numpy as np
import cv2 as cv
import astropy as astro
from astropy.io import fits

# test = fits.open('test.fits')
test = fits.open('new_image.fits')
im = test[0].data   

# im = cv.imread('test.fits')
# imgray = cv.cvtColor(im, cv.COLOR_BGR2GRAY)




ret, thresh = cv.threshold(im, 3443, 255, cv.THRESH_TOZERO)
# tresh = thresh.astype(np.uint8)
scaled = np.divide(thresh,2**8)
imscale = scaled.astype(np.uint8)
contours, hierarchy = cv.findContours(imscale, cv.RETR_TREE, cv.CHAIN_APPROX_NONE)

print(len(contours))

# contourslist = []
# for i in range(len(contours)):
#     print(i)
#     if len(contours[i]) > 9:
#         contourslist.append(contours[i])
        

# im = im.astype(np.uint8)

# ret, thresh = cv.threshold(im, 1, 255, cv.THRESH_TOZERO)
# tresh = thresh.astype(np.uint8)
# contours, hierarchy = cv.findContours(thresh, cv.RETR_TREE, cv.CHAIN_APPROX_NONE)

#%%
contourslist = []
for i in range(len(contours)):
    print(i)
    if len(contours[i]) > 7:
        contourslist.append(contours[i])
        
#%%
len(contourslist)

# cv.drawContours(thresh, contours, -1, (30000,255,0), 1)

# hdu = fits.PrimaryHDU(thresh)   
hdu = fits.PrimaryHDU(thresh) 
# hdu.writeto('contour2.fits')   

#%%
# test = fits.open('test.fits')
# im = test[0].data   

# cv.drawContours(im, contourslist, -1, (30000,255,0), 1)

# hdu2 = fits.PrimaryHDU(im) 
# hdu2.writeto('contour4.fits')  

#%%
galaxies = []
for cnt in contourslist:
    
  print(len(contourslist),len(galaxies))
  area = cv.contourArea(cnt)
  # if area < 2:
  #   continue
  x, y, w, h = cv.boundingRect(cnt)
  # 1. X-coordinate
  x_coord = x + w/2.0
  # 2. Y-coordinate
  y_coord = y + h/2.0
  # 3. brightness
  star_mask = np.zeros(im.shape,np.uint8)
  cv.drawContours(star_mask, [cnt], 0, 255, -1)
  mean_val = cv.mean(im, mask=star_mask)[0]
  min_val, max_val, min_loc, max_loc = cv.minMaxLoc(im, mask=star_mask)
  bright = (mean_val-3419) * area
  # 4. radius
  radius = np.sqrt(area/(2*np.pi))
  galaxies.append([x_coord,
                y_coord,
                bright,
                mean_val,
                max_val,
                radius])
  
  
    # %%
    
f = open('data_new.txt', 'wb')
f.write(galaxies)


# %%
import matplotlib.pyplot as plt
brightneslist=[]
m=[]
for i in range(len(galaxies)):
    if galaxies[i][2] > 0:
        brightneslist.append(galaxies[i][2])
        m.append(25.3 - 2.5*np.log10(galaxies[i][2]))
 

    
def length_of_list(list_of_numbers, number):
     x = [i for i in list_of_numbers if i < number]
     return len(x)
 
b=0.2
x=[]
y=[]
for i in range(110):
    print(i)
    x.append(5+i*b)
    y.append(length_of_list(m, 5+ i*b))





# %%



stderror = np.sqrt(y)
xer = 0.02
plt.errorbar(x, y, yerr=stderror,xerr=xer,fmt='o',color='darkgreen',capsize = 3,label = 'Experimental data N(<m)')

plt.plot(np.arange(5,23,0.1),0.09*np.exp(0.6*np.arange(5,23,0.1)),'-',color='red',label = 'log(N) = 0.6m + const')


plt.xlabel('Magnitude limit m')
plt.ylabel('Number of Galaxies N(m)')
plt.yscale("log") 
# plt.xscale("log")
plt.legend()
plt.show()


# %%
import csv

with open('Galaxy_catalogue.txt', 'w') as f:
    csv_writer = csv.writer(f)
    csv_writer.writerow({'1.x_coordinates', '2.y_coordinates', '3.brightness', '4.mean_brightness', '5.max_brightness', '6.radius'})
    csv_writer.writerows(galaxies)







