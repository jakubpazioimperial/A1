#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  2 15:13:41 2021

@author: jakubpazio
"""

import numpy as np
import matplotlib.pyplot as plt

import astropy as astro
from astropy.io import fits
from scipy.optimize import curve_fit

#%%

mosaic = fits.open('mosaic.fits')

image = mosaic[0].data                           #This is a 2D array so will add all the values to a 1d array
image_list = np.ndarray.tolist(image)

#%%
# image_list = []

# for i in range(len(image)):
    # for j in range(len(image[i])):
        # image_list.append(image[i][j])
    
#%%

# image_list = np.ravel(image)
#%%

# bin_width = 1

# ydata=plt.hist(image_list, bins = int(2**16/bin_width))
# # plt.hist(image_list, bins = 2**14)
# plt.xlabel('Pixel Brightness')
# plt.ylabel('Occurence')

# xdata=np.arange(0,bin_width*len(ydata[0]),bin_width)

# # plt.plot(xdata+bin_width/2, ydata[0],'-',color='red')
# plt.show()


# %%

def gaussian(x, amp, cen, wid):
    return amp * np.exp(-(x-cen)**2 / wid)
# left=int(3420/bin_width-100)
# right=int(3420/bin_width+100)
# bin_range = xdata[left:right]
# y_range = ydata[0][left:right]
# initialGuess = [10**6,3420,100]  

# popt, pcov = curve_fit(gaussian, bin_range+bin_width/2 , y_range, initialGuess)
# # popt, pcov = curve_fit(gaussian, bin_range , Am_241_range, initialGuess)
# print(popt)
# print(pcov)



# param=[0,0]
# param[0]=popt[1]
# param[1]=np.sqrt(popt[2]/2)

# def PlotGaussians():

#     plt.plot(xdata,gaussian(xdata,*popt),color = 'red',label='fit params: cen=%5.3f , wid=%5.3f' % tuple(param))
    
    
#     plt.title('Fitted Gaussians to the Brightness Histogram',weight='bold')
#     plt.xlabel('Pixel Brightness', weight = 'semibold')
#     plt.ylabel('Occurence', weight = 'semibold')
#     plt.legend(title='Scattering Angle',title_fontsize=15,fontsize='small')
#     plt.show()
    
    
    
# PlotGaussians()
    
    
# %%
'''
This section deletes the spike in the middle of the image 
'''

    
for j in range(4611):
    for i in range (23):
        if image[j][1425+i] > 5200:
            
            image[j][1425+i] = 3320
    
        
        
# hdu = fits.PrimaryHDU(image)        
        
# hdu.writeto('new_image.fits')       
        

# %%

'''
This section deletes noise 
'''
#53 87
#10 257
bin_width=1 
for j in range(87): #87
    for i in range (10): #10
        image_list=[]
        for k in range(53):     
            image_list=np.append(image_list, np.ravel(image[53*j+k][i*257:(i+1)*257]))
        binnumber=max(image_list)-min(image_list)+100
        print('j',j,'i',i)
        print(binnumber)
        
        
        ydata=np.histogram(image_list, bins = int(binnumber))
        # plt.hist(image_list, bins = 2**14)
        # plt.xlabel('Pixel Brightness')
        # plt.ylabel('Occurence')
        
        xdata=ydata[1]
        
        # left=int(3420/bin_width-100)
        # right=int(3420/bin_width+100)
        # bin_range = xdata[left:right]
        # y_range = ydata[0][left:right]
        initialGuess = [10**3,3421,1]  
        
        popt, pcov = curve_fit(gaussian, xdata[:-1] , ydata[0], initialGuess)
        print(popt)
        print(pcov)
        


        noise_brightnes = popt[1]+np.sqrt(popt[2]/2) # value must be tahen from histogram
        
        
        if image[53*j][253*i] < noise_brightnes: # noise_brightnes - variable derived from gaussian
            
            image[j][i] = 0
    
        
        
     
        
# %%
hdu = fits.PrimaryHDU(image)   
hdu.writeto('new_image.fits')   


# %%
image_list=np.append(image_list, np.ravel(image[2345][0:257]))


print(image_list)


mosaic.close()

# mosaic[0].header
# MAGZPT  =            2.530E+01 / Photometric ZP (mags) for default extinction   

    
    
    
