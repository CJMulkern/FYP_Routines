# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 12:25:22 2025

@author: Cjmul
"""
import pandas as pd
import numpy as np
#import matplotlib.pyplot as plt
from astropy.io import fits 
from astropy.modeling import models, fitting
#%%
def Gaussian_func(path, x_0, y_0, x_range, y_range):
    
    image = fits.open(path)
    
    hdu_list_new = image
    #hdu_list_new.info()
    image_data = hdu_list_new[0].data
    
    hdu_list_new.close()

    ps = image_data[y_0 - y_range:y_0 + y_range,x_0 - x_range: x_0 + x_range]

    #plt.figure()
    #ps_zoom = plt.imshow(ps,origin ='lower',  cmap='hot', vmin = np.min(image_data)/255, vmax =np.max(image_data)/255)
    #plt.colorbar()
    #plt.show(ps_zoom)
    maxval = np.array(ps).ravel()[np.argmax(ps)]

    coords = np.unravel_index(np.argmax(ps), 
                              shape = np.array(ps).shape)
    fitter = fitting.LevMarLSQFitter()

    y_shape = np.shape(ps)[1]
    x_shape = np.shape(ps)[0]
    y,x = np.mgrid[:y_shape,:x_shape]
    z = ps
    p_init = models.Gaussian2D(amplitude =maxval, x_mean=coords[1], y_mean=coords[0], x_stddev=4, y_stddev=4, theta=None)
    p = fitter(p_init, x, y, z)
    
    x_prime = p.x_mean[0]
    y_prime = p.y_mean[0]
    return x_prime, y_prime

#%%
file = pd.read_excel('C:/Users/Cjmul/OneDrive/Desktop/College, Baby/Final Year/Final Year Project/Python_Files/Pandas_path.xlsx')

for i in range(13):
    # row1 = file.iloc[i]
    row1 = file.iloc[183+i]
    #print(row1)
    path = row1.Path
    x_0 = row1.x_0
    y_0 = row1.y_0
    x_range = row1.x_range
    y_range = row1.y_range
    


    A = Gaussian_func(path, x_0, y_0, x_range, y_range)
    
    x_f = x_0 - x_range + A[0]
    y_f = y_0 - y_range + A[1]
    
    #file.loc[i,'x_f'] = x_f
    #file.loc[i,'y_f'] = y_f
    file.loc[183+i,'x_f'] = x_f
    file.loc[183+i,'y_f'] = y_f
    
    
    print( "<x= {0:.3f}, y={1:.3f}>".format(x_f,y_f))
    
file.to_excel('C:/Users/Cjmul/OneDrive/Desktop/College, Baby/Final Year/Final Year Project/Python_Files/Pandas_path.xlsx', index=False)

#%%
file = pd.read_excel('C:/Users/Cjmul/OneDrive/Desktop/College, Baby/Final Year/Final Year Project/Python_Files/Pandas_path.xlsx')
row1 = file.iloc[195]
print(row1)
#%%
path = row1.Path
x_0 = int(row1.x_0)
y_0 = int(row1.y_0)
x_range = int(row1.x_range)
y_range = int(row1.y_range)

diffy = y_0 - y_range
#%%

A = Gaussian_func(path, x_0, y_0, x_range, y_range)

#%%
x_f = x_0 - x_range + A[0]
y_f = y_0 - y_range + A[1]
#%%
#%%  
file = pd.read_excel("C:/Users/Cjmul/OneDrive/Desktop/College, Baby/Final Year/Final Year Project/Python_Files/abc.xlsx")

for i in range(8):
     row1 = file.iloc[i]
     #print(row1)
     path = row1.Path
     x_0 = row1.x_0
     y_0 = row1.y_0
     x_range = row1.x_range
     y_range = row1.y_range
     


     A = Gaussian_func(path, x_0, y_0, x_range, y_range)
     
     x_f = x_0 - x_range + A[0]
     y_f = y_0 - y_range + A[1]
     
     file.loc[i,'x_f'] = x_f
     file.loc[i,'y_f'] = y_f
     
     print( "<x= {0:.3f}, y={1:.3f}>".format(x_f,y_f))
     
file.to_excel("C:/Users/Cjmul/OneDrive/Desktop/College, Baby/Final Year/Final Year Project/Python_Files/abc.xlsx", index=False)
