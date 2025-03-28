# -*- coding: utf-8 -*-
"""
Created on Fri Mar  7 09:32:20 2025

@author: Cjmul

"""
import pandas as pd
import photutils.aperture as pu
from astropy.io import fits 
import numpy as np
#%%#

def App_Photom(path, x_ps, y_ps):
    # Open the FITS file
    image = fits.open(path)
    image_data = image[0].data
    image.close()

    # Define the original aperture
    aper = pu.CircularAperture([x_ps, y_ps], r=8)
    flux = np.array(pu.aperture_photometry(image_data, aper)['aperture_sum'])

    # Convert the coordinates to polar (r, θ)
    x_rel = x_ps - 511.5
    y_rel = y_ps - 511.5
    r = np.sqrt(x_rel**2 + y_rel**2)  # Radial distance
    theta = np.arctan2(y_rel, x_rel)  # Angle (in radians)

    # Initial antipodal coordinates (same r, same θ)
    x_a = r * np.cos(theta) + 511.5
    y_a = r * np.sin(theta) + 511.5

    # Perform aperture photometry at the mirrored position
    aper_a = pu.CircularAperture([x_a, y_a], r=8)
    flux_a = np.array(pu.aperture_photometry(image_data, aper_a)['aperture_sum'])
    
    
    # If the flux at the opposite point is more than half the original flux, rotate the coordinates by 90 degrees
    for i in range(1000):
        actual = flux - flux_a
        if flux_a <= 0.8 * abs(flux) and flux_a > 0:
            break
        else:
                # Rotate the angle by some angle
                random_rotation = np.random.uniform(0, 2 * np.pi)
                
                theta = theta + random_rotation
                
                # Update the antipodal coordinates with the rotated angle
                x_a = r * np.cos(theta) + 511.5
                y_a = r * np.sin(theta) + 511.5
                
                # Perform aperture photometry again at the new position
                aper_a = pu.CircularAperture([x_a, y_a], r=8)
                flux_a = np.array(pu.aperture_photometry(image_data, aper_a)['aperture_sum'])
                
                # Skipping iteration beacause flux_a is too large
                if flux_a >= flux:
                    continue
                # Trying to ensure flux difference is +ve
                if actual <0:
                    continue

    # Calculate the difference in fluxes
    actual = flux - flux_a

    return flux, flux_a, actual, x_a, y_a

#%%

file = pd.read_excel("C:/Users/Cjmul/OneDrive/Desktop/College, Baby/Final Year/Final Year Project/Python_Files/Apperture_path.xlsx")

for i in range(13):
    row1=file.iloc[183+i]
    path = row1.Path
    x_ps = row1.x_ps
    y_ps = row1.y_ps
    
    A = App_Photom(path, x_ps, y_ps)
    
    Source_flux = A[0]
    Antipodal_flux = A[1]
    Actual_flux = A[2]
    x_a = A[3]
    y_a = A[4]
    
    file.loc[183+i,'Source_flux'] = Source_flux
    file.loc[183 + i,'Antipodal_flux'] = Antipodal_flux
    file.loc[183 + i,'Actual_flux'] = Actual_flux
    #file.loc[i,'x_a'] = x_a
    #file.loc[i,'y_a'] = y_a
    
    print("<flux val = {0:} >".format(Actual_flux))
    
file.to_excel("C:/Users/Cjmul/OneDrive/Desktop/College, Baby/Final Year/Final Year Project/Python_Files/Apperture_path.xlsx",index =False)
#%%
file = pd.read_excel("C:/Users/Cjmul/OneDrive/Desktop/College, Baby/Final Year/Final Year Project/Python_Files/Apperture_path.xlsx")
row1=file.iloc[183]
print(row1)
#%%


path = 'C:/Users/Cjmul/OneDrive/Desktop/College, Baby/Final Year/Final Year Project/Python_Files/data/output(180).fits'
x_ps = 514.2114974
y_ps = 481.2713808
    
A = App_Photom(path, x_ps, y_ps)
    
Source_flux = A[0]
Antipodal_flux = A[1]
Actual_flux = A[2] 
x_a = A[3]
y_a = A[4]

print(A)
#%%
path = "C:/Users/Cjmul/OneDrive/Desktop/College, Baby/Final Year/Final Year Project/Python_Files/star_removed.fits"
    
A = App_Photom(path, x_ps, y_ps)
    
Source_flux = A[0]
Antipodal_flux = A[1]
Actual_flux = A[2] 
x_a = A[3]
y_a = A[4]

print(A)
#%%
path = 'C:/Users/Cjmul/OneDrive - National University of Ireland, Galway/Cian_FYP_Data/Gaia DR2 5854897321965963264/Gaia_DR2_5854897321965963264_2022-02-14_I_tot.fits'
    
A = App_Photom(path, x_ps, y_ps)
    
Source_flux = A[0]
Antipodal_flux = A[1]
Actual_flux = A[2] 
x_a = A[3]
y_a = A[4]

print(A)