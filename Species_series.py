# -*- coding: utf-8 -*-
"""
Created on Fri Mar 14 17:20:41 2025

@author: Cjmul
"""

import matplotlib.pyplot as plt
import pandas as pd
from species import SpeciesInit
from species.data.database import Database
from species.read.read_isochrone import ReadIsochrone

#%%
SpeciesInit()
database = Database()
database.add_isochrones(model='ames')
read_iso = ReadIsochrone(tag='ames-dusty')
#%%
def est_mass( dist, star, del_mag):
    
    mass = read_iso.contrast_to_mass(age=5.,distance=dist,
                                       filter_name='H',
                                       star_mag=star,
                                       contrast=[del_mag],
                                       use_mag=True,atmospheric_model = 'ames-dusty')
    mass_plus = read_iso.contrast_to_mass(age=8.,distance=dist + dist_err,
                                       filter_name='H',
                                       star_mag=star,
                                       contrast=[del_mag],
                                       use_mag=True,atmospheric_model = 'ames-dusty')
    mass_err = mass - mass_plus
    return mass, mass_err
#%%

file = pd.read_excel("C:/Users/Cjmul/OneDrive/Desktop/College, Baby/Final Year/Final Year Project/Python_Files/Species_path.xlsx")

for i in range(190):
    row = file.iloc[i]
    dist = row.Dist
    star = row.Mag_star
    del_mag = row.Del_mag    
    dist_err = row.Dist_Uncertainty
    
    A = est_mass(dist, star, del_mag)
    
    Mass = A[0]
    Error = A[1]
    
    file.loc[i,'Mass_Mj'] = A[0]
    file.loc[i,'Mass_error'] = A[1]
    
    print("<mass is: {0:.3f}>".format(float(A[0])))
    
file.to_excel("C:/Users/Cjmul/OneDrive/Desktop/College, Baby/Final Year/Final Year Project/Python_Files/Species_path.xlsx")

#%%
'''
file = pd.read_excel("C:/Users/Cjmul/OneDrive/Desktop/College, Baby/Final Year/Final Year Project/Python_Files/Species_path.xlsx")
row = file.iloc[189]
print(row)
'''
#%%
'''
row = file.iloc[189]
dist = row.Dist
star = row.Mag_star
del_mag = row.Del_mag    
dist_err = row.Dist_Uncertainty
    
A = est_mass(dist, star, del_mag)
    
Mass = A[0]
Error = A[1]

print(A)
'''