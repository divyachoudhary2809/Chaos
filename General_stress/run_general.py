# -*- coding: utf-8 -*-
"""
Created on Thu Mar  9 10:30:39 2023

@author: some4879
"""

import main_file_general as run
import random

Kcat= 20.  ##Kcat
Kactivation =0.1 ## Kactivation
concentration = 150. ## toxin concentration
duration = 2000# ## duration of simulation in minutes
TrenchLength = 25.### trench length in um
treated=50


### if you want to choose randomised initial conditions
## for having a randomised structure of 'N' = number of trenches.  just supply N in the line below
numberoftrench = 3
trenches=[[]]* numberoftrench
for ntrench in range(numberoftrench):
    trenches[ntrench] =[] ;  
    len_=0
    for i in range(0, random.randint(3,8)):
        if len_ < TrenchLength:
            lenadd= random.randint(2010, 3990)/1000.
            trenches[ntrench] .append(lenadd )
            len_+=lenadd#'''
            
            
## for having initial structure of your choice you can use uncomment the line below and add the values you like: 
#trenches=[[2.241, 2.342, 3.043, 2.244,3.067,3.241]]  ### structure of cells at the start
        
        
run.main_run( Kcat, Kactivation, concentration, duration, trenches, TrenchLength,treated)
