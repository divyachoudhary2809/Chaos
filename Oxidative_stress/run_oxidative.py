# -*- coding: utf-8 -*-
"""
Created on Thu Mar  9 10:51:41 2023

@author: some4879
"""
import main_file as run 
import random

 
duration = 6000 ; #total duration
time_of_treat= 800
trenchLength =25. ;   ### trench length in um
g0 = 0.042;  
concentration=260. ### treatment concentration in uM  


### if you want to choose randomised initial conditions
## for having a randomised structure of 'N' = number of trenches.  just supply N in the line below
numberoftrench = 3
trenches=[[]]* numberoftrench
for ntrench in range(numberoftrench):
    trenches[ntrench] =[] ;  
    len_=0
    for i in range(0, random.randint(3,8)):
        if len_ < trenchLength:
            lenadd= random.randint(2010, 3990)/1000.
            trenches[ntrench] .append(lenadd )
            len_+=lenadd#'''
            
            
## for having initial structure of your choice you can use uncomment the line below and add the values you like: 
#trenches=[[2.241, 2.342, 3.043, 2.244,3.067,3.241]]  ### structure of cells at the start
        
             
        
df = run.main_run( trenchLength, g0, trenches, duration, concentration,time_of_treat) 
 
