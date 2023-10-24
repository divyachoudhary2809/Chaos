# -*- coding: utf-8 -*-
"""
Created on Thu Mar  9 10:51:41 2023

@author: some4879
"""
import main_file as run 
import random

conx=10.
 
duration = 1200 ; trenchLength =25. ;   ### trench length in um
g0 = 0.042;  concentration=conx ### treatment concentration in uM  
time_of_treat= 800

numberoftrench = 1
trenches=[[]]* numberoftrench
for ntrench in range(numberoftrench):
    trenches[ntrench] =[] ;  
    len_=0
    for i in range(0, random.randint(3,8)):
        if len_ < trenchLength:
            lenadd= random.randint(2010, 3990)/1000.
            trenches[ntrench] .append(lenadd )
            len_+=lenadd#'''

df = run.main_run( trenchLength, g0, trenches, duration, concentration,time_of_treat) 
 