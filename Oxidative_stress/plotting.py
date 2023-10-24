# -*- coding: utf-8 -*-
"""
Created on Sat Oct  1 15:36:33 2022

@author: some4879
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from mycolorpy import colorlist as mcp
 
                            
def plotlen (dfallcells,numberoftrench,treated,duration):
        colors=   mcp.gen_color(cmap="nipy_spectral",n=numberoftrench+1)
        image,ax =plt.subplots(figsize=(3,2), dpi=300)
        for concvalues in np.unique(dfallcells.concentration.values):
                dfx = dfallcells[dfallcells.concentration.values == concvalues]
                
                for trench in range(   numberoftrench) [0:1]:
                    
                    df2=  dfx[dfx.TrenchNumber== trench] 
                    
                    for postrench in range(0,1):#int(max(df.max_cells))):
                            df_mother = df2[df2.PositionClosedEnd ==  postrench]
                            df_mother=df_mother 
                            df_mothergr =  df_mother.groupby('Frame').mean()
                          
                            t=df_mothergr.index 
                  
                            sx2= df_mothergr.length.values
                            ax.plot(t-treated, sx2 , lw=1,  label='', color=colors[trench])
 
                            ax.set_xlabel('Time', fontsize=15)
                            ax.set_ylabel('Length', fontsize=15 )                            
def plotgrxa(dfallcells,numberoftrench,treated,duration):
    
        colors=   mcp.gen_color(cmap="nipy_spectral",n=numberoftrench+1)
       
        for concvalues in np.unique(dfallcells.concentration.values):
                dfx = dfallcells[dfallcells.concentration.values == concvalues]
                
                for trench in range(   numberoftrench): 
                    
                    df2=  dfx[dfx.TrenchNumber== trench] 
                    if len(df2)>0:
                        image,ax =plt.subplots(figsize=(2,2), dpi=300)
                        for postrench in range(0,1):#int(max(df.max_cells))):
                                df_mother = df2[df2.PositionClosedEnd ==  postrench]
                                if treated >300:
                                    df_mother=df_mother[df_mother.Time>treated -200]
                                else:
                                    df_mother=df_mother
                                df_mothergr =  df_mother.groupby('Frame').mean()
                              
                                t=df_mothergr.index 
                      
                                sx2= df_mothergr.GrxA.values
                                ax.plot(t-treated, sx2 , lw=1,  label='GrxA' , color=colors[postrench])
                         
                                ax.set_xlabel('Time', fontsize=15)
                                ax.set_ylabel('GrxA', fontsize=15 )  
                                ax.set_title('Cell number = %s'%trench)



    
def plotpoincare(dfallcells,numberoftrench,treated,duration):
    
        colors=   mcp.gen_color(cmap="nipy_spectral",n=numberoftrench+1)
        #image,ax =plt.subplots(figsize=(2,2), dpi=300)
        for concvalues in np.unique(dfallcells.concentration.values):
                dfx = dfallcells[dfallcells.concentration.values == concvalues]
                
                for trench in range(   numberoftrench):
                    
                    df2=  dfx[dfx.TrenchNumber== trench]
                    
                    if len(df2)>0:
                        image,ax =plt.subplots(figsize=(2,2), dpi=300)
                        for postrench in range(0,1):
                                df_mother = df2[df2.PositionClosedEnd ==  postrench]
                                df_mothergr =  df_mother.groupby('Frame').mean()
                              
                                t=df_mothergr.index 
                      
                                sx2= df_mothergr.GrxA.values
                                sx2= sx2[treated+500:]
                                ax.plot(sx2[10:], sx2[:-10] , lw=1,  label='GrxA' , color=colors[postrench])
                         
                                ax.set_xlabel('GrxA at t+10', fontsize=15)
                                ax.set_ylabel('GrxA', fontsize=15 ) 
                                ax.set_title('Cell number = %s'%trench)

def plotgrxagrad(dfallcells,numberoftrench,treated,duration):
   colors=   mcp.gen_color(cmap="nipy_spectral",n=numberoftrench+1)
   image,ax =plt.subplots(figsize=(2,2), dpi=300)
   for concvalues in np.unique(dfallcells.concentration.values)[0:1]:
               df2 = dfallcells[dfallcells.concentration.values == concvalues]
               for postrench in range(1,6):
                       df_mother = df2[df2.PositionClosedEnd ==  df2.max_cells - postrench]
                       df_mother=df_mother
                       df_mothergr =  df_mother.groupby('Frame').mean()
                     
                       t=df_mothergr.index 
                       sx2= df_mothergr.GrxA.values

                       ax.plot(t[5:]-treated, sx2[5:] , lw=1, color = colors[postrench] )
                       ax.set_xlabel('Time', fontsize=15)
                       ax.set_ylabel('GrxA', fontsize=15 )
