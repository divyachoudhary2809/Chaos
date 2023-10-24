import pandas as pd
import math 
import numpy as np
from scipy.integrate import odeint
import imagesFromModel_function_1001 
import abstract_response_model
import plotting 
import random
import matplotlib.pyplot as plt





def main_run (    kcats, Kenzyme , conx,duration  , trenches_,limm,treated):
    modelfunc= imagesFromModel_function_1001.main()
    um =   (10**-6); 
    kgi= .02/um ; cI= 100.*um
     
    kgic= kgi ; cIc=cI  
    
    rc = 1.15/2.  ;w= 1.2
    
    bH = 0.00000 ; bKatG= 0.00; 
    um =   (10**-6); 
    #treated= 50
    ah2= 1.  
                      
    def tps(ll,timedur , trenches,yggput, c0init,kcat,KKatGx) :
              t=[[] for i in range(len(trenches))]
              cellnn = [[ 0 for i in range(len(trenches))]for i in range(len(timedur)-2)]
              cellnumberarray,alllengths,alllgrowths,  h2o2_ , katgall, ygg_cells,conc_cells,cellidx, idxchange, beforeElement,extraa,conc_t_lx,cellbirth= [ [[ [] for i in range(len(trenches))]for i in range(len(timedur)-2)] for i in range(13)]
                
              for lx in range(len(timedur)-2):
                  l =  timedur[lx]
                  celln =[0. for i in range(len(trenches))]  
                  totlengths, totlengths2, totgrowths2,   h2o2_2,   katgall2,  cellnn,t1 = [[[] for i in range(len(trenches))] for i in range(7)]
                
                  for trenchNum in range(len(trenches)):      
                         lenofcells,t1[trenchNum],totlengths[trenchNum] ,celln[trenchNum] ,lengthss =  [ ll[trenchNum], [],0.,0,0.00]
                         ygg_cells[lx][trenchNum]=[yggput[-1]]* len(lenofcells)
                         
                         if lx==0:
                             cellbirth[0][trenchNum]=[20.]* len(lenofcells)
                         if lx>0:
                             ygg_cells[lx][trenchNum]=[[]]* len(lenofcells)
                            
        
                         nfr=2
                         tp, timee ,prevtlen = [np.array([timedur[lx2]/1. for lx2 in range(lx,lx+2)]), np.linspace(lx,lx+1,nfr) , len(alllengths[lx-1][trenchNum] )]
        
                         for cellprev in range(prevtlen-1, -1,-1): 
                            vcell= alllengths[lx-1][trenchNum][cellprev]  /10.  
                            diff=1.*10**-9  
                             
                            diff= diff*1.
                            lambda_ = np.sqrt((((w*w) -(3.14* rc*rc))  *diff * 10**-12)/ ((2*3.14*rc*(10**-6)*1.6* (10**-5))   ))  
                            lenT, rcd, wd = [vcell *um,rc*um,w*um ]
                            gams= 1+  ((2.1* (rcd**3))/((wd*wd*lenT) - (3.14*rcd*rcd*lenT)))
                          
                            lambda_=lambda_ * math.sqrt(gams)     ; alpha_ , beta = [ lenT /lambda_ , 0.]
                          
                            if cellprev == prevtlen-1:
                                c0=c0init[lx]
                                c = c0 * ((np.e**-beta) +(np.e** beta)) / ((np.e**-alpha_) +(np.e** alpha_)) 
                                c0=c
                                conc_t_lx[lx-1][trenchNum]+=[c ]
                            else:
                                c = c0 * ((np.e**-beta) +(np.e** beta)) / ((np.e**-alpha_) +(np.e** alpha_)) 
                                c0=c
                                conc_t_lx[lx-1][trenchNum]+=[c ]
                                
                                
                                
                         for cell in range(len(lenofcells)):#-1, -1,-1): 
                                cellx = len(lenofcells)-cell-1  
                                if lx>0:
                                         c=conc_t_lx[lx-1][trenchNum][ len(lenofcells)-cell-1]  
                                         ygg_cells[lx][trenchNum][cell] = odeint(abstract_response_model.comb_model3_growth,ygg_cells[lx-1][trenchNum][cell] ,timee,args=(tp,[c,c], gI , kgi, cI,kcat,KKatGx ))
                                else:
                                        c = c0init[lx]
                                       
                                        ygg_cells[lx][trenchNum][cell] = odeint(abstract_response_model.comb_model3_growth,ygg_cells[lx][trenchNum][cell],timee,args=(tp,[c,c] ,gI , kgi, cI,kcat,KKatGx))
         
                                inix= ygg_cells[lx][trenchNum][cell] [-1]
                                
                                g = np.mean(gI-(  gI / (1 + 10**(-kgic*((ygg_cells[lx][trenchNum][cell][:,0]*1.)-cIc)))   ))
                                
                                ygg_cells [lx][trenchNum][cell]=inix
                                if lx>0:
                                    cellbirth[lx][trenchNum]=cellbirth[lx-1][trenchNum] 
                               
                                xmega= [h2o2_2[trenchNum], katgall2[trenchNum] ]
                                for i in range(len(xmega)): xmega[i] +=[inix[i]]
                                h2o2_2[trenchNum], katgall2[trenchNum] = xmega
                                totgrowths2[trenchNum]  += [g]
                               
                                lenadded =  (lenofcells[cell]* g); lenofcells[cell] +=   lenadded  ;  length = lenofcells[cell]; birth = cellbirth[lx][trenchNum] [cell]
                             
                                if length -birth <=20.:  
                                    length1 = length
                                    if   totlengths[trenchNum] < limm:
                                        lengthss += (length1); totlengths2[trenchNum]  += [length1]; celln[trenchNum]+=1; cellnn[trenchNum]=[celln[trenchNum]]
                                        totlengths[trenchNum]  = lengthss; t1[trenchNum]+=[totlengths[trenchNum], length1]
        
                                else:        
                                    length1 = length*0.5  
                                    cellbirth[lx][trenchNum][cell] = length1
                                    
                                    beforeElement[lx][trenchNum] += [cell+1 ]; extraa[lx][trenchNum]+=[cell];   cellidx[lx][trenchNum] +=[1]; idxchange[lx][trenchNum] +=[1]
                                                               
                                    for k in range(0,2):   
                                            
                                            if totlengths[trenchNum] < limm:
                                             
                                                lengthss +=  length1;  celln[trenchNum]+=1; cellnn[trenchNum]=[celln[trenchNum]]
                                                totlengths2[trenchNum]  += [length1]; totlengths[trenchNum]  = lengthss
                                                t1[trenchNum]+=[totlengths[trenchNum], length1]   
                                                if k==0:
                                                    if lx>0:
                                                        xmega= [h2o2_2[trenchNum], katgall2[trenchNum] ]
                                                        for i in range(len(xmega)): xmega[i] +=[inix[i]]
                                                        h2o2_2[trenchNum] ,katgall2[trenchNum]   = xmega
                                                        totgrowths2[trenchNum]  += [g]
                                                    else:
                                                        xmega= [h2o2_2[trenchNum] ,katgall2[trenchNum] ]
                                                        for i in range(len(xmega)): xmega[i] +=[inix[i]]
                                                        h2o2_2[trenchNum] ,katgall2[trenchNum] = xmega
                                                        totgrowths2[trenchNum]  += [g]         
                         if cellx==0    and len(cellidx[lx][trenchNum])>0 :
                                            for elementadd in range(len(cellidx[lx][trenchNum] )):
                                                list1 = ygg_cells[lx][trenchNum]
                                                list1.insert(beforeElement[lx][trenchNum][elementadd], ygg_cells[lx][trenchNum][extraa[lx][trenchNum][elementadd]]) 
                                                ygg_cells[lx][trenchNum] = list1   
                                                
                                                list1 = cellbirth[lx][trenchNum]
                                                list1.insert(beforeElement[lx][trenchNum][elementadd], cellbirth[lx][trenchNum][extraa[lx][trenchNum][elementadd]]) 
                                                cellbirth[lx][trenchNum] = list1  
                                                
          
        
                         t[trenchNum]+=[t1[trenchNum]]
                         cellnumberarray [l][trenchNum]  = [ i for i in range (cellnn[trenchNum][0])]
                         if totlengths[trenchNum] >2.:
                                
                                lx, ll[trenchNum], sum__= [l , [], 0]
                                for mm in range(int(len(t[trenchNum][lx] )/2)):
                                    sum__ += t[trenchNum][lx] [mm*2+1] 
                                    if sum__ <limm:
                                        ll[trenchNum] += [t[trenchNum][lx] [mm*2+1]]
                           
                  alllengths[l],alllgrowths[l], h2o2_[l]  ,katgall[l]  =[ totlengths2,totgrowths2, h2o2_2, katgall2 ] 
                 
              return  cellnn ,alllengths,cellnumberarray    , alllgrowths, h2o2_, katgall 
    
     
    limm= limm*10.#
    gI = 0.012*1.*1.; numberoftrench =len(trenches_)

    conx = conx*um
     
    trenches=[]
    for tr in trenches_:
        tr_ = [ i*10.  for i in tr]
        trenches+=[tr_]

    kc =      kcats  
    conc1 = np.array([bH for i in range(treated)]); conc2 = np.array([conx for i in range(duration-treated)]) ;concallt = np.concatenate((conc1,conc2))
    nfr=2; timee = np.linspace(0,1,nfr) 
    y0x=[bH, bKatG ] 
    yggput = odeint(abstract_response_model.comb_model3_growth,y0x,timee,args=([i for i in range(0,nfr)],concallt[0:nfr] ,gI, kgi, cI,kc,Kenzyme))  
    gett_, alllengths,cellnumberarray ,alllgrowths, h2o2_, katgall = tps(trenches , [i for i in range(duration)], [i * 150 for i in range(0,numberoftrench) ],yggput,concallt,kc,Kenzyme )
    gett_2 = modelfunc.tps2([trenches ] , duration, [i * 150 for i in range(0,numberoftrench) ] , gett_,alllengths,cellnumberarray,alllgrowths ,h2o2_,katgall,numberoftrench)

    arrays=[[]]
    for time in range((duration)-2):
        for eachtrenchg in range(len(alllengths[time])):
            for j in range(0,  len( alllengths[time][eachtrenchg]) ):
                arrays+=[[time, eachtrenchg, j, alllgrowths[time][eachtrenchg][ j]   ,   katgall  [time][eachtrenchg][ j]       ,h2o2_ [time][eachtrenchg][ j]    ,  len( katgall[time][eachtrenchg]) ,alllengths[time][eachtrenchg][ j] , conx]]
    df = pd.DataFrame(arrays, columns =['Time','TrenchNumber','PositionClosedEnd', 'Elongation', 'Enzyme', 'Toxin','max_cells','length','concentration'])
    df['Frame'] = df['Time']
    
    plotting.plotenzyme(df,numberoftrench,treated,duration); 
    if duration > treated+600:
        plotting.plotpoincare(df,numberoftrench,treated,duration); 
    
         
  
