
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 27 19:59:46 2023

@author: divyachoudhary
"""
from numba import jit
import numpy as np
import matplotlib.pyplot as plt


bH =0.02 ; OxyRtot = 1.000; bGrxA = 0.00; bKatG= 0.00;
ah2= 1.2  ; kh2= 5900. ; koxx=2583.
HAhpC=0.1 ; HKatG =0.6 *3. ;  HGrxA = 0.1

bAhpC = 0.01*1.;
kcat =  660.*60.    ; kcat2 =490000.*60.  ;
Kox =10.*1 ; Kred = 8.*60.*1
KAhpC=0.2 *1; KKatG =1.5*1 ; KGrxA =0.1*1
 
 
gI = 0.042/3.5 ; um =   (10**-6);
kgi= .02/um ; cI= 100.*um

dt =0.000005  # Time step.
T =  1.# Total time.
n = int(T / dt)  # Number of time steps.
t = np.linspace(0., T, n)  # Vector of times.

 

@jit(nopython=True)
def runmain(lx,c,gI , h,o,gr,k,a, sigmaa):

    H = np.zeros(n)  ; OxyR= np.zeros(n) ; GrxA= np.zeros(n) ; KatG= np.zeros(n);AhpC= np.zeros(n) ;gammas= np.zeros(n)  
    
    H[0] =h; OxyR[0] =o ;GrxA[0] = gr; KatG[0] = k; AhpC[0] = a; gammas[0] =gI
    sigma = [ 0.00000  , sigmaa /100.,sigmaa , sigmaa, sigmaa]## no noise for h2o2
    
    for conx in [c]:
        treated = int(0/dt);
        conc1 = np.array([bH for i in range(treated)]); conc2 = np.array([conx for i in range(n-treated)]) ;concallt = np.concatenate((conc1,conc2))
        KH = concallt[0:n]
        for i2 in range(1, n ):
     
            i=i2-1
            H[i+1] = H[i]  + (dt) * (  ( KH [i] + bH  ) - (kcat*AhpC[i]*H[i]*(1./(H[i] +ah2)) ) - (kcat2  *H[i] *KatG[i]*(1./( H[i]+kh2)) )  )   
            
            gamma =  gI    -(  gI / (1 + 10**(-kgi*(H[i+1]-cI)))   )
            tau = 1/gI
            dt_nd = dt/tau
            dW = np.random.normal(0, np.sqrt(dt_nd))  
            gammas[i+1] =gamma
            OxyR[i+1] = OxyR[i] + dt_nd*(1/gI)* ( -   Kox* (OxyR[i])     *H[i]  + Kred *GrxA[i]* (OxyRtot- OxyR[i]) *(1./((OxyRtot- OxyR[i]) +koxx))   )    + sigma[1]   *dW 
            GrxA[i+1] = GrxA[i]  + dt_nd*(1/gI)*( bGrxA + (KGrxA * ((OxyRtot - OxyR[i]) /((OxyRtot- OxyR[i] ) + HGrxA )) - gamma*GrxA[i]))    + sigma[2]   *dW 
            KatG[i+1] = KatG[i]  +  dt_nd*(1/gI)*(bKatG + (KKatG * (OxyRtot - OxyR[i]) /((OxyRtot- OxyR[i] ) + HKatG )) - gamma * KatG[i]  )     + sigma[3]   *dW 
            AhpC[i+1] = AhpC[i]   + dt_nd*(1/gI)*( bAhpC + ((KAhpC) * ((OxyRtot -OxyR[i]) /((OxyRtot- OxyR[i] ) + HAhpC )) ) - gamma * AhpC[i])    + sigma[4]   *dW  

            
    return    H,   OxyR, GrxA, KatG, AhpC
 
