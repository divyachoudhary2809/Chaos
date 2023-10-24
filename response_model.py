# -*- coding: utf-8 -*-
"""
@author: some4879
"""
import numpy as np
#from numba import jit
# units in uM / minute   
 
bH =0.02
 
OxyRtot = 1.
bGrxA = 0.00;
bKatG= 0.00;
bAhpC = 0.01*1.;

kcat =  660.*60. 
ah2= 1.2  
kh2= 5900. 
kcat2 =490000.*60. 

Kox =10.
Kred = 8.*60. 
koxx=2583.

 
HAhpC=0.1
KAhpC=2.*0.1
HKatG =0.6 *3.
KKatG =5. *3.*0.1
HGrxA = 0.1
KGrxA =1.*0.1

 
#@jit(nopython=True)
def comb_model3_growth(x,t,u,KHa,gI, kgi, cI):    

    KH = KHa[0]
 
    dHdt = KH + bH- kcat*x[4]*x[0]*(1./(x[0] +ah2))  - kcat2  *x[0] *x[3]*(1./( x[0]+kh2))   
    g = gI-(  gI / (1 + 10**(-kgi*(x[0]-cI)))   ) 
    gamma =  g 
    dOxyReddt =  -   Kox* (x[1])     *x[0]  + Kred *x[2]* (OxyRtot- x[1]) *(1./((OxyRtot- x[1]) +koxx))    
    dGrxAdt = bGrxA + (KGrxA * ((OxyRtot - x[1]) /((OxyRtot- x[1] ) + HGrxA )) )- gamma * x[2]  
    dKatGdt = bKatG + (KKatG * ((OxyRtot - x[1]) /((OxyRtot- x[1] ) + HKatG )) )- gamma * x[3]  
    dAhpcCdt =  bAhpC + ((KAhpC) * ((OxyRtot - x[1]) /((OxyRtot- x[1] ) + HAhpC )) ) - gamma * x[4] 
 
    dzdt = [dHdt,dOxyReddt,dGrxAdt,dKatGdt, dAhpcCdt] 
    return dzdt  
