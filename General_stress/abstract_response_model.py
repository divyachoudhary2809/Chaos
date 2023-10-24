import numpy as np

bH = 0.00
benzyme= 0.00;
Henzyme =0.1
ah2= 1.  
 
def comb_model3_growth(x,t,u,KHa,gI, kgi, cI,kcat,Kenzyme):
 
    g = gI-(  gI / (1 + 10**(-kgi*(x[0]-cI)))   ) 
    gamma =  g  
    KH=bH
 
    dHdt = KH + bH   - kcat  *x[0] *x[1]*(1./(x[0] +ah2))    
    denzymedt =  benzyme + (Kenzyme * ((x[0])**1/((x[0])**1+ Henzyme**1)) )- gamma * x[1] 
    dzdt = [dHdt, denzymedt ]
    return dzdt  


