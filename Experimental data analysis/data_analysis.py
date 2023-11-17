import numpy as np
import pandas as pd  
from matplotlib import pyplot as plt
from sklearn.preprocessing import normalize
from statsmodels.tsa import stattools
from scipy.signal import find_peaks 


files = [  ['File1', [Switchtime1]]  , ['File2', [Switchtime2]]  ] #### Folder names of BACMMAN output files, followed by the frame of switch to treatment
path_= "Path to Directory" ### Path to the BACMMAN folder

for file_i in  range (len(files)): 
    
    fileID= files[file_i] ;  tswitch =fileID[1]

    #####Loading the data sets and computing different features of a cell
    
    fn ='/%s/%s_1.csv'%(fileID[0], fileID[0])
    data_allxx2cfp= pd.read_csv( path_ + fn[:-5] +'2.csv', sep=';')  ### CFP fluorescence data 
    data_allxx2cfp['channel'] =[-1]*len(data_allxx2cfp)
    for dd in range(len(data_allxx2cfp)):
        data_allxx2cfp.Position.values[dd] = int(data_allxx2cfp.Position.values[dd].split('_xy')[1]) 
        data_allxx2cfp.channel.values[dd] = int(data_allxx2cfp.Indices.values[dd].split('-')[1])
     
    
    data_allxx2rfp= pd.read_csv( path_ + fn[:-5] +'1.csv', sep=';')  ### mKate2 segmentation and fluorescence data 
    data_allxx2rfp['channel'] =[-1]*len(data_allxx2rfp)
    for dd in range(len(data_allxx2rfp)):
        data_allxx2rfp.Position.values[dd] = int(data_allxx2rfp.Position.values[dd].split('_xy')[1])
        data_allxx2rfp.channel.values[dd] = int(data_allxx2rfp.Indices.values[dd].split('-')[1])
    

    poss_x2 =[] ; frame_2= [] ; channel_2= [] ; max_cells2 = []
    for xf in np.unique(data_allxx2rfp.Position.values):
        data_before = data_allxx2rfp[data_allxx2rfp.Position == xf]   
        data_before = data_before.reset_index(drop=True)
        for fr in np.unique(data_before.Frame):
            data_beforex= data_before[data_before.Frame.values==fr] 
            for ch in np.unique(data_before.channel):
                data_beforey= data_beforex[data_beforex.channel.values==ch] 
                if len(data_beforey)>0:
                    poss_x2 +=[xf] ; frame_2 +=[fr] ; channel_2 +=[ch] ;  max_cells2 +=[max(data_beforey.Idx)]
    max_cells_2= pd.DataFrame(np.array([poss_x2,frame_2,channel_2,max_cells2]).T)
    max_cells_2.columns= ['Position','Frame','channel','max_cells2']
    max_cells_2['combined'] = max_cells_2['Position'].map(str) + '-' + max_cells_2['Frame'].map(str) + '-' + max_cells_2['channel'].map(str)
    data_allxx2rfp['ind_new'] = ['']*len(data_allxx2rfp)
    for dd in range(len(data_allxx2rfp)):
        data_allxx2rfp['ind_new'].values[dd] = (data_allxx2rfp.Indices.values[dd].split('-')[0]) +'-' +(data_allxx2rfp.Indices.values[dd].split('-')[1])
    data_allxx2rfp['combined'] = data_allxx2rfp['Position'].map(str) + '-' + data_allxx2rfp['ind_new'].map(str)  
    data_allxx2rfp.sort_index(inplace=True)
    dictt2 = pd.Series(max_cells_2.max_cells2.values,index=max_cells_2.combined).to_dict()
    data_allxx2rfp['max_cells'] = data_allxx2rfp['combined'].map(dictt2)
    data_allxx2rfp['idx2'] = data_allxx2rfp.max_cells - data_allxx2rfp.Idx  
 


    ### 
    data_all_feat= pd.DataFrame() ## This dataframe contains all the features of the cells (morphology and fluorescence intensity) combined 
    for p in np.unique(data_allxx2rfp.Position.values):
                data_allx = data_allxx2rfp[data_allxx2rfp.Position.values ==p]
                data_all2x = data_allxx2cfp[data_allxx2cfp.Position.values ==p]
                dictt12 = pd.Series(data_all2x.CFPMeanIntensity.values,index=data_all2x.BacteriaIndices.values).to_dict()
                data_allx['MeanIntensity_CFP'] =[-1] * len(data_allx)
                data_allx['MeanIntensity_CFP'] = data_allx['Indices'].map(dictt12)    
                data_all_feat = pd.concat([data_all_feat, data_allx])
 
  
    
  
    





    ####  ACF of DPgrxA for single mother cell traces overlaid with mean ACF trace
    
    plt.figure(figsize=(2,2), dpi=300)
    tim= max(data_all_feat.Frame.values) - tswitch[0] - 20 ; taft_tr=20 
    df = data_all_feat 
    autos= [] ;   n=0
    for p in np.unique(data_all_feat.Position.values)    :
        dfp = data_all_feat[data_all_feat.Position.values ==p]
        for c in np.unique(dfp.channel.values)   :

                dfpx = dfp[dfp.channel.values ==c] ;  dfpm = dfpx[dfpx.Idx.values==0]
                dfpm['MeanIntensity_CFP_relativetobasal'] =dfpm['MeanIntensity_CFP'] - np.mean(dfpm [dfpm.Frame.values< tswitch[0]].MeanIntensity_CFP.values)
                dfpm['MeanIntensity_CFP_relativetobasal'] =dfpm['MeanIntensity_CFP_relativetobasal'].rolling(3).mean()
                dfpm['MeanIntensity_CFP_relativetobasal'] =dfpm['MeanIntensity_CFP_relativetobasal'].rolling(3).mean()
                dfpm['MeanIntensity_CFPdiff'] =dfpm['MeanIntensity_CFP_relativetobasal'] - dfpm['MeanIntensity_CFP_relativetobasal'].shift(1)
                dfpm = dfpm.dropna(subset= ['MeanIntensity_CFPdiff']) ; grrates = np.sort( dfpm.GrowthRateLength.values ) 
                if len(dfpm) > 0 and len(dfpm.GrowthRateLength.values[:-10]) >0 and np.mean(grrates [0: int(len(grrates)/4)])>0.00012*0.5 : 
 
                
                    dfpmtt= dfpm[dfpm.Frame.values > tswitch[0]+taft_tr] 
                    dfpmtt= dfpmtt[dfpmtt.Frame.values < tswitch[0]+taft_tr+tim] 
                    x = dfpmtt.MeanIntensity_CFPdiff.values     
                 
                    if len(x) >0:
                        autocorr = stattools.acf( x, nlags=tim-25)
                        if len(autocorr) >=tim-26:
                            n+=1
                            if n<100:
                                plt.plot([i*3 for i in range(len(autocorr[0:tim-26]))], autocorr[0:tim-26], color='k', lw=0.05)
                                autos+=[autocorr[0:tim-26] ]  
    autoall =   np.mean(autos, axis=0) 
    plt.plot([i*3 for i in range(len(autoall))],autoall,  color='darkcyan', label='Mean ACF')
    plt.ylabel('PgrxA ACF [a.u.]', fontsize=15)
    plt.xlabel('Lags [min]', fontsize=15)
    plt.xlim(0,200)
    plt.legend(frameon=False)
 
    
    
    
    
    
    
    
    
    
    
    ##### Mean Elongation rate for cells at different positions in the growth channel

   
    cx= [ 'k','darkgreen','green','springgreen','gold']
    data_all_feat2 = data_all_feat[data_all_feat.Frame.values < tswitch[0]+140]
    image,ax =plt.subplots(figsize=(2,2), dpi=300)
    dab=pd.DataFrame()
    for lpos in range(4):
        tt= []
        elongs_= []
        for ppp in np.unique(data_all_feat2['Position'].values)   :
            data_allxx2a= data_all_feat2[data_all_feat2['Position'].values == ppp]
            for thi in  np.unique(data_allxx2a.TrackHeadIndices.values): 
                     data_allxx2ab = data_allxx2a[data_allxx2a['TrackHeadIndices'].values == thi]
                     data_allxx2ab['elonrate_'] = [-1.00000000]* len(data_allxx2ab )
                     for d in range(1, len(data_allxx2ab)):
                         if (data_allxx2ab.Frame.values[d] - data_allxx2ab.Frame.values[d-1]) !=0 :
                            data_allxx2ab['elonrate_'].values[d] =  (  np.log(data_allxx2ab.Length.values[d]) -  np.log(data_allxx2ab.Length.values[d-1]))/(data_allxx2ab.Frame.values[d] - data_allxx2ab.Frame.values[d-1])
                     data_allxx2ab= data_allxx2ab[data_allxx2ab.elonrate_.values!=-1] ; data_allxx2abx = data_allxx2ab [data_allxx2ab.idx2.values >lpos-1] ; data_allxx2abx = data_allxx2abx [data_allxx2abx.idx2.values <lpos+1]
                     if len(data_allxx2abx) >2: 
                            tt+=[data_allxx2abx.Frame.values[1:]- tswitch[0]]
                            elongs_ += [ data_allxx2abx.elonrate_.values[1:] ]
        time_=[] 
        for to in tt:
            for tx in to:
                time_+= [tx]
        
        elongs=[] 
        for e in elongs_:
            for ee in e:
                elongs+= [ee]
        
        df_elong= pd.DataFrame()
        df_elong['time'] = time_; df_elong['Elongation_rate_%s'%lpos] = elongs
        df_elong_grouped=df_elong .groupby('time').mean()  ; df_elong_grouped.index = (df_elong_grouped.index  )*3
 
        plt.plot(df_elong_grouped.index,df_elong_grouped['Elongation_rate_%s'%lpos], color = cx[lpos], lw=1, label =lpos)#, s=1)"
        plt.ylabel('Elongation \nrate [min$^{-1}$]', fontsize=15 )
        plt.xlabel('Time [min]', fontsize=15)
        plt.xlim(-50,300)
        plt.ylim( -0.003,0.085 )
        
    plt.legend(title='Number of\nbarrier cells', frameon=False , bbox_to_anchor=[1.3, 0.5], loc='center')
             
     
    
    
    
    
    
    
    ##### Mean PgrxA expression values for cells at different positions in the growth channel
    
 
    cx= ['gold','springgreen','green','darkgreen','k']
    image,ax =plt.subplots(figsize=(2,2), dpi=300)
    for colname in  ['MeanIntensity_CFP'] : 
         
                    data_allxx2abx = data_all_feat 
                    dtut = data_allxx2abx[data_allxx2abx.Frame.values<tswitch[0]]
                    dtut = dtut.dropna(subset=['MeanIntensity_CFP'])
                    dtutcfp = np.mean(dtut.MeanIntensity_CFP.values)
                    for tpos in range(0,4):
                        index_col = 'idx2'  
                        data_allxx2abx2 = data_allxx2abx[data_allxx2abx['%s'%index_col] == tpos]
                        data_allxx2abx_plot = data_allxx2abx2.groupby('Frame').mean()
                        plt.plot((data_allxx2abx_plot.index- tswitch[0])*3, data_allxx2abx_plot['%s'%colname]-dtutcfp, color=cx[tpos], linewidth=1, label=tpos)
 
    plt.xlim(-50,450)  
    plt.ylabel('PgrxA [a.u.]', fontsize=15 )
    plt.xlabel('Time [min]', fontsize=15)
    plt.legend(title='Number of\nbarrier cells', frameon=False , bbox_to_anchor=[1.3, 0.5], loc='center')
    
    
    
    
    
    #### Plot of individual GrxA expression traces of mother cells before and after treatment
    
    colors =['darkblue', 'dodgerblue', 'teal', 'k', 'cyan']
    n=0 ; fig, ax = plt.subplots(figsize=(5,2), dpi = 300)
    for p in np.unique(data_all_feat.Position.values)  :
        
        
        for ch in range(26):
            df_ = data_all_feat[data_all_feat.Position.values==p]
            df_2= df_[df_.channel.values==ch]
            df_2 = df_2[df_2.Frame.values< tswitch[0] +200]
            df_2m = df_2[df_2.Idx.values==0]

            grrates = np.sort( df_2m.GrowthRateLength.values )   
            if len(df_2m)>200 and np.mean(grrates [0: int(len(grrates)/4)])>0.00012*0.5  and len(np.unique(df_2m.ParentTrackHeadIndices.values)) <4  : 
                n+=1
                if n<5:
                    ax.plot((df_2m.Frame.values - tswitch[0])/20 , df_2m.MeanIntensity_CFP-1410, lw=1, color= colors[n-1])#, color='k')
                    ax.set_xlabel('Time [hours]', fontsize=15)
                    ax.set_ylabel('PgrxA [a.u.]', fontsize=15)
                    
 
    
 
    
    
    
    ######Plotting PgrxA traces for alive and dead cells

    dfchaosa=pd.DataFrame() ; dfchaosa_poincare=pd.DataFrame() 
    dfchaosd=pd.DataFrame() ; dfchaosd_poincare=pd.DataFrame() 
    
    tim=210
    taft_tr=50 
    
    alive=[]
    dead=[]
    autos= []
    n=0

    for p in np.unique(data_all_feat.Position.values)    :
        for c in np.unique(data_all_feat.channel.values)   :
            for i in range(1):
                dfp = data_all_feat[data_all_feat.Position.values ==p]
                dfpx = dfp[dfp.channel.values ==c]
                dfpm = dfpx[dfpx.Idx.values==i]
                dfpm= dfpm[dfpm.Frame.values < tswitch[0]+taft_tr+tim] 
                dfpmut = dfpm[dfpm.Frame<tswitch[0]]
                if len(dfpm) >tswitch[0]+tim    and len(np.unique(dfpm.ParentTrackHeadIndices.values)) <2  and max(dfpm.Length.values) <15. :
                    n+=1
                    dfpmtt= dfpm[dfpm.Frame.values > tswitch[0]+taft_tr] 
           
                    x = dfpmtt.MeanIntensity_CFP .rolling(3).mean()    [2:].values
                    t = dfpmtt.Frame.rolling(3).mean()    [2:].values

                    x2 = dfpm.MeanIntensity_CFP.values  
                    t2 = dfpm.Frame.values  
                    dfn =pd.DataFrame() 

                    grrates = np.sort( dfpmtt.GrowthRateLength.values )
                    if len(dfpmtt.GrowthRateLength.values[tswitch[0] +taft_tr:tswitch[0]+taft_tr+tim])>0 and np.mean(grrates [0: int(len(grrates)/4)])>0.00012*0.5:#min(dfpmtt.GrowthRateLength.values[tswitch[0] +taft_tr:tswitch[0]+taft_tr+tim] )>0.00012*0.5 :
                            
                            if len(x2)>0:
                                alive+=[n-1]
                                dfn =pd.DataFrame()  ; dfn_poincare =pd.DataFrame()  
                                dfn['cfp-%s-%s'%(p,c)] =x2    ; dfn_poincare['cfp-%s-%s'%(p,c)] =x    [~np.isnan(x)]                           
                                dfchaosa=pd.concat([dfchaosa,dfn] , axis=1); dfchaosa_poincare=pd.concat([dfchaosa_poincare,dfn_poincare] , axis=1)
                                
                                dfn =pd.DataFrame()  ; dfn_poincare =pd.DataFrame()  
                                dfn['cfpt-%s-%s'%(p,c)] =t2     ; dfn_poincare['cfpt-%s-%s'%(p,c)] =t     [~np.isnan(t)]                       
                                dfchaosa=pd.concat([dfchaosa,dfn] , axis=1); dfchaosa_poincare=pd.concat([dfchaosa_poincare,dfn_poincare] , axis=1)
                            
                    else:
                            
                            if len(x2)>0:
                                dead+=[n-1]
                                dfn =pd.DataFrame()  ; dfn_poincare =pd.DataFrame()  
                                dfn['cfp-%s-%s'%(p,c)] =x2    ; dfn_poincare['cfp-%s-%s'%(p,c)] =x [~np.isnan(x)]                 
                                dfchaosd=pd.concat([dfchaosd,dfn] , axis=1)  ; dfchaosd_poincare=pd.concat([dfchaosd_poincare,dfn_poincare] , axis=1)
                                
                                dfn =pd.DataFrame()  ; dfn_poincare =pd.DataFrame()  
                                dfn['cfpt-%s-%s'%(p,c)] =t2     ; dfn_poincare['cfpt-%s-%s'%(p,c)] =t     [~np.isnan(t)]                      
                                dfchaosd=pd.concat([dfchaosd,dfn] , axis=1); dfchaosd_poincare=pd.concat([dfchaosd_poincare,dfn_poincare] , axis=1)
                                
                            
                
     
    plt.figure(figsize=(4,2), dpi = 300)
    totalplot=20
    for n_cells in range(totalplot ) :
        plt.plot([ i*3- tswitch[0]*3 for i in range (len(dfchaosa.iloc[:,n_cells*2+1]))], dfchaosa.iloc[:,n_cells*2]-1410, lw=1)
        
    plt.xlabel(' Time [min]', fontsize=15)
    plt.ylabel(' PgrxA [a.u.]', fontsize=15)
    plt.title('Alive cells')
         
    plt.figure(figsize=(4,2), dpi = 300)
 
    for n_cells in range( int ( len(dfchaosd.columns) *0.5) -1):
        plt.plot([ i*3- tswitch[0]*3 for i in range (len(dfchaosd.iloc[:,n_cells*2+1]))],  dfchaosd.iloc[:,n_cells*2]-1410, lw=1)   
    
    plt.xlabel(' Time [min]', fontsize=15)
    plt.ylabel(' PgrxA [a.u.]', fontsize=15)
    plt.title('Dead cells')
    
    
    
    ###Poincare plots of alive and dead cells for nth cell to n+10th cell
    
    n=2
    fig, axs = plt.subplots(2,5, figsize=(15, 7) , dpi = 300)
    fig.subplots_adjust(hspace = .5, wspace=.5)
    
    axs = axs.ravel()
    
    for i in range(n, n+10):
        if len(dfchaosa_poincare.iloc[10:,i*2]) >0:
            axs[i-n].plot( dfchaosa_poincare.iloc[:-10,i*2]-1410 ,  dfchaosa_poincare.iloc[10:,i*2]-1410, color='teal')
            axs[i-n].set_xlim(min (dfchaosa_poincare.iloc[10:,i*2]-1410) -100,800)
            axs[i-n].set_ylim(min (dfchaosa_poincare.iloc[10:,i*2]-1410) -100,800)
            axs[i-n].set_xlabel(' PgrxA at t [a.u.]', fontsize=15)
            axs[i-n].set_ylabel(' PgrxA at t+30 min [a.u.]', fontsize=15)
            if i ==2:
                axs[i].set_title('Alive cells', fontsize=20)
        
        
        
    fig, axs = plt.subplots(2,5, figsize=(15, 7) )
    fig.subplots_adjust(hspace = .5, wspace=.5)
    
    axs = axs.ravel()
    
    for i in range(n, n+10):
        if len(dfchaosd_poincare.iloc[10:,i*2]) >0:
            axs[i-n].plot( dfchaosd_poincare.iloc[:-10,i*2]-1410 ,  dfchaosd_poincare.iloc[10:,i*2]-1410, color='crimson')
            axs[i-n].set_xlim(0,800)
            axs[i-n].set_ylim(0,800)
            axs[i-n].set_xlabel(' PgrxA at t [min]', fontsize=15)
            axs[i-n].set_ylabel(' PgrxA at t+30 min [a.u.]', fontsize=15)
            if i ==2:
                axs[i].set_title('Alive cells', fontsize=20) 
        
    
                
    ### Interdivision time versus GrxA ACF peak time for individual cells at steady state
    
    allpeaks= []
    alldivs = []
 
    df= data_all_feat 
    df['divT'] = df['NextDivisionFrame'] - df['PreviousDivisionFrame']

    tim=max(df.Frame.values)-tswitch[0]
    taft_tr=20

    autos= []; n=0;m=0
    print(max(df.Frame.values)-tswitch[0])
    
    for p in np.unique(df.Position.values)    :
        dfp = df[df.Position.values ==p]
        for c in np.unique(df.channel.values)   :
            
                dfp = df[df.Position.values ==p]
                dfpx = dfp[dfp.channel.values ==c]
                dfpm = dfpx[dfpx.Idx.values==0]
                dfpm['MeanIntensity_CFP'] =dfpm['MeanIntensity_CFP'] - np.mean(dfpm [dfpm.Frame.values< tswitch[0]].MeanIntensity_CFP.values)
                dfpm['MeanIntensity_CFP'] =dfpm['MeanIntensity_CFP'].rolling(3).mean()
                dfpm['MeanIntensity_CFP'] =dfpm['MeanIntensity_CFP'].rolling(3).mean()
                dfpm['MeanIntensity_CFPdiff'] =dfpm['MeanIntensity_CFP'] - dfpm['MeanIntensity_CFP'].shift(1)
                dfpm = dfpm.dropna(subset= ['MeanIntensity_CFPdiff'])
                dfpmtt= dfpm[dfpm.Frame.values > tswitch[0]+taft_tr]  ; dfpmtt= dfpmtt[dfpmtt.Frame.values < tswitch[0]+taft_tr+tim] 
                grrates = np.sort( dfpmtt.GrowthRateLength.values )
               
                if len(dfpmtt) > 20  and np.mean(grrates [0:int(len(grrates)/4)])>0.00012*0.5  :
                    x = dfpmtt.MeanIntensity_CFPdiff.values     
                    x = x[~np.isnan(x)]
                    if len(x) >0:
                        autocorr = stattools.acf( x, nlags=tim-25)
                        if len(autocorr) >= 0:
                            autos+=[autocorr]
                            peaks, _ = find_peaks(autocorr, prominence=0.2, width =1)
                            dfpmtt = dfpmtt.dropna(subset = ['divT'])
                            if len(peaks) >0 and [ i *3. for i in range( len(autocorr))][peaks[0]] <100:
                                if str(np.median(dfpmtt.divT.values)) !='nan':
                                    allpeaks+=[[ i *3. for i in range( len(autocorr))][peaks[0]]]
                                    alldivs+=[np.mean(dfpmtt.divT.values)*3.]
                                    
    allpeaks = np.array(allpeaks) ; alldivs = np.array(alldivs)
    plt.figure(figsize=(2,2), dpi = 300)
    plt.errorbar(np.mean(allpeaks), np.mean(alldivs), xerr=np.std(allpeaks), yerr=np.std(alldivs), fmt="o", color="darkorange")
    plt.xlim(-10,160)
    plt.ylim(-10,160)
    plt.xlabel(' ACF peak time [min]', fontsize=15)
    plt.ylabel(' Interdivision \n time [min]', fontsize=15)
    
     
 
    
 

 
