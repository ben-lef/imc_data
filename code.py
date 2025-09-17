# -*- coding: utf-8 -*-

"""
This code contains functions for thorough analysis of the intermetallic compound
layer. From the corresponding database, it computes micro, macro and discontinuity
descriptors for the desired sample.
"""

# Benjamin Leflon (INSA Lyon, Tohoku University, ELyTMaX)
# Contact: benjamin.leflon96@gmail.com
# Last update: 17/09/2025


#%% Modules

from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.interpolate import make_smoothing_spline
from matplotlib.ticker import FixedLocator, FixedFormatter, FormatStrFormatter

#%% Reading database

code_dir = Path(__file__).parent
rel_path = 'data.csv'

data_path = (code_dir / rel_path).resolve()

db = pd.read_csv(data_path)


#%% Function

def full_analysis(db, conds, exp_label, section, plot=True):
    
    """
    Inputs:
        db: DataFrame containing the data
        conds: label of the welding conditions (i.e. 'i140v170')
        exp_label: label of the experiment (1, 2, 3)
        section: name of the cross section to analyze ('T30', 'L1' or 'L2')
        plot: if True, plots the evolution of several descriptors
        
    Output:
        a dict containing every computed local, global, micro and macro descriptors of
        the analyzed layer, alongside the list of corresponding positions
        
        
    Use case:
        to analyze the layer obtained for the 1st realization of the i140v140
        welding condition, for the transverse cross section at z = 30 mm, the
        command is:  
            full_analysis(db, 'i140v140', 1, 'T30')
    """
    
    dat = db.loc[(db.loc[:,'label']==conds) & (db.loc[:,'exp_label']==exp_label) & (db.loc[:,'section']==section),:]
    
    if dat.empty:
        raise Exception('There is no data corresponding to these conditions')
    
    echelle = dat.loc[(dat.loc[:,'image_label'] == 1), 'scale'].max()
    
    mean_positions = dat.loc[:,['image_label','global_position']].groupby(by='image_label').agg('mean')
    mean_positions['global_position'] = mean_positions['global_position']*echelle/1E3
    
    pos = mean_positions['global_position'].to_numpy()
    

    
    
    
    non_nan = dat.loc[dat.loc[:,'dense_thickness'].isna() == False,:]    
    
    overlap = 0.2
    len_pics = dat['local_position'].max()
    
    non_overlap = non_nan.loc[non_nan.loc[:,'local_position']<=int((1-overlap)*len_pics),'dense_thickness'].to_numpy()
    
    if non_overlap[0]==0:
        disc = 1
    else:
        disc = 0
        
    for i in range(1,len(non_overlap)):
        if non_overlap[i-1] !=0:
            if non_overlap[i] ==0:
                disc = disc+1
        
    
    is_imc = non_nan.loc[non_nan.loc[:,'dense_thickness']!=0,:]
                           
    is_imc.loc[:,'dense_thickness'] = echelle*is_imc.loc[:,'dense_thickness']
    is_imc.loc[:,'envelope_thickness'] = echelle*is_imc.loc[:,'envelope_thickness']
    is_imc.loc[:,'global_position'] = echelle*is_imc.loc[:,'global_position']/1E3
    
    
    is_IMC_pos = is_imc['global_position'].to_numpy()
    
    indicators = {}
    
    for i in set(dat.loc[:,'image_label']):
        indicators[i] = {}
        
            
        if len(is_imc.loc[is_imc.loc[:,'image_label']==i])==0:
            indicators[i]['MAD'] = np.nan
            indicators[i]['C'] = np.nan
            indicators[i]['RaFe'] = np.nan
            indicators[i]['RaAl'] = np.nan
            indicators[i]['smFe'] = np.nan
            indicators[i]['smAl'] = np.nan
            indicators[i]['SfFe'] = np.nan
            indicators[i]['SfAl'] = np.nan
            
            
        else:

                
            dense = is_imc.loc[is_imc.loc[:,'image_label']==i,'dense_thickness']
            envelope = is_imc.loc[is_imc.loc[:,'image_label']==i,'envelope_thickness']
            
            MAD = (dense-dense.mean()).abs().mean()
            indicators[i]['MAD'] = MAD
    
            C = 1 - ((envelope-dense)/envelope).mean()
            indicators[i]['C'] = C
        
            #local_positions = is_imc.loc[is_imc.loc[:,'image_label']==i,'local_position']
            
            Al = is_imc.loc[is_imc.loc[:,'image_label']==i,'al_interface']
            Fe = is_imc.loc[is_imc.loc[:,'image_label']==i,'fe_interface']
        
            # regAl = stats.linregress(local_positions, Al)
            # meanlineAl = regAl.intercept + regAl.slope*local_positions
            meanlineAl = Al.mean()
            RaAl = (Al - meanlineAl).abs().mean()
            
            # regFe = stats.linregress(local_positions, Fe)
            # meanlineFe = regFe.intercept + regFe.slope*local_positions
            meanlineFe = Fe.mean()
            RaFe = (Fe - meanlineFe).abs().mean()
            
            indicators[i]['RaAl'] = RaAl
            indicators[i]['RaFe'] = RaFe
            
            crossings = []
            for j in range(len(Al)-1):
                if (Al.iloc[j+1] >=meanlineAl) & (Al.iloc[j] <meanlineAl):
                    crossings.append(j)
                    
            crossings = np.array(crossings)
            si = np.diff(crossings)
            smAl = np.mean(si)
            
            indicators[i]['smAl'] = smAl
            indicators[i]['SfAl'] = RaAl/smAl
            
            
            crossings = []
            for j in range(len(Al)-1):
                if (Fe.iloc[j+1] >=meanlineFe) & (Fe.iloc[j] <meanlineFe):
                    crossings.append(j)
                    
            crossings = np.array(crossings)
            si = np.diff(crossings)
            smFe = np.mean(si) 
            
            indicators[i]['smFe'] = smFe
            indicators[i]['SfFe'] = RaFe/smFe
      
    RaAl = []
    RaFe = []
    smAl = []
    smFe = []
    SfAl = []
    SfFe = []
    C = []
    MAD=[]
        
    for i in indicators.keys():
        RaAl.append(echelle*indicators[i]['RaAl'])
        RaFe.append(echelle*indicators[i]['RaFe'])
        smAl.append(echelle*indicators[i]['smAl'])
        smFe.append(echelle*indicators[i]['smFe'])
        SfAl.append(indicators[i]['SfAl'])
        SfFe.append(indicators[i]['SfFe'])
        C.append(indicators[i]['C'])
        MAD.append(indicators[i]['MAD'])
    
    if plot:
        fig, ax = plt.subplots()
        ax.plot(pos,RaAl,label='$R_{a}^{Al}$')
        ax.plot(pos,RaFe,label='$R_{a}^{Fe}$')
        ax.plot(pos,MAD,label='MAD')
        ax.legend() 
        ax.set_xlabel('Transverse position (mm)', fontdict={'size':14})
        ax.set_ylabel('Interface $R_{a}$, thickness MAD (µm)', fontdict={'size':14})
        
        fig, ax = plt.subplots()
        ax.plot(pos,smAl,label='$s_{m}^{Al}$')
        ax.plot(pos,smFe,label='$s_{m}^{Fe}$')
        ax.legend() 
        ax.set_xlabel('Transverse position (mm)', fontdict={'size':14})
        ax.set_ylabel('Interface $s_{m}$ (µm)', fontdict={'size':14})
        
        fig, ax = plt.subplots()
        ax.plot(pos,SfAl,label='$S_{f}^{Al}$')
        ax.plot(pos,SfFe,label='$S_{f}^{Fe}$')
        ax.legend() 
        ax.set_xlabel('Transverse position (mm)', fontdict={'size':14})
        ax.set_ylabel('Interface $S_{f}$', fontdict={'size':14})
                
        fig, ax = plt.subplots()
        ax.plot(pos,C,label='C')
        ax.legend() 
     
     
    boxprops = dict(linestyle='-', linewidth=1, color='black')
    medianprops = dict(linestyle='none', linewidth=1, color='black')
    meanlineprops = dict(marker='s', linewidth=1, markerfacecolor='None', 
                         markeredgecolor='firebrick',markersize=3)
    whiskerprops = dict(linestyle='-', linewidth=1, color='black')
    boxwidth = dat.loc[(dat.loc[:,'image_label'] == 1), 'local_position'].max()*echelle/1E3
    
             
    dense_thickn = dat.loc[:,['image_label','dense_thickness','global_position']]
    dense_thickn['dense_thickness'] = dense_thickn['dense_thickness']*echelle
    dense_thickn['global_position'] = dense_thickn['global_position']*echelle/1E3
    dense_thickn.loc[dense_thickn.loc[:,'dense_thickness']==0,'dense_thickness'] = np.nan
    
    
    fig, ax =plt.subplots()
    
    
    
    box = dense_thickn.boxplot(column='dense_thickness',by='image_label',widths=boxwidth, 
                      meanline=False, showmeans=True, boxprops=boxprops,
                      medianprops=medianprops, meanprops=meanlineprops,
                      whis=(0,100),
                      whiskerprops=whiskerprops,return_type='dict', ax = ax,
                      positions=pos,grid=False)
    
    ax.xaxis.set_major_locator(FixedLocator(np.linspace(dat.loc[:,'global_position'].min()*echelle,dat.loc[:,'global_position'].max()*echelle/1E3,7)))
    ax.xaxis.set_major_formatter(FixedFormatter(np.linspace(dat.loc[:,'global_position'].min()*echelle,dat.loc[:,'global_position'].max()*echelle/1E3,7)))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
    ax.set_xlabel('Transverse length (mm)', fontdict={'size':14})
    ax.set_ylabel('Intermetallic thickness (µm)',fontdict={'size':14})
    ax.tick_params(axis='x',labelsize=14, direction='in')
    ax.tick_params(axis='y',labelsize=14, direction='in')
    
    if not plot:
        plt.close()
    
    
    positionals = box.iloc[0]
    
    mu = [item.get_ydata()[0] for item in positionals['means']]
    Xmin, Xmax = ([item.get_ydata()[1] for item in positionals['whiskers'][0::2]],
                               [item.get_ydata()[1] for item in positionals['whiskers'][1::2]])
    Q1, Q3 = ([item.get_ydata()[0] for item in positionals['boxes']],
              [item.get_ydata()[2] for item in positionals['boxes']])
    
    
    xnew = np.linspace(is_imc.loc[:,'global_position'].min(), is_imc.loc[:,'global_position'].max(), 300) 
    
    NaN = []
    for i in range(len(mu)):
        if np.isnan(mu[i]) :
            NaN.append(i)
    
    cleared_mu = np.array([mu[i] for i in range(len(mu)) if i not in NaN])
    cleared_pos = np.array([pos[i] for i in range(len(pos)) if i not in NaN])
    
    to_smooth = pd.DataFrame({'pos':cleared_pos,'mu':cleared_mu})
    to_smooth = to_smooth.sort_values('pos',ascending=True)
    
    
    # BSpline 
    # spl = make_interp_spline(to_smooth['pos'], to_smooth['xb'], k=3)  # type: BSpline
    # power_smooth = spl(xnew)
    # ax.plot(xnew, power_smooth)
    
    # Polynomial
    # coefs = np.polyfit(to_smooth['pos'], to_smooth['xb'], 6) # first and last picture ignored
    # p = np.poly1d(coefs)
    # ax.plot(xnew,p(xnew))
    
    
    # Spline
    
    spl = make_smoothing_spline(to_smooth['pos'], to_smooth['mu'], lam=0.5)
    smoothed = spl(xnew)
    
    if plot:
        ax.plot(xnew,smoothed,label=r'$\tilde{\mu}$')
        ax.legend()
    
    W = smoothed.max()
    
    D = []
    for i in range(len(xnew)-1):
        absolute_der = abs((smoothed[i+1]-smoothed[i])/(xnew[i+1]-xnew[i]))
        D.append(absolute_der)
    
    D = np.mean(D)
    
    L = abs(is_IMC_pos[0] - is_IMC_pos[-1])
    
    
    return {'RaAl':np.array(RaAl), 'RaFe':np.array(RaFe), 'smAl':np.array(smAl), 'smFe':np.array(smFe), 'SfAl': np.array(SfAl), 
            'SfFe':np.array(SfFe), 'C': np.array(C), 'MAD':np.array(MAD), 'pos': np.array(pos), 'mu': np.array(mu),
            'Q1':np.array(Q1),'Q3':np.array(Q3),'Xmin':np.array(Xmin),'Xmax':np.array(Xmax),'W':W,'D':D,'L':L,'delta':disc}


def smoothing(pos,loc_descriptor,lam):
    xnew = np.linspace(pos.min(), pos.max(), 300) 

    NaN = []
    for i in range(len(loc_descriptor)):
        if np.isnan(loc_descriptor[i]) :
            NaN.append(i)

    cleared_loc = np.array([loc_descriptor[i] for i in range(len(loc_descriptor)) if i not in NaN])
    cleared_pos = np.array([pos[i] for i in range(len(pos)) if i not in NaN])

    to_smooth = pd.DataFrame({'pos':cleared_pos,'loc':cleared_loc})
    to_smooth = to_smooth.sort_values('pos',ascending=True)
        
    spl = make_smoothing_spline(to_smooth['pos'], to_smooth['loc'], lam=lam)
    smoothed = spl(xnew)
    
    return xnew, smoothed
