###############################################################################
#                       FUNCTIONS TO PREPARE THE DATA                         #
###############################################################################

import sys
import format_functions
import general_functions
from general_functions import *
from numpy import *
from scipy.signal import savgol_filter
from matplotlib.pylab import *

rcParams['figure.figsize'] = (10,6)
rcParams['font.size'] = 18

components=['all', 'd', 's', 'g']

###############################################################################

# DEFINING brho and the slopes s=alpha and sb

def define_brho(gl,polyorder=3,sigma = 11,mode= 'interp',double_smooth=False,rlim=[-2.,0.],use_fangzhou_Rvir=True,verbose=False,D200=False):
    for ss in gl: 
        a = array(ss['a'])
        sim = ss['sim']
        if use_fangzhou_Rvir:
            Rvir=get_fangzhou_radii(sim,array([a]),get_all=False,D200=D200)[2][0]
        else:
            Rvir = ss['Rvir']

        r = array(ss['all']['r'])
        r_range=where((log10(r/Rvir)>=rlim[0])&(log10(r/Rvir)<rlim[1]))
        rmin=r_range[0]
        rmax=r_range[-1]
        
        for c in components:
            n = array(ss[c]['n'])
            M = array(ss[c]['M'])
            dM = M/sqrt(cumsum(n))
            rho=array(ss[c]['rho'])
            drho=rho/sqrt(n)
            
            # Mean density profile
            brho = M/(4.*pi/3.*r**3)
            dbrho = brho/sqrt(cumsum(n))
            
            # Slopes of the densities
            logrho_smooth=nan*ones(size(r))
            s  = nan*ones(size(r))            
            try:
                logrho_smooth[r_range]= savgol_filter(interpolate_nan(log10(rho[r_range])),sigma,polyorder,deriv=0,mode=mode,delta=diff(log10(r))[0])
                if double_smooth:
                    logrho_smooth[r_range]= savgol_filter(logrho_smooth[r_range],sigma,polyorder,deriv=0,mode=mode,delta=diff(log10(r))[0])
                s[r_range] = -savgol_filter(logrho_smooth[r_range],sigma,polyorder,deriv=1,mode=mode,delta=diff(log10(r))[0])
            except:
                if verbose:
                    print 'Warning: logrho_smooth and s could not be defined'            

            logbrho_smooth=nan*ones(size(r))
            bs = nan*ones(size(r))
            try:
                logbrho_smooth[r_range]= savgol_filter(interpolate_nan(log10(brho[r_range])),sigma,polyorder,deriv=0,mode=mode,delta=diff(log10(r))[0])
                if double_smooth:
                    logbrho_smooth[r_range]= savgol_filter(logbrho_smooth[r_range],sigma,polyorder,deriv=0,mode=mode,delta=diff(log10(r))[0])
                bs[r_range] = -savgol_filter(logbrho_smooth[r_range],sigma,polyorder,deriv=1,mode=mode,delta=diff(log10(r))[0])
            except:
                if verbose:
                    print 'Warning: logbrho_smooth and bs could not be defined'   
            
            # Other quantities
            Rmax = r[-1]
            Mmax = M[-1]
            Rvir = ss['all']['Rvir']
            Mvir = ss[c]['Mvir']
            Rthr = nan
            
            ss[c].update({'dM':dM,
                          'drho':drho, 
                          'brho':brho,
                          'dbrho':dbrho,
                          's':s, 
                          'bs':bs,
                          'Rmax':Rmax, 
                          'Mmax':Mmax, 
                          'Rthr':Rthr})
            
    return gl
   
# REDUCE THE RANGE OF THE DIFFERENT ARRAYS INSIDE gl

def get_fitrange(gl,use_fangzhou_Rvir=True,component='d',D200=False):
    fitrange={'all':[],'d':[],'s':[],'g':[]}

    for c in components:
        fitrange_c=[]
        for ss in gl: 
            a = array(ss['a'])
            sim = ss['sim']
            if use_fangzhou_Rvir:
                Rvir=get_fangzhou_radii(sim,array([a]),get_all=False,D200=D200)[2][0]
                Rstar=get_fangzhou_radii(sim,array([a]),get_all=False,D200=D200)[1][0]
                if isnan(Rvir):
                    Rvir = ss['Rvir']
                if isnan(Rstar):
                    Rstar = 0.15*ss['Rvir']
            else:
                Rvir = ss['Rvir']
                Rstar = 0.15*ss['Rvir']

            r = array(ss['all']['r'])
            eps=ss[c]['eps']
            if component=='s':
                outer = r <= 0.15*Rvir
            else:
                outer = r <= Rvir
            conv  = r > 0.01*Rvir    
            soft  = r >= eps
            fitrange_c.append(conv & soft & outer)
            
        fitrange[c]=fitrange_c
        
    return fitrange

def reduce_range_gl(gl,fitrange,verbose=False):
    print 'Reducing the range of gl'
    for (i,ss) in zip(range(size(gl)),gl): 
                
        for c in components:
            sizei=size(ss[c]['r'])
            fitrange_i=fitrange[c][i]
            for key in ss[c]:
                if size(ss[c][key])==size(fitrange_i):
                    ss[c][key]=ss[c][key][fitrange_i]
            for key in ss['flowdata'][c]['in']:
                if size(ss['flowdata'][c]['in'][key])==size(fitrange_i):
                    ss['flowdata'][c]['in'][key]=ss['flowdata'][c]['in'][key][fitrange_i]
                    ss['flowdata'][c]['out'][key]=ss['flowdata'][c]['out'][key][fitrange_i]
            sizef=size(ss[c]['r'])
            if verbose:
                if c=='d':
                    print 'i=%i (%s): from %i to %i'%(i,c,sizei,sizef)
    return gl

def reduce_range_Treal(Treal,fitrange,c='d',verbose=False):
    print 'Reducing the range of Treal'
    Treal_b=[]
    for i in range(shape(Treal)[0]):
        [r,Rvir,Tr]=Treal[i]
        sizei=size(r)
        fitrange_i=fitrange[c][i]
        rb=r[fitrange_i]
        Trb=Tr[fitrange_i]
        sizef=size(rb)
        Treal_b.append([rb,Rvir,Trb])
        if verbose:
            print 'i=%i (%s): from %i to %i'%(i,c,sizei,sizef)
    return Treal_b
        
def plot_quantities(gl,k,c):
    ssc=gl[k][c]
    r=ssc['r']
    Rvir=ssc['Rvir']
    x=log10(r/Rvir)
    for key in ssc:
        figure()
        plot(x,ssc[key]*ones(size(x)),'r')
        xlabel(r'$\log(r/R_{vir})$')
        ylabel(r'%s'%key)
        show()
        
    
