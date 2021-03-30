###############################################################################
#          FUNCTIONS PERTAINING TO THE SLOPES ALPHA, BETA, GAMMA              #
###############################################################################

import sys
from numpy import *
from scipy.signal import savgol_filter
from matplotlib.pylab import *

import general_functions
from general_functions import *

components=['all', 'd', 'g','s']

###############################################################################

def derive_slopes(gl,polyorder=3,sigma = 11,mode= 'interp',double_smooth=False,rlim=[-2.,0.],use_fangzhou_Rvir=True,D200=False,linearize=False,betanull=False,verbose=False):
    '''
    Warning: When using savgol_filter, delta=diff(log10(r))[0] means that the 
    radius r has to be spaced logarithmically! 
    '''
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
        x=log10(r/Rvir)
        
        for c in components:   
            n = array(ss[c]['n'])
            M = array(ss[c]['M'])
            Mall=array(ss['all']['M'])
            dM = M/sqrt(cumsum(n))  
            
            # DEFINE ALPHA
            rho=array(ss[c]['rho'])
            logrho=log10(rho)
            logrho_smooth=nan*ones(size(logrho))
            alpha=nan*ones(size(logrho))
            try:
                logrho_smooth[r_range]= savgol_filter(interpolate_nan(logrho[r_range]),sigma,polyorder,deriv=0,mode=mode,delta=diff(log10(r))[0])
                if double_smooth:
                    logrho_smooth[r_range]= savgol_filter(logrho_smooth[r_range],sigma,polyorder,deriv=0,mode=mode,delta=diff(log10(r))[0])
                alpha[r_range] = -savgol_filter(logrho_smooth[r_range],sigma,polyorder,deriv=1,mode=mode,delta=diff(log10(r))[0])
                if linearize:
                    p=polyfit(x[r_range],alpha[r_range],1)
                    alpha[r_range]=p[0]*x[r_range]+p[1]
            except:
                if verbose:
                    print 'Warning: logrho_smooth and alpha could not be defined'            
                    
            # DEFINE BETA
            beta = noinf(array(ss[c]['beta']))
            beta_smooth=nan*ones(size(beta))
            try:
                beta_smooth[r_range]= savgol_filter(interpolate_nan(beta[r_range]),sigma,polyorder,deriv=0,mode=mode,delta=diff(log10(r))[0])
                if double_smooth:
                    beta_smooth[r_range]= savgol_filter(beta_smooth[r_range],sigma,polyorder,deriv=0,mode=mode,delta=diff(log10(r))[0])
                if linearize:
                    p=polyfit(x[r_range],beta_smooth[r_range],1)
                    beta_smooth[r_range]=p[0]*x[r_range]+p[1]
                if betanull:
                    beta_smooth[r_range]=0.*x[r_range]
            except:
                if verbose:
                    print 'Warning: beta_smooth could not be defined'

            # DEFINE GAMMA
            sigmar = array(ss[c]['vr_disp'].in_units('km s^-1'))
            logsigmar = log10(sigmar)
            logsigmar2 = 2*log10(sigmar)
            logsigmar_smooth=nan*ones(size(logsigmar))
            logsigmar2_smooth=nan*ones(size(logsigmar2))
            gamma=nan*ones(size(logsigmar2))
            try:
                logsigmar_smooth[r_range]= savgol_filter(interpolate_nan(logsigmar[r_range]),sigma,polyorder,deriv=0,mode=mode,delta=diff(log10(r))[0])
                logsigmar2_smooth[r_range]= savgol_filter(interpolate_nan(logsigmar2[r_range]),sigma,polyorder,deriv=0,mode=mode,delta=diff(log10(r))[0])
                if double_smooth:
                    logsigmar_smooth[r_range]= savgol_filter(logsigmar_smooth[r_range],sigma,polyorder,deriv=0,mode=mode,delta=diff(log10(r))[0])
                    logsigmar2_smooth[r_range]= savgol_filter(logsigmar2_smooth[r_range],sigma,polyorder,deriv=0,mode=mode,delta=diff(log10(r))[0])
                gamma[r_range] = -savgol_filter(logsigmar2_smooth[r_range],sigma,polyorder,deriv=1,mode=mode,delta=diff(log10(r))[0])
                if linearize:
                    p=polyfit(x[r_range],gamma[r_range],1)
                    gamma[r_range]=p[0]*x[r_range]+p[1]
            except:
                if verbose:
                    print 'Warning: sigmar_smooth, sigmar2_smooth and gamma could not be defined'

            ss[c].update({'M':noinf(M),
                          'Mall':noinf(Mall),
                          'rho':10**noinf(logrho),
                          'logrho':noinf(logrho),
                          'rho_smooth':10**noinf(logrho_smooth),
                          'logrho_smooth':noinf(logrho_smooth),
                          'sigmar':noinf(sigmar),
                          'logsigmar':noinf(logsigmar),
                          'logsigmar2':noinf(logsigmar2),
                          'logsigmar_smooth':noinf(logsigmar_smooth),
                          'logsigmar2_smooth':noinf(logsigmar2_smooth),
                          'sigmar_smooth':10**noinf(logsigmar_smooth),
                          'sigmar2_smooth':10**noinf(logsigmar2_smooth),
                          'alpha':noinf(alpha),
                          'beta_smooth':noinf(beta_smooth),
                          'gamma':noinf(gamma)})

    return gl

###############################################################################

def plot_q(gl,k,quantities=['logrho'],c='d',rlim=[-2.,0.]):
    ss=gl[k][c]
    r=ss['r']
    Rvir=gl[k]['Rvir']
    x=log10(r/Rvir)
    r_range=(x>=rlim[0])&(x<rlim[1])
    
    for quantity in quantities:
        figure()
        plot(x[r_range],ss[quantity][r_range],'r')
        if quantity+'_smooth' in ss.keys():
            plot(x[r_range],ss[quantity+'_smooth'][r_range],'k--')
        xlabel(r'$\log(r/R_{vir})$')
        label=q_label(quantity)
        ylabel(r'$%s$'%label)
        ax=gca()
        ax.set_xlim(rlim)  
        
    
def q_label(quantity):
    if quantity=='logrho':
        return '\\log(\\rho)'
    if quantity=='alpha':
        return '\\alpha'
    if quantity=='beta' or quantity=='beta_smooth':
        return '\\beta'
    if quantity=='gamma':
        return '\\gamma'
