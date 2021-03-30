###############################################################################
#                         MODEL: EVOLVE PROFILES                              #
###############################################################################

import sys
from numpy import *
from lmfit import *
from scipy.signal import savgol_filter
import profiles as prf
import fitting as fit

# gravitational constant [kpc^3 Gyr^-2 Msun^-1]
G = 4.499753324353496e-06 

###############################################################################

'''
methods:
- halo: the main toy model's method
- halo_noUex: halo method without potential energy from outer shells
- shell: lagrangian shell method
- shell_halo: paricle inside a fixed halo
- adiabatic: adiabatic flow model
    
difference between r and ri:
r - fixed range of radii for halo, halo_noUex methods
ri - radii that change with the evolution of lagrangian shells
    

'''

def evolve(r, ri, Mi, pi, m, alphai=0.,alphaf=0.,gammai=0., gammaf=0.,betai=0.,betaf=0.,Ttype='jeans',Mcen=0, w=1, model='an', method='halo',add_params=[],Rvirf=0.,Mvirf=0.,component='d'):

    gamma=0.
    
    if method == 'halo' or method == 'halo_noUex':
        if model == 'an':
            (ci, ai, bi, gi, Rviri, Mviri) = pi
            if Rvirf==0.: Rvirf=Rviri
            if Mvirf==0.: Mvirf=Mviri 
            
            paramsf = Parameters()
            paramsf.add('c', value= ci, min=1e-16, vary=True)
            paramsf.add('a', value= ai, vary=True)
            paramsf.add('b', value= bi, vary=False)
            paramsf.add('g', value= gi, vary=False)
            paramsf.add('Rvir', value= float(Rvirf), vary=False)
            if component=='d':
                        paramsf.add('Mvir', value= float(Mvirf), vary=False)
            else:
                        paramsf.add('Mvir', value= float(Mvirf), vary=True)
            if Ttype=='jeans-gamma':
                paramsf.add('gamma', value= float(gamma), min=-1., max=1.,vary=True)
        result = minimize(min_E_diff, paramsf, args=(r, pi, m,alphai,alphaf,gammai, gammaf,betai,betaf,Ttype,w, model, method,add_params))
        if Ttype=='jeans-gamma':
            gamma=result.values['gamma']
        pf = fit.result2p(result, model)

        #energy error calc
        ere = E_diff(r, pi, pf, m,alphai,alphaf,gammai, gammaf,betai,betaf,Ttype,model, method,add_params)[1]
        errterms = [log10(r/Rviri) < -1.33, #inner region
                       (log10(r/Rviri) >= -1.33) & (log10(r/Rviri) <= -0.67), #middle region
                       log10(r/Rviri) > -0.67, #outer region
                       r == r] #all 
        ererms = [sqrt(mean(ere[errterm]**2)) for errterm in errterms]
        
        #cosider mass bath at center - important for a sequence
        Mcenf = Mcen + m
        rf = r
        rf0 = prf.inv_M(Mi, pf)
        Mf = prf.M(rf, pf, model)
    
    return {'ri':ri,#initial radii
            'rf0':rf0, #final radii
            'rf':rf, #final radii, sorted
            'Mf':Mf, #final mass profile
            'Mcenf':Mcenf, #final mass at center
            'pi':pi, #initial profile parameters
            'pf':pf, #final profile parameters
            'gammaf':gamma, # gamma
            'm':m, #added mass
            'model':model, #profile model used
            'method':method, #evolution method used
            'ere':ere, #energy relative error array
            'ererms':ererms} #rms of energy relative error    

def evolve_constrained(r, ri, Mi, pi, m, constraint=(),alphai=0.,alphaf=0.,gammai=0., gammaf=0.,betai=0.,betaf=0.,Ttype='jeans',Mcen=0, w=1, model='an', method='halo',add_params=[],Rvirf=0.,Mvirf=0.):

    gamma=0.
    print Rvirf, Mvirf
    
    if method == 'halo' or method == 'halo_noUex':
        if model == 'an':
            (ci, ai, bi, gi, Rviri, Mviri) = pi
            if Rvirf==0.: Rvirf=Rviri
            if Mvirf==0.: Mvirf=Mviri 
            
            paramsf = Parameters()
            paramsf.add('c', value= ci, min=1e-16, vary=True)
            paramsf.add('a', value= ai, vary=True)
            paramsf.add('b', value= bi, vary=False)
            paramsf.add('g', value= gi, vary=False)
            paramsf.add('Rvir', value= float(Rvirf), vary=False)
            paramsf.add('Mvir', value= float(Mvirf+Mcen), vary=False)
            print paramsf
            if Ttype=='jeans-gamma':
                paramsf.add('gamma', value= float(gamma), min=-1., max=1.,vary=True)
        
        result = minimize(min_E_diff, paramsf, args=(r, pi, m,alphai,alphaf,gammai, gammaf,betai,betaf,Ttype,w, model, method,add_params),method='SLSQP',constraints=constraint,options={'disp': True})
        if Ttype=='jeans-gamma':
            gamma=result.values['gamma']
        pf = fit.result2p(result, model)
        
        #energy error calc
        ere = E_diff(r, pi, pf, m,alphai,alphaf,gammai, gammaf,betai,betaf,Ttype,model, method,add_params)[1]
        errterms = [log10(r/Rviri) < -1.33, #inner region
                       (log10(r/Rviri) >= -1.33) & (log10(r/Rviri) <= -0.67), #middle region
                       log10(r/Rviri) > -0.67, #outer region
                       r == r] #all 
        ererms = [sqrt(mean(ere[errterm]**2)) for errterm in errterms]
        
        #consider mass bath at center - important for a sequence
        Mcenf = Mcen + m
        rf = r
        rf0 = prf.inv_M(Mi, pf)
        Mf = prf.M(rf, pf, model)
    
    return {'ri':ri,#initial radii
            'rf0':rf0, #final radii
            'rf':rf, #final radii, sorted
            'Mf':Mf, #final mass profile
            'Mcenf':Mcenf, #final mass at center
            'pi':pi, #initial profile parameters
            'pf':pf, #final profile parameters
            'gammaf':gamma, # gamma
            'm':m, #added mass
            'model':model, #profile model used
            'method':method, #evolution method used
            'ere':ere, #energy relative error array
            'ererms':ererms} #rms of energy relative error  

###############################################################################

# AUXILIARY FUNCTIONS

def min_E_diff(paramsf, ri, pi, m,alphai=0,alphaf=0,gammai=0., gammaf=0.,betai=0.,betaf=0.,Ttype='jeans', w=1, model='an', method='halo',add_params=[]): 
    #minization function for the energy difference
    pf = fit.params2p(paramsf, model)
    
    return w*E_diff(ri, pi, pf, m, alphai,alphaf,gammai, gammaf,betai,betaf,Ttype,model, method,add_params)[0]

def E_diff(ri, pi, pf, m, alphai=0., alphaf=0.,gammai=0., gammaf=0.,betai=0.,betaf=0., Ttype='jeans',model='an', method='halo',add_params=[],do_smooth=False):
    #returns the energy difference of the mass enclosing shell between the before and after states. the difference should be ideally zero according to the model
    
    if Ttype[-5:]=='Mreal':
        Mi = array(add_params)
    else:
        Mi = prf.M(ri, pi,model)
    
    rf = prf.inv_M(Mi, pf, model) #the new radii that enclose the same masses
    Ui = prf.U(ri, pi, model)
    Uf = prf.U(rf, pf, model)
    
    if Ttype=='gamma-Treal':
        Treal=add_params
        gammai=1.5*G*Mi/ri/Treal-prf.alpha_Dekel(ri,pi)
        gammaf=gammai
    
    Ti=get_T(ri,[alphai,betai,gammai,pi],m=0.,add_params=add_params,Ttype=Ttype,do_smooth=do_smooth)
    Tf=get_T(rf,[alphaf,betaf,gammaf,pf],m=m ,add_params=add_params,Ttype=Ttype,do_smooth=do_smooth)
    
    Ei = Ui - G*m/ri + Ti #0.5*G*Mi/ri
    Ef = Uf - G*m/rf + Tf #0.5*G*Mi/rf    

    return array([Ef-Ei, (Ef-Ei)/abs(Ei)],dtype=float)

def get_T(r,params,m=0.,add_params=[],Ttype='jeans-alpha',do_smooth=False,polyorder=3,sigma = 21,mode= 'interp',rlim=[-2,0]):
    alpha, beta, gamma, p = params
    if Ttype[-5:]=='Mreal':
        M = array(add_params)+array(m)
    else:
        M = prf.M(r, p,model='an')+m
    
    Rvir=p[-2]

    if Ttype=='jeans-alpha' or Ttype=='jeans' or Ttype=='jeans-Mreal' or Ttype=='alpha-p-Mreal':
        denominator=get_denominator(r,params,Ttype=Ttype)
        T=0.5*(3.-2*beta)/denominator*G*M/r
    elif Ttype=='alpha' or Ttype=='alpha-Mreal' or Ttype=='betazero' or Ttype=='gamma-Treal':
        denominator=get_denominator(r,params,Ttype=Ttype)
        T=1.5/denominator*G*M/r
    elif Ttype=='zero':
        T=zeros_like(r)
    elif Ttype=='virial' or Ttype=='virial-Mreal':
        T=0.5*G*M/r
    elif Ttype=='jeans-smooth':
        T=0.5*(3.-2*beta)/denominator*G*M/r
        do_smooth=True
    elif Ttype=='Tdekel':
        # T from hydrostatic equilibrum (sigmar)
        (c, a, b, g, Rvir, Mvir) = p
        x=r/Rvir*c
        sigmar2=prf.sigmar2_dekel_m(x,Mvir,Rvir,c,a,m=m,mtype='center')
        T=1.5*sigmar2
    elif Ttype=='Tmulti':
        (c, a, b, g, Rvir, Mvir) = p
        [Mratio,n]=add_params
        x=r/Rvir*c
        T=prf.K_Mratio(x,Mvir,Rvir,c,a,Mratio,n,m)
            
    if do_smooth:
        r_range=where((log10(r/Rvir)>=rlim[0])&(log10(r/Rvir)<rlim[1]))
        T_smooth=nan*ones(size(T))
        T_smooth[r_range]= savgol_filter(T[r_range],sigma,polyorder,deriv=0,mode=mode,delta=diff(log10(r))[0])
        T=T_smooth
    
    return T

def get_denominator(r,params,Ttype='jeans-alpha'):
    alpha, beta, gamma, p = params
    
    if Ttype=='jeans-alpha' or Ttype=='alpha-p-Mreal':
        denominator=prf.alpha_Dekel(r,p)+gamma-2*beta
    elif Ttype=='alpha' or Ttype=='alpha-Mreal':
        denominator=prf.alpha_Dekel(r,p)
    elif Ttype=='betazero' or Ttype=='gamma-Treal':
        denominator=prf.alpha_Dekel(r,p)+gamma   
    elif Ttype=='Tdekel':
        denominator=1.
    else:
        denominator=alpha+gamma-2*beta
   
    denominator=redress_denominator(denominator)
    return denominator
        
def redress_denominator(denominator):
    # Prevent denominator=alpha+gamma-2beta to be zero or negative
    if size(where(denominator<=0)[0])>0:
        imin=where(denominator<=0)[0][-1]+1
        denominator[:imin]=denominator[imin]*ones(size(denominator[:imin]))
    return denominator
    