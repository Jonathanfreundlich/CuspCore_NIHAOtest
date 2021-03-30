###############################################################################
#                            FITTING PROFILES                                 #
###############################################################################

import sys
from numpy import *
from lmfit import *
from bisect import bisect
from scipy.integrate import simps

import profiles as prf

###############################################################################

# FIT PROFILES

def do_fits(gl,rvir_fangzhou,mvir_fangzhou,rmax_fit,rmin_fit,m_constraint,components=['d'],irange=[]):

    for (i,ss) in zip(range(size(gl)),gl):
            for c in components:
                a = array(ss['a'])
                r = array(ss[c]['r'])
                brho = array(ss[c]['brho'])
                
                if size(r)>0:
                    Rvir=rvir_fangzhou[i]
                    Mvir=mvir_fangzhou[i]
                
                    #weights array
                    maxr = rmax_fit*Rvir
                    minr = rmin_fit*Rvir
                    wcore = 0*(r<=minr)+1*((r>minr)&(r<=maxr))+0*(r>maxr)

                    #fit brho: b=2, g=3, constrained alpha > -3 sqrt(rmin/rs)
                    params = Parameters()
                    params.add('c', value= 10, min=1.e-16, vary=True)
                    params.add('a', value= 1,vary=True)
                    params.add('b', value= 0.5, vary=False)
                    params.add('g', value= 3, vary=False)
                    params.add('Rvir', value= float(Rvir), vary=False)
                    if c=='d':
                        params.add('Mvir', value= float(Mvir), vary=False)
                    else:
                        params.add('Mvir', value= float(Mvir), vary=True)
                    rmin=r.min()
                    exec("constraint=({'type': 'ineq', 'fun': lambda x:  x[1]+3.5*sqrt(x[0]*10**(-%.12f))})"%(float(m_constraint))) # if the condition applies to s(rmin)
                    #exec("constraint=({'type': 'ineq', 'fun': lambda x:  x[1]+3.*sqrt(x[0]*10**(-%.12f))})"%(float(m_constraint))) # if the condition applies to bs(rmin)
                    fit_output=fit_prof_constrained(r, brho, params, constraint, wcore, model='an', y='brho')
                    ss[c].update({'lsfit_brho_b2_g3_constrained':fit_output})

                    #fit brho: b=2, g=3, unconstrained
                    params = Parameters()
                    params.add('c', value= 10, min=0, vary=True)
                    params.add('a', value= 1,vary=True)
                    params.add('b', value= 0.5, vary=False)
                    params.add('g', value= 3, vary=False)
                    params.add('Rvir', value= float(Rvir), vary=False)
                    if c=='d':
                        params.add('Mvir', value= float(Mvir), vary=False)
                    else:
                        params.add('Mvir', value= float(Mvir), vary=True)
                    fit_output=fit_prof_ls(r, brho, params, wcore, model='an', y='brho')
                    ss[c].update({'lsfit_brho_b2_g3_unconstrained':fit_output})
                else: 
                    fit_output=fit_prof_empty()
                    ss[c].update({'lsfit_brho_b2_g3_constrained':fit_output})
                    ss[c].update({'lsfit_brho_b2_g3_unconstrained':fit_output})
                    
    return gl

def fit_prof_ls(r, data, params, w=1, model='an', y='brho', xi1=0.015, xi2=1.000):
    #least squares fitting
    result = minimize(min_res, params, args=(r, data, w, model, y))
    p = result2p(result, model)
    Rvir = p[-2]
    Mvir = p[-1]
    brho = prf.brho(r, p, model)
    rho = prf.rho(r, p, model)
    bs = prf.bs(r, p, model)
    s = prf.s(r, p, model)
    M = prf.M(r, p, model)
    V = prf.V(r, p, model)
    if y == 'brho':
        rms = rmslog(brho, data)
    if y == 'rho':
        rms = rmslog(rho, data)
    try: c2 = Rvir/r[bisect(bs,2)-1] 
    except: c2 = 0
    return {'p':p,
            'brho':brho,
            'rho':rho,
            'bs':bs,
            's':s,
            'M':M,
            'V':V,
            's1':prf.s(xi1*Rvir, p, model),
            's2':prf.s(xi2*Rvir, p, model),
            'bs1':prf.bs(xi1*Rvir, p, model),
            'bs2':prf.bs(xi2*Rvir, p, model),
            'c2':c2,
            'rms':rms,
            'Merr':(prf.calc_total_mass(r, rho, Rvir)-Mvir)/Mvir}

def fit_prof_constrained(r, data, params, constraint=(), w=1, model='an', y='brho', xi1=0.015, xi2=1.000):
    #least squares fitting
    result = minimize(min_res, params, args=(r, data, w, model, y),method='SLSQP',constraints=constraint)
    p = result2p(result, model)
    Rvir = p[-2]
    Mvir = p[-1]
    brho = prf.brho(r, p, model)
    rho = prf.rho(r, p, model)
    bs = prf.bs(r, p, model)
    s = prf.s(r, p, model)
    M = prf.M(r, p, model)
    V = prf.V(r, p, model)
    if y == 'brho':
        rms = rmslog(brho, data)
    if y == 'rho':
        rms = rmslog(rho, data)
    try: c2 = Rvir/r[bisect(bs,2)-1] 
    except: c2 = 0
    return {'p':p,
            'brho':brho,
            'rho':rho,
            'bs':bs,
            's':s,
            'M':M,
            'V':V,
            's1':prf.s(xi1*Rvir, p, model),
            's2':prf.s(xi2*Rvir, p, model),
            'bs1':prf.bs(xi1*Rvir, p, model),
            'bs2':prf.bs(xi2*Rvir, p, model),
            'c2':c2,
            'rms':rms,
            'Merr':(prf.calc_total_mass(r, rho, Rvir)-Mvir)/Mvir}

def fit_prof_empty():
     return {'p':(nan,nan,nan,nan,nan,nan),
            'brho':[],
            'rho':[],
            'bs':[],
            's':[],
            'M':[],
            'V':[],
            's1':[],
            's2':[],
            'bs1':[],
            'bs2':[],
            'c2':[],
            'rms':[],
            'Merr':[]}
    
###############################################################################

def result2p(result, model):
#converts a fit result object to a tuple of the parameters' values
    if model == 'an':
        c = result.values['c']
        a = result.values['a']
        b = result.values['b']
        g = result.values['g']
        Rvir = result.values['Rvir']
        Mvir = result.values['Mvir']
        p = (c, a, b, g, Rvir, Mvir)
        
    if model == 'dbl':
        f = result.values['f']
        c1 = result.values['c1']
        c2 = result.values['c2']
        a = result.values['a']
        b = result.values['b']
        g = result.values['g']
        Rvir = result.values['Rvir']
        Mvir = result.values['Mvir']
        p = (f, c1, c2, a, b, g, Rvir, Mvir)   
        
    if model == 'nfw':
        c = result.values['c']
        Rvir = result.values['Rvir']
        Mvir = result.values['Mvir']
        p = (c, Rvir, Mvir) 
    
    if model == 'enfw':
        c = result.values['c']
        a = result.values['a']
        b = result.values['b']
        g = result.values['g']
        Rvir = result.values['Rvir']
        Mvir = result.values['Mvir']
        p = (c, a, b, g, Rvir, Mvir)
    
    if model == 'ein':
        c = result.values['c']
        n = result.values['n']
        Rvir = result.values['Rvir']
        Mvir = result.values['Mvir']
        p = (c, n, Rvir, Mvir)
    
    return p

def params2p(params, model):
#converts a Parameters object to a tuple of the parameters' values
    if model == 'an':
        c = params['c'].value
        a = params['a'].value
        b = params['b'].value
        g = params['g'].value
        Rvir = params['Rvir'].value
        Mvir = params['Mvir'].value
        p = (c, a, b, g, Rvir, Mvir)
        
    if model == 'dbl':
        f = params['f'].value
        c1 = params['c1'].value
        c2 = params['c2'].value
        a = params['a'].value
        b = params['b'].value
        g = params['g'].value
        Rvir = params['Rvir'].value
        Mvir = params['Mvir'].value 
        p = (f, c1, c2, a, b, g, Rvir, Mvir)
        
    if model == 'nfw':
        c = params['c'].value
        Rvir = params['Rvir'].value
        Mvir = params['Mvir'].value
        p = (c, Rvir, Mvir)
        
    if model == 'enfw':
        c = params['c'].value
        a = params['a'].value
        b = params['b'].value
        g = params['g'].value
        Rvir = params['Rvir'].value
        Mvir = params['Mvir'].value
        p = (c, a, b, g, Rvir, Mvir)
        
    if model == 'ein':
        c = params['c'].value
        n = params['n'].value
        Rvir = params['Rvir'].value
        Mvir = params['Mvir'].value
        p = (c, n, Rvir, Mvir)     
        
    return p

def Delta_area(logr,logdata1,logdata2,w=1,xlimits=[-2,-1.5],ymin=2.):
    # Warning: all entries must have the same size
    index=where((logr>xlimits[0]) & (logr<xlimits[1]))
    area1=simps(logdata1[index]-ymin,logr[index])
    area2=simps(logdata2[index]-ymin,logr[index])
    Darea=simps(abs(logdata2[index]-logdata1[index]),logr[index])
    return Darea/mean([area1,area2])

def Darea(logr,logdata1,logdata2,w=1,xlimits=[-2,-1.5]):
    # Warning: all entries must have the same size
    index=where((logr>xlimits[0]) & (logr<xlimits[1]))
    Darea=simps(logdata2[index]-logdata1[index],logr[index])
    return Darea
     
###############################################################################

# AUXILIARY FUNCTIONS

def rmslog(fit, data,w=1):
    res = log10(data)-log10(fit)
    res = res*w
    for i in range(len(res)):
        if isnan(res[i]) or isinf(res[i]): res[i] = 0
    return sqrt(mean(res**2))

def min_res(params, r, data, w=1, model='an', y='brho'):
    p = params2p(params, model)  
        
    if y == 'brho':    
        fun = prf.brho(r, p, model)
    if y == 'rho':    
        fun = prf.rho(r, p, model)
        
    return w*(log10(fun) - log10(data))
