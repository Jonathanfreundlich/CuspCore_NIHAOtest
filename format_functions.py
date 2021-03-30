###############################################################################
#      FORMATTING THE PROFILE AND FLOW DATA FROM THE SIMULATION OUTPUTS       #
###############################################################################

import os
import sys
from numpy import *
import pickle
import warnings
import pynbody
import general_functions
from general_functions import *

G = pynbody.array.SimArray(6.67e-11)
G.units = 'm**3 kg**-1 s**-2'
snapshots = ['.' + str(i).zfill(5) for i in range(16, 1024+16, 16)]
 
###############################################################################

# LOAD GL FILE IF IT EXISTS, CREATE IT OTHERWISE
def load_or_create_gl(sim,centering_com='s',directory='/cs/sci/freundlich/CUSPCORE/NIHAO_data/',name='NIHAO',remove=False,use_fangzhou_Rvir=False,get_vcirc=False,D200=False, get_flowdata=True): 
    if remove: 
        remove_files(sim,directory,name)
     
    if os.path.exists(directory+name+'-%s.pickle'%(sim[1:])):
        print ' '
        print 'Loading existing file %s-%s.pickle'%(name,sim[1:])
        gl=load_gl(sim,directory,name)
    
    else:
        print ' ' 
        print 'Creating file %s-%s.pickle...'%(name,sim[1:])
        
        # Load temporary file if it exists
        if os.path.exists(directory+name+'-%s-temp.pickle'%(sim[1:])):
            gl=load_gl_temp(sim,directory,name)
            i_start=size(gl)           
        else:
            gl=[]
            i_start=0
        
        # Load remaining snapshots
        for (ss, i) in zip(snapshots, range(len(snapshots))):
            if i<i_start:
                print sim, 'snap=',ss[1:],'output=',i, '-- already defined'
            else:
                print sim, 'snap=',ss[1:], 'output=',i
                
                data=getdata(sim,ss,centering_com=centering_com,use_fangzhou_Rvir=use_fangzhou_Rvir,get_vcirc=get_vcirc, D200=D200)
                if get_flowdata:
                    if not isnan(data['a']):
                        data=derive_SFRsnapshot(sim,ss,data)
                        data=derive_flowdata(sim,ss,data)
                    
                gl.append(data)
                save_gl_temp(sim,gl,directory,name)
                flush()
    
        save_gl(sim,gl,directory,name)
    return gl

def load_gl(sim,directory,name='NIHAO'):
    with open(directory+name+'-%s.pickle'%(sim[1:])) as f:
        gl = pickle.load(f)
    return gl

def save_gl(sim,gl,directory,remove_file=True,name='NIHAO'):
    if remove_file:
        try:
            os.remove(directory+name+'-%s.pickle'%(sim[1:]))
        except:
            'No file to remove'
    with open(directory+name+'-%s.pickle'%(sim[1:]), 'w') as f:
        pickle.dump(gl, f)

def load_gl_temp(sim,directory,name='NIHAO'):
    with open(directory+name+'-%s-temp.pickle'%(sim[1:])) as f:
        gl = pickle.load(f)
    return gl   

def save_gl_temp(sim,gl,directory,name='NIHAO'):
    with open(directory+name+'-%s-temp.pickle'%(sim[1:]), 'w') as f:
        pickle.dump(gl, f)

def flush():  
    return sys.stdout.flush()

def remove_files(sim,directory,name='NIHAO'):
    try:
        os.remove(directory+name+'-%s-temp.pickle'%(sim[1:]))
        print '%s-%s-temp.pickle was removed'%(name,sim[1:])
    except: 
        print 'No file named %s-%s-temp.pickle'%(name,sim[1:])
    try:
        os.remove(directory+name+'-%s.pickle'%(sim[1:]))
        print '%s-%s.pickle was removed'%(name,sim[1:])
    except: 
        print 'No file named %s-%s.pickle'%(name,sim[1:])   


###############################################################################

# GET PROFILE DATA FROM THE SIMULATION OUTPUTS
def getdata(sim,ss,centering_com='s',nbins=150,use_fangzhou_Rvir=False,get_vcirc=False,D200=False):
  print '   Run getdata...'
  
  try:
    s=load_s(sim,ss,centering_com)
    
    a = s.properties['a']
    z = 1/a - 1
    t = pynbody.array.SimArray(float(s.properties['time'].in_units('Gyr')))
    t.units = 'Gyr'

    # GET GENERAL PARAMETERS
    if use_fangzhou_Rvir:
        print '   Using Rvir from Fangzhou...'
        ok_fangzhou,r12_fangzhou,rvir_fangzhou,mvir_fangzhou=get_fangzhou_radii(sim,array([a]),get_all=False,D200=D200)
        Rvir=rvir_fangzhou[0]
    else:
        h1=s.halos()[1]
        Rvir = float(max(h1['r']))
        
    rho_crit = pynbody.analysis.cosmology.rho_crit(s, unit='Msol kpc^-3')
    fbar = 0.15
    try:
        h1=s.halos()[1]
        fsv = float(sum(h1.s['mass'])/(sum(h1['mass'])*fbar))
        er = float(h1.properties['Rmax']/max(h1['r']))
        Mvir = float(sum(h1['mass'].in_units('Msol')))
    except:
        fsv = nan
        er = nan
        Mvir = nan 

    # MSTAR, SFR AND SSFR
    print '   Calculating SFR...'
    Rmax = 0.15*Rvir
    Mstar_tot = pynbody.array.SimArray(sum(s.s[pynbody.filt.Sphere(Rmax)]['mass']))
    SFR_tot=get_SFR(s,Rmax,age_min=40.,age_max=80.,age_step=0.2)
    SFR_Ha_tot=get_SFR(s,Rmax,age_min=30.,age_max=60.,age_step=0.2)
    SFR_UV_tot=get_SFR(s,Rmax,age_min=80.,age_max=120.,age_step=0.2)
    SSFR_tot=SFR_tot/Mstar_tot*1.e9    # In Gyr-1
    
    # GET BINNED PROFILES
    print '   Get profiles...'
    Rmin = s['eps'].min()
    Rmax = 10*Rvir
    p = pynbody.analysis.profile.Profile(s, min=Rmin, max=Rmax, type='log', nbins=nbins, ndim=3)
    pg = pynbody.analysis.profile.Profile(s.g, min=Rmin, max=Rmax, type='log', nbins=nbins, ndim=3)
    pd = pynbody.analysis.profile.Profile(s.d, min=Rmin, max=Rmax, type='log', nbins=nbins, ndim=3)
    try:
        ps = pynbody.analysis.profile.Profile(s.s, min=Rmin, max=Rmax, type='log', nbins=nbins, ndim=3)
    except:
        ps = p
    
    try:
        seps = min(s.s['eps'])
    except:
        seps = 0
        
    if get_vcirc:
        try:
            vcirc=p['v_circ']
        except: 
            vcirc=nan*ones_like(p['rbins'])
            print 'Warning: vcirc not calculated for output %s.%s'%(sim,ss)
        # If we want to calculate v_circ for all components:
        #for ci in ['','g','d','s']:
        #    exec('pi=p%s'%ci)
        #    try:
        #        exec("vcirc%s=p%s['v_circ']"%(ci,ci))
        #    except:
        #        exec("vcirc%s=nan*ones_like(p%s['rbins'])"%(ci,ci))

    #pack everything in an array
    print '   Pack data...'             
    data = {'sim':sim, 
            'a':a,
            'z':z,
            't':t,
            'Rvir':Rvir,
            'rho_crit':rho_crit, 
            'fsv':fsv, 
            'er':er,
            #'flowdata':flowdata,
            'Mstar':Mstar_tot,
            'SFR':SFR_tot,
            'SFR_Ha':SFR_Ha_tot,
            'SFR_UV':SFR_UV_tot,
            'all':{'r':p['rbins'].in_units('kpc'), 
                   'M':p['mass_enc'].in_units('Msol'), 
                   'n':cumsum(p['n']),
                   'Rvir':Rvir, 
                   'Mvir':Mvir, 
                   'eps':min(s['eps']),
                   'rho':p['density'].in_units('Msol kpc^-3'),
                   'vr_disp':p['vr_disp'].in_units('km s^-1'),
                   'vt_disp':p['vt_disp'].in_units('km s^-1'),
                   'beta':p['beta'],
                   'vphi':p['vphi']},                 
            'd':{  'r':pd['rbins'].in_units('kpc'), 
                   'M':pd['mass_enc'].in_units('Msol'), 
                   'n':cumsum(pd['n']),
                   'Rvir':Rvir, 
                   'Mvir':Mvir, 
                   'eps':min(s.d['eps']),
                   'rho':pd['density'].in_units('Msol kpc^-3'),
                   'vr_disp':pd['vr_disp'].in_units('km s^-1'),
                   'vt_disp':pd['vt_disp'].in_units('km s^-1'),
                   'beta':pd['beta'],
                   'vphi':pd['vphi']},
            'g':{  'r':pg['rbins'].in_units('kpc'), 
                   'M':pg['mass_enc'].in_units('Msol'), 
                   'n':cumsum(pg['n']),
                   'Rvir':Rvir, 
                   'Mvir':Mvir, 
                   'eps':min(s.g['eps']),
                   'rho':pg['density'].in_units('Msol kpc^-3'),
                   'vr_disp':pg['vr_disp'].in_units('km s^-1'),
                   'vt_disp':pg['vt_disp'].in_units('km s^-1'),
                   'beta':pg['beta'],
                   'vphi':pg['vphi']},
            's':{  'r':ps['rbins'].in_units('kpc'), 
                   'M':ps['mass_enc'].in_units('Msol'), 
                   'n':cumsum(ps['n']),
                   'Rvir':Rvir, 
                   'Mvir':Mvir, 
                   'eps':seps,
                   'rho':ps['density'].in_units('Msol kpc^-3'),
                   'vr_disp':ps['vr_disp'].in_units('km s^-1'),
                   'vt_disp':ps['vt_disp'].in_units('km s^-1'),
                   'beta':ps['beta'],
                   'vphi':ps['vphi']}}
   
    if get_vcirc:
        data['all'].update({'vcirc':vcirc})
  
  except:
    print '   The output could not be loaded'
    data = {'sim':sim, 
            'a':nan,
            'z':nan,
            't':nan,
            'Rvir':nan,
            'rho_crit':nan, 
            'fsv':nan, 
            'er':nan,
            #'flowdata':flowdata,
            'Mstar':nan,
            'SFR':nan,
            'SFR_Ha':nan,
            'SFR_UV':nan,
            'all':{'r':nan, 
                   'M':nan, 
                   'n':nan,
                   'Rvir':nan, 
                   'Mvir':nan, 
                   'eps':nan,
                   'rho':nan,
                   'vr_disp':nan,
                   'vt_disp':nan,
                   'beta':nan,
                   'vphi':nan},                 
            'd':{  'r':nan, 
                   'M':nan, 
                   'n':nan,
                   'Rvir':nan, 
                   'Mvir':nan, 
                   'eps':nan,
                   'rho':nan,
                   'vr_disp':nan,
                   'vt_disp':nan,
                   'beta':nan,
                   'vphi':nan},
            'g':{  'r':nan, 
                   'M':nan, 
                   'n':nan,
                   'Rvir':nan, 
                   'Mvir':nan, 
                   'eps':nan,
                   'rho':nan,
                   'vr_disp':nan,
                   'vt_disp':nan,
                   'beta':nan,
                   'vphi':nan},
            's':{  'r':nan, 
                   'M':nan, 
                   'n':nan,
                   'Rvir':nan, 
                   'Mvir':nan, 
                   'eps':nan,
                   'rho':nan,
                   'vr_disp':nan,
                   'vt_disp':nan,
                   'beta':nan,
                   'vphi':nan}}    
  
  return data

def load_s(sim,ss,centering_com='s'):
    print '   --- load_s %s%s'%(sim,ss)
    s = pynbody.load('/vol/sci/astro/cosmo/nas2/Data/nihao/' + sim + '/' + sim + ss)
    s.physical_units()
    
    try: 
        h = s.halos()
        h1 = h[1]
    
        # CENTERING
        if centering_com=='d':
            pynbody.analysis.halo.center(h1.d,mode='ssc',shrink_factor=0.95)
        if centering_com=='s':
            try:
                pynbody.analysis.halo.center(h1.s,mode='ssc',shrink_factor=0.95)
            except:
                pynbody.analysis.halo.center(h1.d,mode='ssc',shrink_factor=0.95)
    
        pynbody.analysis.angmom.faceon(h1)
    except: 
        print '   --- Warning: no centering for output %s.%s'%(sim,ss)
        
    return s        

def get_SFR(self,Rmax,age_min=40.,age_max=80.,age_step=0.2):
    "Determine the SFR from stars of age between 40 and 80 Myr (20-40 for Halpha)"
    "Instead of taking sum(masses[40<ages<80])/40 we take the mean of sum(masses[ages<t])/t"
    "Unit: Msol/yr"
    ages=get_age(self.s[pynbody.filt.Sphere(Rmax)]).in_units('Myr')
    masses=self.s[pynbody.filt.Sphere(Rmax)]['mass'].in_units('Msol')
    age_series = arange(age_min, age_max + age_step/2., age_step)
    sfr_series = zeros((len(age_series)))
    for j, age_j in enumerate(age_series):
        sfr_series[j] = sum( masses[ages<= age_j] / (age_j * 1.e6) )
    SFR=mean(sfr_series)
    return SFR

def get_age(self):
    """Stellar age determined from formation time and current snapshot time"""
    return self.properties['time'].in_units(self['tform'].units, **self.conversion_context()) - self['tform']

###############################################################################

# COUNT THE STELLAR MASS FORMED BETWEEN TWO SNAPSHOTS
# AND THE ASSOCIATED SFR, SSFR

def derive_SFRsnapshot(sim,ss,data,centering_com='s'):
    print '   Run derive_SFRsnapshot...'

    ss_prev=snapshots[snapshots.index(ss)-1]

    s=load_s(sim,ss,centering_com)    
    t = pynbody.array.SimArray(float(s.properties['time'].in_units('Gyr')))
    
    try: 
        t0=get_time(sim,ss_prev)
        h1=s.halos()[1]
        M = sum(h1.s['mass']).in_units('Msol')
        try:
            massform = sum(h1.s['massform'][(h1.s['tform'] > t0) & (h1.s['tform'] <= t)]).in_units('Msol')
        except:
            massform = 0.
        SFR_snapshot = (massform/(t-t0))
        SSFR_snapshot= (SFR_snapshot/M)
    except: 
        massform=nan
        SFR_snapshot=nan
        SSFR_snapshot=nan
        print '   --- Warning: no SFR for output %s.%s'%(sim,ss)
    
    data.update({'DM_snapshot':massform,
                 'SFR_snapshot':SFR_snapshot,
                 'SSFR_snapshot':SSFR_snapshot})
    return data

def get_time(sim,ss):
    s=load_s(sim,ss)
    t = pynbody.array.SimArray(float(s.properties['time'].in_units('Gyr')))
    return t
    
###############################################################################

# GET THE FLOW DATA BETWEEN TWO SNAPSHOTS

def derive_flowdata(sim,ss,data,centering_com='s'):
    print '   Run derive_flowdata...'

    s=load_s(sim,ss,centering_com)
    
    ss_prev=snapshots[snapshots.index(ss)-1]
    
    r = array(data['all']['r'])
    Rvir = array(data['Rvir'])
    
    try:
        sprev=load_s(sim,ss_prev,centering_com)
        if sprev != []:
            flowdata = getflowseq(sprev, s, r, Rvir, sim)
        else:
            flowdata = []
    except:
        flowdata = []
    
    data.update({'flowdata':flowdata})
    return data
    
def getflowseq(s1, s2, r, Rvir, sim, dmo=False):
    print '   --- getflowseq'
    #returns a dictionary with arrays of the inflow and outflow data through multiple spheres with radii given in r
    flowdata = {'all':{'in':{},'out':{}},
              'd':{'in':{},'out':{}},
              'g':{'in':{},'out':{}},
              's':{'in':{},'out':{}}}

    a1 = s1.properties['a']
    z1 = 1/a1 - 1
    t1 = s1.properties['time'].in_units('Gyr')

    a2 = s2.properties['a']
    z2 = 1/a2 - 1
    t2 = s2.properties['time'].in_units('Gyr')
    
    dt = pynbody.array.SimArray(t2-t1)
    dt.units = 'Gyr'

    for rt in r:
        rt = pynbody.array.SimArray(rt)
        rt.units = 'kpc'
        fdt = getflow(s1,s2,rt,Rvir,t1,sim,dt,dmo)
        for comkey in fdt:
            for inoutkey in fdt[comkey]:
                if flowdata[comkey][inoutkey] == {}: 
                    flowdata[comkey][inoutkey] = {key:[fdt[comkey][inoutkey][key]] for key in fdt[comkey][inoutkey]}
                else:
                    for key in fdt[comkey][inoutkey]:
                        flowdata[comkey][inoutkey][key] = append(flowdata[comkey][inoutkey][key], fdt[comkey][inoutkey][key])

        flowdata.update({'r':r,
                        'z1':z1,
                        'a1':a1,
                        't1':t1,
                        'z2':z2,
                        'a2':a2,
                        't2':t2,
                        'dt':dt})

    return flowdata

def getflow(s1, s2, r, Rvir, t, sim, dt, dmo=False):
    #print '   --- Run getflow'
    #returns a dictionary with details of the inflow and outflow through a sphere of radius r
    
    #print '    Get haloes'
    #search the inflowing/outflowing particles in a big sphere instead of the entire simulation, to reduce runtime
    sf1 = s1[pynbody.filt.Sphere(3*Rvir)]
    sf2 = s2[pynbody.filt.Sphere(3*Rvir)]
    
    #draw_halo(sf2, Rvir, qty='rho',filename='ev_rho_' + sim + '_{:.02f}'.format(float(t))+'_95.png')
    #draw_halo(sf2, Rvir, qty='rho_zoom',filename='ev_rho_zoom_' + sim + '_{:.02f}'.format(float(t))+'_95.png')
    
    #print '    Filtering'
    sp = s1[pynbody.filt.Sphere(r)]
    M = pynbody.array.SimArray(sum(sp['mass']))
    M.units = 'Msol'
    vcirc = sqrt(G*M/r).in_units('km s**-1')
    #dr = min((dt*vcirc).in_units('kpc'), r)
    #tdyn = r/vcirc

    #identify flow using a Bridge
    #print '    Define bridge'
    b = pynbody.bridge.OrderBridge(sf1, sf2)
    outer1 = sf1[pynbody.filt.HighPass('r', r)]
    inner1 = sf1[pynbody.filt.LowPass('r', r)]
    gotin = b(outer1)[pynbody.filt.LowPass('r', r)]
    gotout = b(inner1)[pynbody.filt.HighPass('r', r)]    
    
    
    res = {'all':{'in':{},'out':{}},
          'd':{'in':{},'out':{}},
          'g':{'in':{},'out':{}},
          's':{'in':{},'out':{}}}
    
    #print '    Create data array'
    for (pcls, inout) in zip([gotin, gotout], ['in','out']):
        if dmo: comps = zip([pcls],['all'])
        else: comps = zip([pcls, pcls.d, pcls.g, pcls.s],['all','d','g','s'])
        for (pclsc, c) in comps:
            vr = ((pclsc['r']-b(pclsc)['r'])/dt).in_units('km s**-1')
            ad = abs(vr) <= vcirc
            imp = abs(vr) > vcirc
            res[c][inout]['vr'] = mean(vr) #flow velocity
            res[c][inout]['vr_ad'] = mean(vr[ad]) #adiabatic part's velocity
            res[c][inout]['vr_imp'] = mean(vr[imp]) #impulsive part's velocity
            res[c][inout]['vcirc'] = vcirc #circular velocity
            res[c][inout]['m'] = sum(pclsc['mass']) #flow mass
            res[c][inout]['m_ad'] = sum(pclsc[ad]['mass']) #adiabatic part's mass
            res[c][inout]['m_imp'] = sum(pclsc[imp]['mass']) #impulsive part's mass

    return res

###############################################################################
