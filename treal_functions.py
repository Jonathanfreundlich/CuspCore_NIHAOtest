###############################################################################
#                  DETERMINING TREAL AND SAVING THE ARRAY                     #
###############################################################################

'''
Creates Treal such that for each snapshot i, Treal[i]=[r,Rvir,Tr]

Notes:

- The radii r are taken from the files NIHAO-xxx.pickle containing the profile and flow data

- It is possible to choose Rvir from the simulation haloes or from Fangzhou's data
 
'''

import os
import sys
from numpy import *
import pickle
import pynbody

import format_functions
from format_functions import *
import general_functions
from general_functions import *

from scipy.signal import savgol_filter
from scipy.integrate import simps

# gravitational constant [kpc^3 Gyr^-2 Msun^-1]
G = 4.499753324353496e-06 

from matplotlib.pylab import *
rcParams['figure.figsize'] = (10,6)
rcParams['font.size'] = 18

snapshots = ['.' + str(i).zfill(5) for i in range(16, 1024+16, 16)]

###############################################################################
# LOAD GL FILE IF IT EXISTS, CREATE IT OTHERWISE

def load_or_create_gl(sim,centering_com='s',directory='/cs/sci/freundlich/CUSPCORE/NIHAO_data/',name='Treal',c='d',remove=False,do_plot=False,use_fangzhou_Rvir=True,D200=False): 

    if c<>'d':
        name=name+'-'+c
    
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
            i_start=shape(gl)[0]         
        else:
            gl=[]
            i_start=0
        
        # Load remaining snapshots
        r_all=get_all_r(sim,directory,use_fangzhou_Rvir=use_fangzhou_Rvir,D200=D200)
        for (ss, i) in zip(snapshots, range(len(snapshots))):
            if i<i_start:
                print sim, 'snap=',ss[1:],'output=',i, '-- already defined'
            else:
                print sim, 'snap=',ss[1:], 'output=',i
                
                r=r_all[i]
                r,Rvir,Tr=get_Treal(sim,i,r,c,use_fangzhou_Rvir=use_fangzhou_Rvir,D200=D200)[:3]
                data=[r,Rvir,Tr]
                    
                gl.append(data)
                save_gl_temp(sim,gl,directory,name)
                flush()
                
                if do_plot:
                    plot_treal(data,rlim=[-2,0],title_label='%s.%s'%(sim,ss[1:]))

        save_gl(sim,gl,directory,name=name)
    return gl

def load_gl(sim,directory,name='Treal'):
    with open(directory+name+'-%s.pickle'%(sim[1:])) as f:
        gl = pickle.load(f)
    return gl

def save_gl(sim,gl,directory,remove_file=True,name='Treal'):
    if remove_file:
        try:
            os.remove(directory+name+'-%s.pickle'%(sim[1:]))
        except:
            'No file to remove'
    with open(directory+name+'-%s.pickle'%(sim[1:]), 'w') as f:
        pickle.dump(gl, f)

def load_gl_temp(sim,directory,name='Treal'):
    with open(directory+name+'-%s-temp.pickle'%(sim[1:])) as f:
        gl = pickle.load(f)
    return gl   

def save_gl_temp(sim,gl,directory,name='Treal'):
    with open(directory+name+'-%s-temp.pickle'%(sim[1:]), 'w') as f:
        pickle.dump(gl, f)

def flush():  
    return sys.stdout.flush()

def remove_files(sim,directory,name='Treal'):
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

# GET ALL RADII AT WHICH TO EVALUATE Treal

def get_all_r(sim,directory,use_fangzhou_Rvir=True,D200=False):
    gl=format_functions.load_or_create_gl(sim,directory=directory,name='NIHAO',use_fangzhou_Rvir=use_fangzhou_Rvir,D200=D200)
    r_all=[]
    for k in range(len(snapshots)):
        ss=gl[k]
        r = array(ss['all']['r'])
        r_all.append(r)
    return r_all

###############################################################################

# EVALUATE Treal

def get_Treal(sim,k,r,c='d',polyorder=3,sigma = 11,centering_com='s',enlarge_radii=False,use_fangzhou_Rvir=True,D200=False):
  try:
    ss=snapshots[k]
    s=load_s(sim,ss,centering_com)

    s.physical_units()
    a = s.properties['a']
    
    # get general parameters
    if use_fangzhou_Rvir:
        Rvir=get_fangzhou_radii(sim,array([a]),get_all=False,D200=D200)[2][0]
    else:
        h = s.halos()
        h1 = h[1]
        Rvir = float(max(h1['r']))
    
    if c=='d':
        sc=s.d
    elif c=='s':
        sc=s.s
    elif c=='g':
        sc=s.g
    elif c=='all':
        sc=s
    else:
        print 'Error: wrong component c'

    # ENLARGE RADII
    if enlarge_radii:
        dlogr=diff(log10(r))[0]
        num=int(log10(2.)/dlogr)
        rb=concatenate((r,10**linspace(max(log10(r))+dlogr,max(log10(r))+num*dlogr,num=num)))
    else:
        rb=r
    
    # REAL KINETIC ENERGY

    K=zeros(size(rb))
    M=zeros(size(rb))
    for i in range(size(rb)):
        si=sc[pynbody.filt.Sphere(rb[i])]
        K[i]=0.5*sum(si['mass']*si['v2'])
        M[i]=sum(si['mass'])

    try:
        r2=min(where(log10(rb/Rvir)>=-2.5)[0])
        istart=max(first_nonan(noinf(K))[0],r2)
        Tr=derive_Tr(rb,K,M,istart,sigma,polyorder)
    except:
        print 'Warning: Tr not defined'
        Tr=nan*ones(size(rb))

    # RETRIEVE INITIAL RADII IF enlarge_radii
    Tr=Tr[:size(r)]

  except:
    print '   The output could not be loaded'
    r=array([nan])
    Rvir=nan
    Tr=array([nan])
    rb=array([nan])
    K=array([nan])
    M=array([nan])

  return r,Rvir,Tr,rb,K,M

def derive_Tr(rb,K,M,istart,sigma,polyorder):
    K_smooth=concatenate((nan*ones(istart),savgol_filter(interpolate_nan(K[istart:]),sigma,polyorder,deriv=0)))
    dK=concatenate((nan*ones(istart), savgol_filter(K_smooth[istart:],sigma,polyorder,deriv=1,delta=diff(log10(rb[istart:]))[0])))
    dK_smooth=concatenate((nan*ones(istart),savgol_filter(interpolate_nan(dK[istart:]),sigma,polyorder,deriv=0)))

    M_smooth=concatenate((nan*ones(istart),savgol_filter(interpolate_nan(M[istart:]),sigma,polyorder,deriv=0)))
    dM=concatenate((nan*ones(istart), savgol_filter(M_smooth[istart:],sigma,polyorder,deriv=1,delta=diff(log10(rb[istart:]))[0])))
    dM_smooth=concatenate((nan*ones(istart),savgol_filter(interpolate_nan(dM[istart:]),sigma,polyorder,deriv=0)))
           
    Tr=dK_smooth/dM_smooth
    return Tr

###############################################################################

def Uout_log(ss,r):
    '''
    Second term of the expression of the gravitational potential
    Uout = -int_r^R_vir 4 pi G rho x^2 dlnx 
    r in kpc
    rho in Msol/kpc^3
    G in kpc^3 Gyr^-2 Msun^-1
    > Uout in kpc^2 Gyr^-2
    '''
    rho=ss['all']['rho']
    rr=ss['all']['r']
    Uout=zeros_like(r)
    
    
    if size(r)>1:
        dlnr=log(r[1]/r[0])
        for i in range(size(r)):
            Uout[i]=- sum(np.array(4.*pi*G*rr**2*rho*dlnr)[where(rr>r[i])])
    
    return Uout

def Uout(ss,r):
    '''
    Second term of the expression of the gravitational potential
    Uout = -int_r^R_vir 4 pi G rho x dx 
    - r in kpc
    - rho in Msol/kpc^3
    - G in kpc^3 Gyr^-2 Msun^-1
    - Uout in kpc^2 Gyr^-2
    '''
    rho=ss['all']['rho']
    rr=ss['all']['r']
    Uout=zeros_like(r)
    if size(r)>1:
        for i in range(size(r)):
            indices=where(rr>r[i])
            Uout[i]=-trapz(np.array(4.*pi*G*rr*rho)[indices],x=rr[indices])
            
    return Uout
    
def Uin(ss,r):
    '''
    Uin = GM(<r)/r
    '''
    M=np.array(ss['all']['M'])
    rr=np.array(ss['all']['r'])
    r=np.array(r)
    Uin=zeros_like(r)
    if size(r)>0:
        for i in range(size(r)):
            Uin[i]=-G*M[abs(rr-r[i])==min(abs(rr-r[i]))]/rr[abs(rr-r[i])==min(abs(rr-r[i]))]
    
    return Uin

def Utot(ss,r):
    '''
    Utot=Uin+Uout
    '''
    return Uin(ss,r)+Uout(ss,r)

def Utot_int(ss,r):
    M=np.array(ss['all']['M'])
    rr=np.array(ss['all']['r'])
    r=np.array(r)
    Utot=zeros_like(r)
    if size(r)>1:
        for i in range(size(r)):
            indices=where(rr>r[i])
            Utot[i]=-trapz(G*M[indices]/rr[indices]**2,x=rr[indices])-G*M[-1]/rr[-1]
            #Utot[i]=-simps(G*M[indices]/rr[indices]**2,x=rr[indices])-G*M[-1]/rr[-1]
    return Utot

###############################################################################

# PLOT Treal
def plot_treal(data,rlim=[-2,0],title_label='',subplot=False):
    [r,Rvir,Tr]=data
    x=log10(r/Rvir)
    index=(x>rlim[0])&(x<rlim[1])
    
    if not subplot:
        print 'No subplot'
        figure()
    plot(x[index],Tr[index])
    xlabel(r'$\rm \log(r/R_{vir})$')
    ylabel(r'$\rm T_{real}$')
    title(r'%s'%title_label)
    
    if not subplot:
        show()
