###############################################################################
#                             GENERAL FUNCTIONS                               #
###############################################################################

from numpy import *
import pynbody

def flush():  
    return sys.stdout.flush()

def array_nonan(array1,array2=array([])):
    if size(array2)==0:
        array2=array1
    if size(array1)<>size(array2):
        raise ValueError('The two arrays must have the same size')
    else:
        n=size(array2)-size(where(isnan(array2)))
        new_array=zeros(n)
        j=0
        for i in range(size(array2)):
            if not isnan(array2[i]):
                new_array[j]=array1[i]
                j=j+1
        return new_array

def interpolate_nan(array):
    output=array.copy()
    nans=isnan(output)
    x = lambda z: z.nonzero()[0]
    output[nans]=interp(x(nans), x(~nans), output[~nans])
    return output

def first_nonan(array):
    return next((index,x) for (index,x) in zip(range(size(array)),array) if ~isnan(x))

def noinf(array):
    array[array == inf]=nan
    array[array == -inf]=nan
    return array

def load_s(sim,ss,centering_com='s'):
    print '   --- load_s %s%s'%(sim,ss)
    s = pynbody.load('/vol/sci/astro/cosmo/nas1/nihao/' + sim + '/' + sim + ss)

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

def load_s_DM(sim,ss):
    print '   --- load_s %s.dm%s'%(sim,ss)
    s = pynbody.load('/vol/sci/astro/cosmo/nas2/Data/nihao/dmo2/' + sim + '/' + sim + '.dm'+ ss)
    s.physical_units()
    
    try: 
        h = s.halos()
        h1 = h[1]
    
        # CENTERING
        pynbody.analysis.halo.center(h1.d,mode='ssc',shrink_factor=0.95)
        pynbody.analysis.angmom.faceon(h1)
    except: 
        print '   --- Warning: no centering for output %s.dm.%s'%(sim,ss)
        
    return s  

def get_fangzhou_radii(sim,a_array,get_all=False,get_stars=False,D200=False):
        '''
        Define different radii from Fangzhou Jiang's files. 
        
        Syntax: 
            If get_all=False (default):
            ok,r12,rvir,mvir=get_fangzhou_radii(sim,a)
            
            If get_all=True
            ok,r12,rvir,mvir,r12_coldgas,r12_coldbar,r12_SFR,rs_NFW,r2_Einasto,rc_Dekel,re_Sersic_star=
                get_fangzhou_radii(sim,a,get_all=True)
        '''
        
        # DEFINE R12, RVIR
        if get_all:
            data_r12=genfromtxt('/cs/sci/freundlich/CUSPCORE/Fangzhou/spinhistory_%s_extended.txt'%sim,usecols=(0,1,2,33,34,35,55,68,76,90,95))
        elif get_stars:
            data_r12=genfromtxt('/cs/sci/freundlich/CUSPCORE/Fangzhou/spinhistory_%s_extended.txt'%sim,usecols=(0,1,2,33,5))
        else:
            data_r12=genfromtxt('/cs/sci/freundlich/CUSPCORE/Fangzhou/spinhistory_%s_extended.txt'%sim,usecols=(0,1,2,33))
        
        ok=[]
        r12=[]
        rvir=[]
        mvir=[]
        
        if get_stars:
            mstar=[]
            
        if get_all:
            r12_coldgas=[]
            r12_coldbar=[]
            r12_SFR=[]
            rs_NFW=[]
            r2_Einasto=[]
            rc_Dekel=[]
            re_Sersic_star=[]
        
        a_fangzhou=data_r12[:,0]
        
        for a in a_array:
            i=where(a_fangzhou==round(a,6))
        
            if size(i)==1:
                    ok.append(True)
                    r12.append(data_r12[:,3][i][0])
                    rvir.append(data_r12[:,2][i][0])
                    mvir.append(data_r12[:,1][i][0])
                    #
                    if get_stars:
                        mstar.append(data_r12[:,4][i][0])
                    #
                    if get_all:
                        r12_coldgas.append(data_r12[:,4][i][0])
                        r12_coldbar.append(data_r12[:,5][i][0])
                        r12_SFR.append(data_r12[:,6][i][0])
                        rs_NFW.append(data_r12[:,7][i][0])
                        r2_Einasto.append(data_r12[:,8][i][0])
                        rc_Dekel.append(data_r12[:,9][i][0])   
                        re_Sersic_star.append(data_r12[:,10][i][0])   
            else:
                    ok.append(False)
                    r12.append(nan)
                    rvir.append(nan)
                    mvir.append(nan)
                    #
                    if get_stars:
                        mstar.append(nan)
                    #
                    if get_all: 
                        r12_coldgas.append(nan)
                        r12_coldbar.append(nan)
                        r12_SFR.append(nan)
                        rs_NFW.append(nan)
                        r2_Einasto.append(nan)
                        rc_Dekel.append(nan) 
                        re_Sersic_star.append(nan)   

        if D200:
            rvir=[]
            mvir=[]
            for a in a_array:
                try:
                    catalog='/cs/sci/freundlich/CUSPCORE/catalogs/NIHAO_a%.4f.txt'%a
                    data=genfromtxt(catalog,skip_header =1)
                    i_ID=where(data[:,0]==float(sim[1:]))[0][0]
                    mvir.append(data[i_ID,7])
                    rvir.append(data[i_ID,8])
                except:
                    mvir.append(nan)
                    rvir.append(nan)                     
        
        # SAVE QUANTITIES INTO ARRAYS
        r12=array(r12)
        rvir=array(rvir)
        mvir=array(mvir)
        
        if get_stars:
            mstar=array(mstar)
        
        if get_all:
            r12_coldgas=array(r12_coldgas)
            r12_coldbar=array(r12_coldbar)
            r12_SFR=array(r12_SFR)
            rs_NFW=array(rs_NFW)
            r2_Einasto=array(r2_Einasto)
            rc_Dekel=array(rc_Dekel)
            re_Sersic_star=array(re_Sersic_star)
        
        if get_all:
            return ok,r12,rvir,mvir,r12_coldgas,r12_coldbar,r12_SFR,rs_NFW,r2_Einasto,rc_Dekel,re_Sersic_star
        elif get_stars:
            return ok,r12,rvir,mvir,mstar
        else: 
            return ok,r12,rvir,mvir

def get_fangzhou_mstar(sim,a_array):
        '''
        Get Mstar from Fangzhou Jiang's files. 
        
        Syntax: 
            mstar=get_fangzhou_radii(sim,a)
        '''
        data=genfromtxt('/cs/sci/freundlich/CUSPCORE/Fangzhou/spinhistory_%s_extended.txt'%sim,usecols=(0,5))
        a_fangzhou=data[:,0]
        
        mstar=[]
        for a in a_array:
            i=where(a_fangzhou==round(a,6))
            if size(i)==1:
                mstar.append(data[:,1][i][0])
            else:
                mstar.append(nan)
        mstar=array(mstar)
        return mstar

def redress_denominator(denominator):
    n=size(denominator)
    n2=n/2
    # Prevent denominator=alpha+gamma-2beta to be zero or negative
    if size(where(denominator<=0)[0])>n2:
        denominator=zeros_like(denominator)
    # Lower half
    elif size(where(denominator[:n2]<=0)[0])>0:
        imin=where(denominator[:n2]<=0)[0][-1]+1
        denominator[:imin]=denominator[imin]*ones(size(denominator[:imin]))
    # Upper half
    elif size(where(denominator[n2:]<=0)[0])>0:
        imin=where(denominator[n2:]<=0)[0][0]-1
        denominator[imin:]=denominator[imin]*ones(size(denominator[imin:]))
            
    return denominator

def linearize(beta,x,rlim):
    r_range=where((x>=rlim[0])&(x<rlim[1]))
    beta_linear=nan*ones(size(beta))
    try:
        p=polyfit(x[r_range],beta[r_range],1)
        beta_linear[r_range]=p[0]*x[r_range]+p[1]
    except:
        print 'No linear fit: nan instead'
    return beta_linear

def linearized_function(x_new,x_ref,beta_ref,rlim):
    r_range=where((x_ref>=rlim[0])&(x_ref<rlim[1]))
    beta_new=nan*ones(size(x_new))
    try:
        p=polyfit(x_ref[r_range],beta_ref[r_range],1)
        beta_new=p[0]*x_new+p[1]
    except:
        print 'No linear fit: nan instead'
    return beta_new
        
    