###############################################################################
#                              PROFILES                                       #
###############################################################################

'''         NOTE THAT THE DEKEL-ZHAO PROFILE USES THE CONVENTION 1/b        '''

from numpy import *
from scipy.integrate import quad
from math import factorial
from scipy.special import betainc, gamma
from scipy import stats

# gravitational constant [kpc^3 Gyr^-2 Msun^-1]
G = 4.499753324353496e-06

###############################################################################

# DENSITY
def rho(r, params, model='an'):
    if model == 'an': 
        (c, a, b, g, Rvir, Mvir) = params
        rs = float(Rvir)/c
        x = r/rs
        return (q_bar(c, a, b, g, Mvir, Rvir)*float(3-a)/3*(1+float(3-g)/(3-a)*pow(x,b))/
                array(pow(x,a)*pow((1+pow(x,b)),(g-a)/b+1)))

    if model == 'dbl':
        (f, c1, c2, a, b, g, Rvir, Mvir) = params
        return rho(r, (c1, a, b, g, Rvir, f*Mvir), 'an') + rho(r, (c2, a, b, g, Rvir, (1-f)*Mvir), 'an')
    
    if model == 'nfw':
        (c, Rvir, Mvir) = params
        rs = float(Rvir)/c
        x = r/rs
        qb = Mvir*pow(c,3)/(4*pi*pow(Rvir,3)*(log(1+c)-c/(1+c)))
        return qb/(x*pow((1+x),2)) 
    
    if model == 'enfw':
        (c, a, b, g, Rvir, Mvir) = params
        rs = float(Rvir)/c
        x = r/rs
        raw_rho = 1/(pow(x,a)*pow((1+pow(x,b)),(g-a)/b))
        q = Mvir/calc_total_mass(r,raw_rho,Rvir)
        return q*raw_rho
    
    if model == 'ein':
        (c, n, Rvir, Mvir) = params
        rs = float(Rvir)/c
        h = rs/float(pow(2*n,n))    
        raw_rho = exp(-pow((r/h),pow(n,-1)))
        q = Mvir/calc_total_mass(r,raw_rho,Rvir)
        return q*raw_rho

# MEAN DENSITY    
def brho(r, params, model='an'):
    if model == 'an': 
        (c, a, b, g, Rvir, Mvir) = params
        rs = float(Rvir)/c
        x = r/rs
        return q_bar(c, a, b, g, Mvir, Rvir)/array(pow(x,a)*pow((1+pow(x,b)),(g-a)/b))
        
    if model == 'nfw':
        (c, Rvir, Mvir) = params
        rs = float(Rvir)/c
        x = r/rs
        return M(r, params, model)/(4.*pi/3.*pow(r,3))

# SLOPE OF THE DENSITY
def s(r, params, model='an'):
    '''
    Warning: old definition of b, i.e, b=1/b
    '''
    if model == 'an': 
        (c, a, b, g, Rvir, Mvir) = params
        rs = float(Rvir)/c
        x = r/rs
        return -(3-g)/float(3-a)*b*pow(x,b)/(1+(3-g)/float(3-a)*b*pow(x,b)) + (a+(g+b)*pow(x,b))/array(1+pow(x,b))
    
    if model == 'dbl':
        (f, c1, c2, a, b, g, Rvir, Mvir) = params
        rho_tot = rho(r, (f, c1, c2, a, b, g, Rvir, Mvir), 'dbl')
        rho1 = rho(r, (c1, a, b, g, Rvir, f*Mvir), 'an')
        rho2 = rho(r, (c2, a, b, g, Rvir, (1-f)*Mvir), 'an')
        s1 = s(r, (c1, a, b, g, Rvir, f*Mvir), 'an')
        s2 = s(r, (c2, a, b, g, Rvir, (1-f)*Mvir), 'an')
        return (rho1*s1+rho2*s2)/rho_tot

# SLOPE OF THE MEAN DENSITY    
def bs(r, params, model='an'):
    if model == 'an': 
        (c, a, b, g, Rvir, Mvir) = params
        rs = float(Rvir)/c
        x = r/rs
        return (a+g*pow(x,b))/array(1+pow(x,b))
    
    if model == 'dbl':
        (f, c1, c2, a, b, g, Rvir, Mvir) = params
        brho_tot = brho(r, (f, c1, c2, a, b, g, Rvir, Mvir), 'dbl')
        brho1 = brho(r, (c1, a, b, g, Rvir, f*Mvir), 'an')
        brho2 = brho(r, (c2, a, b, g, Rvir, (1-f)*Mvir), 'an')
        bs1 = bs(r, (c1, a, b, g, Rvir, f*Mvir), 'an')
        bs2 = bs(r, (c2, a, b, g, Rvir, (1-f)*Mvir), 'an')
        return (brho1*bs1+brho2*bs2)/brho_tot         

# ENCLOSED MASS
def M(r, params, model='an'):
    if model == 'an': 
        (c, a, b, g, Rvir, Mvir) = params
        rs = float(Rvir)/c
        x = r/rs
        return Mvir*mu(c, a, b, g)*pow(x,3-a)/array(pow(1+pow(x,b),(g-a)/b))
    
    if model == 'dbl':
        (f, c1, c2, a, b, g, Rvir, Mvir) = params
        return M(r, (c1, a, b, g, Rvir, f*Mvir), 'an') + M(r, (c2, a, b, g, Rvir, (1-f)*Mvir), 'an')
    
    if model == 'nfw':
        (c, Rvir, Mvir) = params
        rs = float(Rvir)/c
        x = r/rs
        return float(Mvir)*(log(1+pow(x,1.))-pow(x,1.)/(1+pow(x,1.)))/(log(1+c)-c/(1+c))

# ORBITAL VELOCITY
def V(r, params, model='an'):
    if model == 'an':     
        (c, a, b, g, Rvir, Mvir) = params
        rs = float(Rvir)/c
        x = r/rs
        Vvir2 = G*Mvir/Rvir
        return sqrt(Vvir2*mu(c, a, b, g)*c*pow(x,2-a)/array(pow(1+pow(x,b),(g-a)/b)))
    
    if model == 'dbl':
        (f, c1, c2, a, b, g, Rvir, Mvir) = params
        return sqrt(V(r, (c1, a, b, g, Rvir, f*Mvir), 'an')**2 + V(r, (c2, a, b, g, Rvir, (1-f)*Mvir), 'an')**2)

# GRAVITATIONAL FORCE
def F(r, params, model='an'):
    if model == 'an':
        (c, a, b, g, Rvir, Mvir) = params
        rs = float(Rvir)/c
        x = r/rs
        Fvir = G*Mvir/Rvir**2
        return Fvir*mu(c, a, b, g)*c**2*pow(x,1-a)/array(pow(1+pow(x,b),(g-a)/b))
    
    if model == 'dbl':
        (f, c1, c2, a, b, g, Rvir, Mvir) = params
        return F(r, (c1, a, b, g, Rvir, f*Mvir), 'an') + F(r, (c2, a, b, g, Rvir, (1-f)*Mvir), 'an')

# POTENTIAL
def U(r, p, model='an'): 
    if model == 'an': 
        (c, a, b, g, Rvir, Mvir) = p
        rs = float(Rvir)/c
        x = r/rs
        Vvir2 = G*Mvir/Rvir
        if (b==1 and g==3): #analytic case
            return c*mu(c,a,b,g)/float(2-a)*Vvir2*(pow(x/array(1+x),2-a)-pow(c/float(1+c),2-a))-Vvir2
        if (b==0.5 and g==3): #analytic case
            chi = pow(x,0.5)/array(1+pow(x,0.5))
            chic = pow(c,0.5)/array(1+pow(c,0.5))
            return -Vvir2-2*c*mu(c,a,b,g)*Vvir2*((pow(chic,2*(2-a))-pow(chi,2*(2-a)))/float(2*(2-a))-(pow(chic,2*(2-a)+1)-pow(chi,2*(2-a)+1))/float(2*(2-a)+1))
        else:
            dU = lambda y: G*M(y,p,model)/array(y**2)
            return array([-Vvir2-quad(dU,rt,Rvir)[0] for rt in r])
    
    if model == 'dbl':
        (f, c1, c2, a, b, g, Rvir, Mvir) = params
        return U(r, (c1, a, b, g, Rvir, f*Mvir), 'an') + U(r, (c2, a, b, g, Rvir, (1-f)*Mvir), 'an')

# POTENTIAL ENERGY
def Ep(Rmax, p, model='an'):
    dE = lambda y: 0.5*rho(y,p,model)*U(y,p,model)*4*pi*y**2
    return quad(dE,0,Rmax)[0]

# KINETIC ENERGY FROM THE VIRIAL THEOREM
def T_fid(r,params):
    '''
    0.5*G*M/r
    '''
    (c, a, b, g, Rvir, Mvir)=params
    Mi=M(r,params)
    return 0.5*G*Mi/r

# DENSITY SLOPE
def alpha_Dekel(r,params):
    ''' 
    alpha=-dln(rho)/dln(r) the local slope
    Note: there was a mistake in Guy's function s(r,params)!!!
    '''
    (c, a, b, g, Rvir, Mvir) = params
    Rs = float(Rvir)/c
    x = r/Rs
    term1 = (a+(g+b)*pow(x,b))/array(1+pow(x,b))
    term2 = (3-g)/float(3-a)*b*pow(x,b)/(1+(3-g)/float(3-a)*pow(x,b))
    return term1-term2

# CONCENTRATION PARAMETER
def cmax(params):
    '''
    Warning: old definition of b, i.e, b=1/b
    '''
    (c, a, b, g, Rvir, Mvir) = params
    return c*((g-2.)/(2.-a))**(1/b)

def c2(params):
    '''
    Warning: old definition of b, i.e, b=1/b
    '''
    (c, a, b, g, Rvir, Mvir) = params
    return c*(1.5/(2.-a))**2
###############################################################################

# GET RADII FROM MASS PROFILE
def inv_M(M, params, model='an'):
    if model == 'an': 
        (c, a, b, g, Rvir, Mvir) = params
        if g==3: #analytic case
            x = pow(1/array(pow(mu(c, a, b, g)*Mvir/M, b/float(3-a))-1),1/float(b))
            rs = float(Rvir)/c
            return rs*x

###############################################################################

# AUXILIARY FUNCTIONS

def q_bar(c, a, b, g, Mvir, Rvir):
    brho_vir = Mvir/(4*pi/3*pow(Rvir,3))
    return brho_vir*pow(c,a)*pow((1+pow(c,b)),(g-a)/b)

def calc_total_mass(r, rho, Rmax):
    return sum(rho[r<=Rmax]*4*math.pi*r[r<=Rmax]**2*concatenate(([r[0]], diff(r[r<=Rmax]))))

def mu(c, a, b, g):
    return pow(1+pow(c,b),(g-a)/b)/float(pow(c,3-a))


###############################################################################

# VELOCITY DISPERSION

def sigma_r_sqr_Dekel_iso(x,Mv,Rv,c,alpha):
    r"""
    radial velocity dispersion squared of a Dekel+17 halo at radius r --
    for isotropic velocity distribution, beta=0
        
    where:
        x: r/rs (scalar)
        Mv: halo mass [Msun] (scalar)
        c: concentration (scalar)
        alpha: inner density slope (scalar)
        z: redshift (scalar, default=0.)    
    """
    Vv_sqr = G * Mv / Rv
    X = chi(x)
    u = 4*(1.-alpha)
    return 2.*Vv_sqr *c/g(c,alpha) *(x**3.5)/(X**(2.*(3.5-alpha))) \
        * ( (1.-X**u)/u - 8.*(1.-X**(u+1.))/(u+1.) \
        + 28.*(1.-X**(u+2.))/(u+2.) - 56.*(1.-X**(u+3.))/(u+3.) \
        + 70.*(1.-X**(u+4.))/(u+4.) - 56.*(1.-X**(u+5.))/(u+5.) \
        + 28.*(1.-X**(u+6.))/(u+6.) - 8.*(1.-X**(u+7.))/(u+7.) \
        + (1.-X**(u+8.))/(u+8.))

def chi(x):
    r"""
    Auxiliary function for Dekel+17 profile
    
        chi := x^0.5 / 1+x^0.5  
    
    Syntax:
        chi(x)
    where 
        x: dimensionless radius r/rs (scalar or array)
    """
    u = x**0.5
    return u/(1.+u)

def g(x,alpha):
    r"""
    Auxiliary function for Dekel+17 profile
    
        g(x;alpha):= chi^2(3-alpha), with chi := x^0.5 / 1+x^0.5  
    
    Syntax:
        g(x,alpha)
    where 
        x: dimensionless radius r/rs (scalar or array)
        alpha: inner density slope (scalar)
    """
    return chi(x)**(2.*(3.-alpha))

def B9(a,x):
    """
    Incomplete beta function B(a,b,x) with b=9
    betainc(a,9,x)
    """
    B=0.
    for i in range(9):
        B+=factorial(8)/factorial(i)*gamma(a)/gamma(a+9-i)*x**(a+8-i)*(1.-x)**i
    return B

def Bfangzhou(a,x):
    """
    Formula by Fangzhou
    """
    B=0.
    for i in range(9):
        B+=(-1)**i*factorial(8)/factorial(i)/factorial(8-i)*(1-x**(a+i))/(a+i)
    return B

def Binc(a,b,x):
    return betainc(a,b,x)*gamma(a)*gamma(b)/gamma(a+b)
    
def sigmar2_dekel(x,Mv,Rv,c,alpha):
    """
    radial velocity dispersion squared of a Dekel+17 halo at radius r --
    for isotropic velocity distribution, beta=0
        
    where:
        x: r/rs (scalar)
        Mv: halo mass [Msun] (scalar)
        c: concentration (scalar)
        alpha: inner density slope (scalar)
        z: redshift (scalar, default=0.)    
    """
    Vv_sqr = G * Mv / Rv
    X = chi(x)
    u = 4*(1.-alpha)
    return 2.*Vv_sqr *c/g(c,alpha) * x**alpha*(1.+x**0.5)**(2.*(3.5-alpha)) \
        * (Binc(u,9,1)-Binc(u,9,X))  # or (B9(u,1)-B9(u,X))

def sigmar2_dekel_B9(x,Mv,Rv,c,alpha):
    """
    radial velocity dispersion squared of a Dekel+17 halo at radius r --
    for isotropic velocity distribution, beta=0
        
    where:
        x: r/rs (scalar)
        Mv: halo mass [Msun] (scalar)
        c: concentration (scalar)
        alpha: inner density slope (scalar)
        z: redshift (scalar, default=0.)    
    """
    Vv_sqr = G * Mv / Rv
    X = chi(x)
    u = 4*(1.-alpha)
    return 2.*Vv_sqr *c/g(c,alpha) * x**alpha*(1.+x**0.5)**(2.*(3.5-alpha)) \
        * (B9(u,1)-B9(u,X))    #(Binc(u,9,1)-Binc(u,9,X))  # or 
    
def sigmar2_dekel_m(x,Mvir,Rvir,c,a,m=0,mtype='center'):
    X = chi(x)
    Xc= chi(c)
    factor=2.*G*Mvir/Rvir* x**a*(1.+x**0.5)**(2.*(3.5-a))
    #
    u = 4*(1.-a)
    term1=factor*c/g(c,a)*(Binc(u,9,1)-Binc(u,9,X))
    #
    term2=0.
    if mtype=='center':
        u=-2.-2.*a
        term2=factor*c*m/Mvir*(Binc(u,9,1)-Binc(u,9,X))
    elif mtype=='uniform':
        u=-2.*a
        term2=factor*m/Mvir*(Binc(u,7,Xc)-Binc(u,7,X))
        u=-2.-2.*a
        term2+=factor*c*m/Mvir*(Binc(u,9,1)-Binc(u,9,Xc))
    elif mtype=='density':
        u=4*(1.-a)
        term2=factor*c*m/Mvir*(Binc(u,9,1)-Binc(u,9,X))
    else:
        print 'mtype not correct'
        
    #
    return term1+term2

def K_Mratio(x,Mvir,Rvir,c,a,Mratio=1.,n=0.,m=0.,mtype='center'):
    """
    local kinetic energy for a Dekel+17 halo at radius r --
    for isotropic velocity distribution, beta=0
    
    where we assume Mtot/Mdm=Mratio*(r/Rvir)**-n
    
    and we take into account the additional mass m
    
    where:
        x: r/rs (scalar)
        Mvir: halo mass [Msun] (scalar)
        Rvir: virial radius [kpc] (scalar)
        c: concentration (scalar)
        a: inner density slope (scalar)
        Mratio: Mtot/Mdm at Rvir (scalar, default=1.)  
        n: exponent of Mtot/Mdm=Mratio*(r/Rvir)**-n (scalar, default=0.)
        m: additional mass (vector, default=0.)
        mtype: string
       
    """
    X = chi(x)
    Xc= chi(c)
    u = 4*(1.-a)-2*n
    b=0.5
    g=3.
    factor=3.*G*Mvir/Rvir*c*mu(c, a, b, g)*x**a*(1.+x**0.5)**(2.*(3.5-a))
    if u>0: # slightly more rapid
        K0=factor*Mratio*c**n*(Binc(u,9+2*n,1)-Binc(u,9+2*n,X))
    else:
        K0=factor*Mratio*c**n*Dbetainc(u,9+2*n,X)#(Binc(u,9+2*n,1)-Binc(u,9+2*n,X))
        
    
    Km=0.
    if mtype=='center':
        u=-2.-2.*a
        if u>0: # slightly more rapid
            Km=factor*g(c,a)*m/Mvir*(Binc(u,9,1)-Binc(u,9,X))
        else:
            Km=factor*g(c,a)*m/Mvir*Dbetainc(u,9,X)#(Binc(u,9,1)-Binc(u,9,X))
    #elif mtype=='uniform':
    #    u=-2.*a
    #    term2=factor*m/Mvir*(Binc(u,7,Xc)-Binc(u,7,X))
    #    u=-2.-2.*a
    #    term2+=factor*c*m/Mvir*(Binc(u,9,1)-Binc(u,9,Xc))
    #elif mtype=='density':
    #    u=4*(1.-a)
    #    term2=factor*c*m/Mvir*(Binc(u,9,1)-Binc(u,9,X))
    else:
        print 'mtype not correct'  
    #
    return K0+Km

def K(r,params, model='an',Mratio=1.,n=0.,m=0.,rm=0.,mtype='center'):
    if model == 'an': 
        (c, a, b, g, Rvir, Mvir) = params
        rs = float(Rvir)/c
        x = r/rs
        X = chi(x)
        Xc= chi(c)
        factor=3.*G*Mvir/Rvir*c*mu(c, a, b, g)*x**a*(1.+x**0.5)**(2.*(3.5-a))
       
        u = 4*(1.-a)-2*n
        if u>0: # slightly more rapid
            K0=factor*Mratio*c**n*(Binc(u,9+2*n,1)-Binc(u,9+2*n,X))
        else:
            K0=factor*Mratio*c**n*Dbetainc(u,9+2*n,X)#(Binc(u,9+2*n,1)-Binc(u,9+2*n,X))
        
    
        Km=0.
        if mtype=='center':
            u=-2.-2.*a
            if u>0: # slightly more rapid
                Km=factor*m/(mu(c, a, b, g)*Mvir)*(Binc(u,9,1)-Binc(u,9,X))
            else:
                Km=factor*m/(mu(c, a, b, g)*Mvir)*Dbetainc(u,9,X)#(Binc(u,9,1)-Binc(u,9,X))
        elif mtype=='uniform':
            xm=rm/rs
            Xm=chi(xm)
            Km=nan*ones_like(r)
            Km[r>=rm]=Km_function(r[r>=rm],params, model='an',m=m,mtype='center')
            # First term
            u=4.-2.*a
            if u>0: 
                DB=(Binc(u,3,Xm)-Binc(u,3,X))# slightly more rapid
            else:
                DB=Dbetainc(u,3,X,Xm)
            Km1=factor*m/(mu(c, a, b, g)*Mvir)/xm**3*DB
            # Second term
            u=-2.-2.*a
            if u>0: 
                DB=(Binc(u,9,1)-Binc(u,9,Xm))# slightly more rapid
            else:
                DB=Dbetainc(u,9,Xm)
            Km2=factor*m/(mu(c, a, b, g)*Mvir)*DB                
            Km[r<rm]=Km1[r<rm]+Km2[r<rm]      
        else:
            print 'mtype not correct'  
    #
    else:
        print 'model not correct'
    return K0+Km

def Km_function(r,params, model='an',m=0.,rm=0,mtype='center'):
    if model == 'an':
        (c, a, b, g, Rvir, Mvir) = params
        rs = float(Rvir)/c
        x = r/rs
        X = chi(x)
        Xc= chi(c)
        factor=3.*G*Mvir/Rvir*c*mu(c, a, b, g)*x**a*(1.+x**0.5)**(2.*(3.5-a))

        Km=0.
        if mtype=='center':
            u=-2.-2.*a
            if u>0: # slightly more rapid
                Km=factor*m/(mu(c, a, b, g)*Mvir)*(Binc(u,9,1)-Binc(u,9,X))
            else:
                Km=factor*m/(mu(c, a, b, g)*Mvir)*Dbetainc(u,9,X)#(Binc(u,9,1)-Binc(u,9,X))
        elif mtype=='uniform':
            xm=rm/rs
            Xm=chi(xm)
            Km=nan*ones_like(r)
            Km[r>=rm]=Km_function(r[r>=rm],params, model='an',m=m,mtype='center')
            # First term
            u=4.-2.*a
            if u>0: 
                DB=(Binc(u,3,Xm)-Binc(u,3,X))# slightly more rapid
            else:
                DB=Dbetainc(u,3,X,Xm)
            Km1=factor*m/(mu(c, a, b, g)*Mvir)/xm**3*DB
            # Second term
            u=-2.-2.*a
            if u>0: 
                DB=(Binc(u,9,1)-Binc(u,9,Xm))# slightly more rapid
            else:
                DB=Dbetainc(u,9,Xm)
            Km2=factor*m/(mu(c, a, b, g)*Mvir)*DB                
            Km[r<rm]=Km1[r<rm]+Km2[r<rm]          
        else:
            print 'mtype not correct'  
    
    return Km

def U(r, p, model='an'): 
    if model == 'an': 
        (c, a, b, g, Rvir, Mvir) = p
        rs = float(Rvir)/c
        x = r/rs
        Vvir2 = G*Mvir/Rvir
        if (b==1 and g==3): #analytic case
            return c*mu(c,a,b,g)/float(2-a)*Vvir2*(pow(x/array(1+x),2-a)-pow(c/float(1+c),2-a))-Vvir2
        if (b==0.5 and g==3): #analytic case
            chi = pow(x,0.5)/array(1+pow(x,0.5))
            chic = pow(c,0.5)/array(1+pow(c,0.5))
            return -Vvir2-2*c*mu(c,a,b,g)*Vvir2*((pow(chic,2*(2-a))-pow(chi,2*(2-a)))/float(2*(2-a))-(pow(chic,2*(2-a)+1)-pow(chi,2*(2-a)+1))/float(2*(2-a)+1))
        else:
            dU = lambda y: G*M(y,p,model)/array(y**2)
            return array([-Vvir2-quad(dU,rt,Rvir)[0] for rt in r])
    
def Um(r, mtype='center',m=0,rm=0): 
    if mtype=='center':
        Um=-G*m/r
    elif mtype=='uniform':
        Um=nan*ones_like(r)
        Um[r>=rm]=-G*m/r[r>=rm]
        Um[r<rm]=-0.5*G*m/rm*(3.-(r[r<rm]/rm)**2)
    return Um
        
def Bintegrand(t,a,b):
    return t**(a-1)*(1-t)**(b-1)

def Dbetainc(a,b,x1,x2=1):
    if size(x1)==1:
        return quad(Bintegrand,x1,x2,args=(a,b,))[0]
    elif size(x1)>1:
        Dbeta=nan*zeros_like(x1)
        for i in range(size(x1)):
            Dbeta[i]=quad(Bintegrand,x1[i],x2,args=(a,b,))[0]
        return Dbeta
    