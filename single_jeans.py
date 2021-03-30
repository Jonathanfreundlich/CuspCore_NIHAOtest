# COMPARE (alpha+gamma-2beta)sigma_r^2 and Vc^2

from matplotlib.pylab import *
G = 4.499753324353496e-06 # gravitational constant [kpc^3 Gyr^-2 Msun^-1]
parsec=3.085677581e16 #m
year=3.1557600e7 #s
kms_to_kpcGyr=1/parsec*year*1e9

fontsize=20
legendsize=18
linewidth=2
linesize=5
component='d'

rcParams['axes.linewidth'] = linewidth
rcParams['xtick.major.size'] = linesize
rcParams['xtick.major.width'] = linewidth
rcParams['xtick.minor.size'] = linesize
rcParams['xtick.minor.width'] = linewidth
rcParams['ytick.major.size'] = linesize
rcParams['ytick.major.width'] = linewidth
rcParams['ytick.minor.size'] = linesize
rcParams['ytick.minor.width'] = linewidth
rcParams['xtick.labelsize'] = fontsize-4
rcParams['ytick.labelsize'] = fontsize-4

figure()
axhline(0,color='gray')

for (i,col,kstring,label) in zip((k,k+1),('blue','red'),('k','k+1'),('before','after')):
    ss=gl[i]
    r= ss[component]['r']
    Rvir=ss['Rvir']
    Mvir=ss[component]['Mvir']
    Kvir=0.5*G*Mvir/Rvir
    logr=log10(ss[component]['r']/ss['Rvir'])
    t=ss['t']    
    logr=log10(r/Rvir)

    # Term1 = (alpha+gamma-2beta)sigma_r^2
    sigma2=(ss[component]['sigmar_smooth']/3.085677581*3.1556952)**2
    
    alpha=ss[component]['alpha']
    beta=ss[component]['beta_smooth']
    gamma=ss[component]['gamma']
    alpha_lin=linearize(alpha,logr,rlim)
    beta_lin=linearize(beta,logr,rlim)
    gamma_lin=linearize(gamma,logr,rlim)
    
    Term1 = (alpha+gamma-2.*beta)*sigma2/Kvir
    Term1_lin = (alpha_lin+gamma_lin-2.*beta_lin)*sigma2/Kvir
    
    # Term2 = Vc2
    vc2=G*ss[component]['Mall']/r
    Term2=vc2/Kvir
    
    line1,=plot(logr,Term1,color=col,linestyle='-',label=r'$\rm K_{1}=(\alpha+\gamma-2\beta)\sigma_r^2$')
    line1lin,=plot(logr,Term1_lin,color=col,linestyle=':')
    
    line2,=plot(logr,Term2,color=col,linestyle='--',label=r'$\rm K_{2} = V_c^2$')

xlabel(r'$\log(r/R_{\rm vir})$',fontsize=fontsize)
ylabel(r'$\rm K/K_{vir}$',fontsize=fontsize) #$\rm [kpc^2 Gyr^{-2}]$
legend((line1,line2),(r'$\rm K_{1}=(\alpha+\gamma-2\beta)\sigma_r^2$',r'$\rm K_{2} = V_c^2$'), fontsize=legendsize,frameon=False,loc='upper left')

ylim(-2,6)
xlim(rlim)

vlines=concatenate((linspace(0.01,0.1,10),linspace(0.2,1,9)))
for xv in vlines:
    axvline(x=log10(xv),color='k',linestyle='-',alpha=0.2)
    
ax=gca()
ax.text(0.03,0.11, r'$\rm blue:$ $k=%i$'%k, transform=ax.transAxes, color='blue',fontsize=fontsize)
ax.text(0.03,0.04, r'$\rm red:$ $k+1=%i$'%(k+1), transform=ax.transAxes, color='red',fontsize=fontsize)
    
############################################

figure()

axhline(0,color='gray')

for (i,col,kstring,label) in zip((k,k+1),('blue','red'),('k','k+1'),('before','after')):
    ss=gl[i]
    r= ss[component]['r']
    Rvir=ss['Rvir']
    Mvir=ss[component]['Mvir']
    Kvir=0.5*G*Mvir/Rvir
    logr=log10(ss[component]['r']/ss['Rvir'])
    t=ss['t']    
    logr=log10(r/Rvir)

    # Term1 = (alpha+gamma-2beta)sigma_r^2
    sigma2=(ss[component]['sigmar_smooth']/3.085677581*3.1556952)**2
    
    alpha=ss[component]['alpha']
    beta=ss[component]['beta_smooth']
    gamma=ss[component]['gamma']
    alpha_lin=linearize(alpha,logr,rlim)
    beta_lin=linearize(beta,logr,rlim)
    gamma_lin=linearize(gamma,logr,rlim)
    
    Term1 = (alpha+gamma-2.*beta)*sigma2/Kvir
    Term1_lin = (alpha_lin+gamma_lin-2.*beta_lin)*sigma2/Kvir
    
    # Term2 = Vc2
    vc2=G*ss[component]['Mall']/r
    Term2=vc2/Kvir
    
    DT=(Term1-Term2)/Term2
    DT_lin=(Term1_lin-Term2)/Term2
    
    line1,=plot(logr,DT,color=col,linestyle='-')
    line1_lin,=plot(logr,DT_lin,color=col,linestyle=':')
    
axhline(0.5,color='gray')
axhline(-0.5,color='gray')
axhline(1,color='gray')
axhline(-1,color='gray')

vlines=concatenate((linspace(0.01,0.1,10),linspace(0.2,1,9)))
for xv in vlines:
    axvline(x=log10(xv),color='k',linestyle='-',alpha=0.2)
    
xlabel(r'$\log(r/R_{\rm vir})$',fontsize=fontsize)
ylabel(r'$\rm (K_{1}-K_{2})/K_{2}$',fontsize=fontsize) #$\rm [kpc^2 Gyr^{-2}]$
legend(fontsize=legendsize,frameon=False,loc='upper left')

ylim(-2.,2.)
xlim(rlim)
ax=gca()
ax.text(0.05,0.9,r'Relative difference',fontsize=fontsize,transform=ax.transAxes)
ax.text(0.03,0.11, r'$\rm blue:$ $k=%i$'%k, transform=ax.transAxes, color='blue',fontsize=fontsize)
ax.text(0.03,0.04, r'$\rm red:$ $k+1=%i$'%(k+1), transform=ax.transAxes, color='red',fontsize=fontsize)

############################################

# COMPARE Kreal and Kjeans = 0.5 (3-2beta)/(alpha+gamma-2beta) Vc^2

figure()
axhline(0,color='gray')

for (i,col,kstring,label) in zip((k,k+1),('blue','red'),('k','k+1'),('before','after')):
    ss=gl[i]
    r= ss[component]['r']
    Rvir=ss['Rvir']
    Mvir=ss[component]['Mvir']
    Kvir=0.5*G*Mvir/Rvir
    logr=log10(ss[component]['r']/ss['Rvir'])
    t=ss['t']    
    logr=log10(r/Rvir)
    p = ss['d'][fitname]['p']

    # Tjeans = 0.5 (3-2beta)/(alpha+gamma-2beta) Vc^2
    vc2=G*ss[component]['Mall']/r
       
    alpha=ss[component]['alpha']
    beta=ss[component]['beta_smooth']
    gamma=ss[component]['gamma']
    alpha_lin=linearize(alpha,logr,rlim)
    beta_lin=linearize(beta,logr,rlim)
    gamma_lin=linearize(gamma,logr,rlim)
    num=redress_denominator(3-2.*beta)
    num_lin=redress_denominator(3-2.*beta_lin)
    den=redress_denominator(alpha+gamma-2.*beta)
    den_lin=redress_denominator(alpha_lin+gamma_lin-2.*beta_lin)
    
    Tjeans=0.5*num/den*vc2/Kvir
    Tjeans_lin=0.5*num_lin/den_lin*vc2/Kvir
    
    # T real
    T_real=treal[i][2]/Kvir
    
    # T_evol
    M_fit = prf.M(r, p)
    M_real=ss[component]['Mall']
    
    add_params=[]
    if Ttype=='jeans-Mreal' or Ttype=='alpha-Mreal' or Ttype=='alpha-p-Mreal':
        M = M_real
        add_params=M_real
    else: 
        M=M_fit
    
    T_evol=get_T(r,[alpha,beta,gamma,p],m=0.,Ttype=Ttype,do_smooth=False,add_params=add_params)/Kvir
    
    line1,=plot(logr,Tjeans,color=col,linestyle='-',label=r'$\rm K_{Jeans}=\frac{3-2\beta}{\alpha+\gamma-2\beta}\frac{GM(r)}{2r}$')
    line1_evol,=plot(logr,T_evol,color=col,linestyle='--',label=r'$\rm K_{model}$')
    line2,=plot(logr,T_real,color=col,linestyle=':',label=r'$\rm K_{real}$')

xlabel(r'$\log(r/R_{\rm vir})$',fontsize=fontsize)
ylabel(r'$\rm K/K_{vir}$',fontsize=fontsize) #$\rm [kpc^2 Gyr^{-2}]$
legend((line1,line1_evol,line2),(r'$\rm K_{Jeans}=\frac{3-2\beta}{\alpha+\gamma-2\beta}\frac{GM(r)}{2r}$',r'$\rm K_{model}$',r'$\rm K_{real}$'),fontsize=legendsize,frameon=False,loc='upper left')

ylim(-2,6)
xlim(rlim)

vlines=concatenate((linspace(0.01,0.1,10),linspace(0.2,1,9)))
for xv in vlines:
    axvline(x=log10(xv),color='k',linestyle='-',alpha=0.2)
    
ax=gca()
ax.text(0.03,0.11, r'$\rm blue:$ $k=%i$'%k, transform=ax.transAxes, color='blue',fontsize=fontsize)
ax.text(0.03,0.04, r'$\rm red:$ $k+1=%i$'%(k+1), transform=ax.transAxes, color='red',fontsize=fontsize)
    
    
    

figure()
axhline(0,color='gray')

for (i,col,kstring,label) in zip((k,k+1),('blue','red'),('k','k+1'),('before','after')):
    ss=gl[i]
    r= ss[component]['r']
    Rvir=ss['Rvir']
    Mvir=ss[component]['Mvir']
    Kvir=0.5*G*Mvir/Rvir
    logr=log10(ss[component]['r']/ss['Rvir'])
    t=ss['t']    
    logr=log10(r/Rvir)

    # Tjeans = 0.5 (3-2beta)/(alpha+gamma-2beta) Vc^2
    vc2=G*ss[component]['Mall']/r
       
    alpha=ss[component]['alpha']
    beta=ss[component]['beta_smooth']
    gamma=ss[component]['gamma']
    alpha_lin=linearize(alpha,logr,rlim)
    beta_lin=linearize(beta,logr,rlim)
    gamma_lin=linearize(gamma,logr,rlim)
    num=redress_denominator(3-2.*beta)
    num_lin=redress_denominator(3-2.*beta_lin)
    den=redress_denominator(alpha+gamma-2.*beta)
    den_lin=redress_denominator(alpha_lin+gamma_lin-2.*beta_lin)
    
    Tjeans=0.5*num/den*vc2/Kvir
    Tjeans_lin=0.5*num_lin/den_lin*vc2/Kvir
    
    # T real
    T_real=treal[i][2]/Kvir

    # T_evol
    M_fit = prf.M(r, p)
    M_real=ss[component]['Mall']
    
    add_params=[]
    if Ttype=='jeans-Mreal' or Ttype=='alpha-Mreal' or Ttype=='alpha-p-Mreal':
        M = M_real
        add_params=M_real
    else: 
        M=M_fit
    
    T_evol=get_T(r,[alpha,beta,gamma,p],m=0.,Ttype=Ttype,do_smooth=False,add_params=add_params)/Kvir
    
    # Relative difference
    DT=(Tjeans-T_real)/T_real
    DT_lin=(Tjeans_lin-T_real)/T_real
    DT_evol=(T_evol-T_real)/T_real
    
    line1,=plot(logr,DT,color=col,linestyle='-')
    line1_lin,=plot(logr,DT_lin,color=col,linestyle=':')
    line2,=plot(logr,DT_evol,color=col,linestyle='--')

xlabel(r'$\log(r/R_{\rm vir})$',fontsize=fontsize)
ylabel(r'$\rm (K_{Jeans}-K_{real})/K_{real}$',fontsize=fontsize) #$\rm [kpc^2 Gyr^{-2}]$

ylim(-2,2)
xlim(rlim)
axhline(0.5,color='gray')
axhline(-0.5,color='gray')
axhline(1,color='gray')
axhline(-1,color='gray')

vlines=concatenate((linspace(0.01,0.1,10),linspace(0.2,1,9)))
for xv in vlines:
    axvline(x=log10(xv),color='k',linestyle='-',alpha=0.2)
    
ax=gca()

ax.text(0.55,0.9,r'Relative difference',fontsize=fontsize,transform=ax.transAxes)
ax.text(0.03,0.11, r'$\rm blue:$ $k=%i$'%k, transform=ax.transAxes, color='blue',fontsize=fontsize)
ax.text(0.03,0.04, r'$\rm red:$ $k+1=%i$'%(k+1), transform=ax.transAxes, color='red',fontsize=fontsize)
       

    





