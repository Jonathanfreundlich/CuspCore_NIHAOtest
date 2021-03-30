# PLOT ALPHA, BETA, GAMMA

figure()

'''
ssi = gl[k]
ri = ssi['d']['r']
pi = ssi['d'][fitname]['p']
(ci, ai, bi, gi, Rviri, Mviri)=pi   
alphai = ssi['d']['alpha']
alphai_model = prf.alpha_Dekel(ri,pi)
betai = ssi['d']['beta_smooth']
gammai = ssi['d']['gamma']

ssf = gl[k+1]
rf = ssf['d']['r']
pf = ssf['d'][fitname]['p']
(cf, af, bf, gf, Rvirf, Mvirf)=pf  
alphaf = ssf['d']['alpha']
alphaf_model = prf.alpha_Dekel(rf,pf)
betaf = ssf['d']['beta_smooth']
gammaf = ssf['d']['gamma']

if linear_slopes:
    alphai=linearize(alphai,log10(ri/Rviri),rlim)
    betai=linearize(betai,log10(ri/Rviri),rlim)
    gammai=linearize(gammai,log10(ri/Rviri),rlim)
    alphaf=linearize(alphaf,log10(rf/Rvirf),rlim)
    betaf=linearize(betaf,log10(rf/Rvirf),rlim)
    gammaf=linearize(gammaf,log10(rf/Rvirf),rlim)

t1 = ssf['flowdata']['t1']
t2 = ssf['flowdata']['t2']
'''

alphai_model = prf.alpha_Dekel(ri,pi)

plot(log10(ri/Rviri),nancore*alphai_evol/4.,'r',linestyle='-',label=r'$\alpha/4$')
plot(log10(ri/Rviri),alphai_evol/4.,'r',linestyle=':')

plot(log10(ri/Rviri),nancore*alphai_model/4.,'r',linestyle='--',dashes=(3, 1))
plot(log10(ri/Rviri),alphai_model/4.,'r',linestyle=':')

plot(log10(rf/Rvirf),nancore*alphaf_evol/4.,'r',linestyle='--')
plot(log10(rf/Rvirf),alphaf_evol/4.,'r',linestyle=':')

plot(log10(ri/Rviri),nancore*betai_evol,'b',linestyle='-',label=r'$\beta$')
plot(log10(ri/Rviri),betai_evol,'b',linestyle=':')

plot(log10(rf/Rvirf),nancore*betaf_evol,'b',linestyle='--')
plot(log10(rf/Rvirf),betaf_evol,'b',linestyle=':')

plot(log10(ri/Rviri),nancore*gammai_evol,'g',linestyle='-',label=r'$\gamma$')
plot(log10(ri/Rviri),gammai_evol,'g',linestyle=':')

plot(log10(rf/Rvirf),nancore*gammaf_evol,'g',linestyle='--')
plot(log10(rf/Rvirf),gammaf_evol,'g',linestyle=':')

axhline(0)

ax = gca()
ax.axis(limits_abg)
ax.text(0.03,0.9, r'$\rm %i \rightarrow %i$: '%(k,k+1), transform=ax.transAxes, color='k',fontsize=textfont)
ax.text(0.03,0.83, r'${:.01f}$'.format(t1) + r'$\rm \rightarrow$' + r'${:.01f}$'.format(t2) + r'$\rmGyr$', transform=ax.transAxes, color='k',fontsize=textfont)
           
legend(loc='lower right')
xlabel(r'$\log(r/R_{\rmvir})$',fontsize=textfont+2)




# PLOT 

figure()

numerator_i=3.-2.*betai
denominator_i=get_denominator(ri,[alphai,betai,gammai,pi],Ttype='jeans')
denominator_model_i=get_denominator(ri,[alphai,betai,gammai,pi],Ttype='jeans-alpha')

numerator_f=3.-2.*betaf
denominator_f=get_denominator(rf,[alphaf,betaf,gammaf,pf],Ttype='jeans')
denominator_model_f=get_denominator(rf,[alphaf,betaf,gammaf,pf],Ttype='jeans-alpha')

plot(log10(ri/Rviri),nancore*numerator_i,'r',linestyle='-',label=r'$3-2\beta_i$')
plot(log10(ri/Rviri),numerator_i,'r',linestyle=':')

plot(log10(ri/Rviri),nancore*denominator_i,'b',linestyle='-',label=r'$\alpha_i+\gamma_i-2\beta_i$')
plot(log10(ri/Rviri),denominator_i,'b',linestyle=':')

plot(log10(ri/Rviri),nancore*denominator_model_i,'b',linestyle='--',dashes=(3, 1),label=r'$\alpha_{\rm Dekel,i}+\gamma_i-2\beta_i$')
plot(log10(ri/Rviri),denominator_model_i,'b',linestyle=':')

plot(log10(rf/Rvirf),nancore*numerator_f,'orange',linestyle='-',label=r'$3-2\beta_f$')
plot(log10(rf/Rvirf),numerator_f,'orange',linestyle=':')

plot(log10(rf/Rvirf),nancore*denominator_f,'green',linestyle='-',label=r'$\alpha_f+\gamma_f-2\beta_f$')
plot(log10(rf/Rvirf),denominator_f,'green',linestyle=':')

plot(log10(rf/Rvirf),nancore*denominator_model_f,'green',linestyle='--',dashes=(3, 1),label=r'$\alpha_{\rm Dekel,f}+\gamma_f-2\beta_f$')
plot(log10(rf/Rvirf),denominator_model_f,'green',linestyle=':')


axhline(0,color='k')

ax = gca()
ax.axis([-2,0,-1,6])
ax.set_xlim([limits_abg[0], limits_abg[1]])
ax.text(0.03,0.9, r'$\rm %i \rightarrow %i$: '%(k,k+1), transform=ax.transAxes, color='k',fontsize=textfont)
ax.text(0.03,0.83, r'${:.01f}$'.format(t1) + r'$\rm \rightarrow$' + r'${:.01f}$'.format(t2) + r'$\rmGyr$', transform=ax.transAxes, color='k',fontsize=textfont)

legend(fontsize=14,ncol=2,loc='lower right')
xlabel(r'$\log(r/R_{\rmvir})$',fontsize=textfont+2)
