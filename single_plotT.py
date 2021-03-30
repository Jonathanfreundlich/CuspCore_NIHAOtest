# PLOT T

figure()

Tviri=0.5*prf.G*Mviri/Rviri
Tvirf=0.5*prf.G*Mvirf/Rvirf
Tfidi=prf.T_fid(rf,pi)
Tfidf=prf.T_fid(rf_evol,pf)
DU=Uf_evol-Ui_evol

plot(log10(rf/Rviri),nancore*DU/Tviri,'b',linestyle='-',label=r'$\Delta U$')
plot(log10(rf/Rviri),DU/Tviri,'b',linestyle=':')
plot(log10(rf/Rviri),-nancore*DU/Tviri,'b',linestyle='--')
plot(log10(rf/Rviri),-DU/Tviri,'b',linestyle=':')

# Treal
for (i,linestyle) in zip((k,k+1),('-','--')):
    r_real=treal[i][0]
    Rvir_real=treal[i][1]
    T_real=treal[i][2]
    Tline,=plot(log10(r_real/Rvir_real),T_real/Tviri,'k',linestyle=linestyle)

Tline.set_label(r'$\rm T_{real}$')

plot(log10(rf/Rviri),nancore*Ti_evol/Tviri,'r',linestyle='-',label=r'$\rm T_{model}$')
plot(log10(rf/Rviri),Ti_evol/Tviri,'r',linestyle=':')
plot(log10(rf_evol/Rvirf),nancore*Tf_evol/Tviri,'r',linestyle='--')
plot(log10(rf_evol/Rvirf),Tf_evol/Tviri,'r',linestyle=':')
plot(log10(rf/Rviri),nancore*Tfidi/Tviri,'g',linestyle='-',label=r'$T_{fid}$')
plot(log10(rf/Rviri),Tfidi/Tviri,'g',linestyle=':')
plot(log10(rf_evol/Rvirf),nancore*Tfidf/Tviri,'g',linestyle='--')
plot(log10(rf_evol/Rvirf),Tfidf/Tviri,'g',linestyle=':')

vlines=concatenate((linspace(0.01,0.1,10),linspace(0.2,1,9)))
for xv in vlines:
    axvline(x=log10(xv),color='k',linestyle='-',alpha=0.2)
    
ax=gca()
ax.axis(limits_T)
ax.text(0.03,0.9, r'$\rm %i \rightarrow %i$: '%(k,k+1), transform=ax.transAxes, color='k',fontsize=textfont)
ax.text(0.3,0.9, r'${:.01f}$'.format(t1) + r'$\rm \rightarrow$' + r'${:.01f}$'.format(t2) + r'$\rmGyr$', transform=ax.transAxes, color='k',fontsize=textfont)

legend()

xlabel(r'$\log(r/R_{\rmvir})$',fontsize=textfont+2)
