# ENERGIES USED TO EVOLVE THE MODEL

# ALL ENERGIES
figure()

# Real kinetic energy
for (i,linestyle) in zip((k,k+1),('-','--')):
    r_real=treal[i][0]
    Rvir_real=treal[i][1]
    T_real=treal[i][2]
    Tline,=plot(log10(r_real/Rviri),T_real/Tviri,'k',linestyle=linestyle)

Tline.set_label(r'$\rm T_{real}$')

# Model kinetic energy
plot(log10(rf/Rviri),nancore*Ti_evol/Tviri,'r',linestyle='-',label=r'$\rm T_{model}$')
plot(log10(rf/Rviri),Ti_evol/Tviri,'r',linestyle=':')
plot(log10(rf_evol/Rvirf),nancore*Tf_evol/Tviri,'r',linestyle='--')
plot(log10(rf_evol/Rvirf),Tf_evol/Tviri,'r',linestyle=':')

# Model potential energy
plot(log10(rf/Rviri),nancore*Ui_evol/Tviri,'b',linestyle='-',label=r'$\rm U_{model}(p)$')
plot(log10(rf/Rviri),Ui_evol/Tviri,'b',linestyle=':')
plot(log10(rf_evol/Rvirf),nancore*Uf_evol/Tviri,'b',linestyle='--')
plot(log10(rf_evol/Rvirf),Uf_evol/Tviri,'b',linestyle=':')
     
# Potential variation
plot(log10(rf/Rviri), nancore*G*m/rf/Tviri,'g',linestyle='-',label=r'$\rm Gm/r$')
plot(log10(rf/Rviri), G*m/rf/Tviri,'g',linestyle=':')    
plot(log10(rf_evol/Rvirf),nancore*G*m/rf_evol/Tviri,'g',linestyle='--')
plot(log10(rf_evol/Rvirf),G*m/rf_evol/Tviri,'g',linestyle=':')

# Total energy
plot(log10(rf/Rviri),nancore*Ei_evol/Tviri,'m',linestyle='-',label=r'$\rm E_{model}=U_{model}-Gm/r+T_{model}$')
plot(log10(rf/Rviri),Ei_evol/Tviri,'m',linestyle=':')
plot(log10(rf_evol/Rvirf),nancore*Ef_evol/Tviri,'m',linestyle='--')
plot(log10(rf_evol/Rvirf),Ef_evol/Tviri,'m',linestyle=':')

vlines=concatenate((linspace(0.01,0.1,10),linspace(0.2,1,9)))
for xv in vlines:
    axvline(x=log10(xv),color='k',linestyle='-',alpha=0.2)

ax=gca()
ax.set_xlim(limits_T[0], limits_T[1])
ax.set_ylim(-12, 4)

legend(fontsize=textfont,loc='lower right')
xlabel(r'$\log(r/R_{\rm vir})$',fontsize=textfont+2)
ylabel(r'$E/T_{\rm vir}$',fontsize=textfont+2)
title(r'Energies from the profile fits')
