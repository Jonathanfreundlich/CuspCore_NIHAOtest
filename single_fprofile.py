# M PROFILE
                    
figure()
plot(log10(rf/Rvir),nancore*m/Mvir,color='b')
plot(log10(rf/Rvir),m/Mvir,color='b',linestyle=':')

plot(log10(rf/Rvir),-nancore*m/Mvir,color='r')
plot(log10(rf/Rvir),-m/Mvir,color='r',linestyle=':')

vlines=concatenate((linspace(0.01,0.1,10),linspace(0.2,1,9)))
for xv in vlines:
    axvline(x=log10(xv),color='k',linestyle='-',alpha=0.2)

ax = gca()
ax.axis([-2,0,0,0.03])
ax.text(0.03,0.9, r'$\rm %i \rightarrow %i$: '%(k,k+1), transform=ax.transAxes, color='k',fontsize=textfont)
ax.text(0.03,0.83, r'${:.01f}$'.format(t1) + r'$\rm \rightarrow$' + r'${:.01f}$'.format(t2) + r'$\rmGyr$', transform=ax.transAxes, color='k',fontsize=textfont)
ax.text(0.03,0.20, r'$\Delta snap={:.2f}$'.format(Dsnapshots[k+1]), fontsize = textfont, transform=ax.transAxes,color='k')
ax.text(0.03,0.13, r'$|f|={:.2f} \pm {:.2f} $'.format(fmean[k+1],fstd[k+1]), fontsize = textfont, transform=ax.transAxes,color='k')
ax.text(0.03,0.06, r'$\Delta fit={:.2f}$'.format(delta[k+1]), fontsize = textfont, transform=ax.transAxes,color='k')

xlabel(r'$\log(r/R_{\rmvir})$',fontsize=textfont+2)
ylabel(r'$\rm |m/M_{\rm vir}|$',fontsize=textfont+2)


# F PROFILE
            
figure()
plot(log10(rf/Rvir),nancore*fprofile,color='b')
plot(log10(rf/Rvir),fprofile,color='b',linestyle=':')

plot(log10(rf/Rvir),-nancore*fprofile,color='r')
plot(log10(rf/Rvir),-fprofile,color='r',linestyle=':')

vlines=concatenate((linspace(0.01,0.1,10),linspace(0.2,1,9)))
for xv in vlines:
    axvline(x=log10(xv),color='k',linestyle='-',alpha=0.2)
axhspan(0, merger_thr, alpha=0.2, color='k')

ax = gca()
ax.axis(limits_fprofile)
ax.text(0.03,0.9, r'$\rm %i \rightarrow %i$: '%(k,k+1), transform=ax.transAxes, color='k',fontsize=textfont)
ax.text(0.03,0.83, r'${:.01f}$'.format(t1) + r'$\rm \rightarrow$' + r'${:.01f}$'.format(t2) + r'$\rmGyr$', transform=ax.transAxes, color='k',fontsize=textfont)
ax.text(0.03,0.20, r'$\Delta snap={:.2f}$'.format(Dsnapshots[k+1]), fontsize = textfont, transform=ax.transAxes,color='k')
ax.text(0.03,0.13, r'$|f|={:.2f} \pm {:.2f} $'.format(fmean[k+1],fstd[k+1]), fontsize = textfont, transform=ax.transAxes,color='k')
ax.text(0.03,0.06, r'$\Delta fit={:.2f}$'.format(delta[k+1]), fontsize = textfont, transform=ax.transAxes,color='k')

xlabel(r'$\log(r/R_{\rmvir})$',fontsize=textfont+2)
ylabel(r'$|f|$',fontsize=textfont+2)