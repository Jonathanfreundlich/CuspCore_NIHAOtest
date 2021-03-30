import matplotlib.patheffects as PathEffects
import pynbody.plot.sph as sph
from scipy import stats

snapshots = ['.' + str(i).zfill(5) for i in range(16, 1024+16, 16)]

def view(sim,k,components=['d'],angles=['faceon'],styles=['cmap'],figsize=(12,6),clim=[],xcircle=0.1):
    s1=load_s(sim,snapshots[k],'s')
    s2=load_s(sim,snapshots[k+1],'s')
    h1=s1.halos()[1]
    h2=s2.halos()[1]
    t1 = pynbody.array.SimArray(float(s1.properties['time'].in_units('Gyr')))
    t2 = pynbody.array.SimArray(float(s2.properties['time'].in_units('Gyr')))
    Rvir1 = rvir_fangzhou[k]
    Rvir2 = rvir_fangzhou[k+1]
    
    for angle in angles:
        if angle=='faceon':
            pynbody.analysis.angmom.faceon(h1)
            pynbody.analysis.angmom.faceon(h2)     
        elif angle=='sideon':
            pynbody.analysis.angmom.sideon(h1)
            pynbody.analysis.angmom.sideon(h2)     
        else:
            print 'The angle is not right'
        
        for style in styles:
            if style=='cmap':
                for component in components: 
                    figure(figsize=figsize)

                    subplot(121)
                    ax1=gca()
                    txt=ax1.text(0.03,0.8, r'$k = %i$'%(k), transform=ax1.transAxes, color='k',fontsize=14)
                    txt.set_path_effects([PathEffects.withStroke(linewidth=5, foreground='w')])
                    if component=='d':
                        if size(clim)==0:
                            im=sph.image(h1.d,qty="rho",units="g cm^-3",width=80,subplot=ax1,show_cbar=False,ret_im=True)
                            [vmin,vmax]=im.get_clim()
                        else:
                            sph.image(h1.d,qty="rho",units="g cm^-3", width=80, subplot=ax1, show_cbar=False, vmin=clim[0], vmax=clim[1])
                        txt=ax1.text(0.03,0.9, r'Dark matter', transform=ax1.transAxes, color='k',fontsize=14)
                    elif component=='s':
                        if size(clim)==0:
                            im=sph.image(h1.s,qty="rho",units="g cm^-3",width=80,subplot=ax1,show_cbar=False,ret_im=True)
                            [vmin,vmax]=im.get_clim()
                        else:
                            sph.image(h1.s,qty="rho",units="g cm^-3", width=80, subplot=ax1, show_cbar=False, vmin=clim[0], vmax=clim[1])
                        txt=ax1.text(0.03,0.9, r'Stars', transform=ax1.transAxes, color='k',fontsize=14)
                    elif component=='g':
                        if size(clim)==0:
                            im=sph.image(h1.g,qty="rho",units="g cm^-3",width=80,subplot=ax1,show_cbar=False,ret_im=True)
                            [vmin,vmax]=im.get_clim()                          
                        else:
                            im=sph.image(h1.g,qty="rho",units="g cm^-3",width=80,subplot=ax1,show_cbar=False,vmin=clim[0], vmax=clim[1])
                        txt=ax1.text(0.03,0.9, r'Gas', transform=ax1.transAxes, color='k',fontsize=14)
                    else:
                        print 'Components not valid'
                    txt.set_path_effects([PathEffects.withStroke(linewidth=5, foreground='w')])
                    
                    u_st = h1['pos'].units.latex()
                    xlabel("$x/%s$" % u_st)
                    ylabel("$y/%s$" % u_st)
                            
                    ax=gca()
                    circle=Circle((0,0),xcircle*Rvir1,color='k',fill=False)
                    ax.add_artist(circle)
                            
                    if angle=='faceon':
                        txt=ax1.text(0.75,0.9, r'Face on', transform=ax1.transAxes, color='k',fontsize=14)
                    elif angle=='sideon':
                        txt=ax1.text(0.75,0.9, r'Side on', transform=ax1.transAxes, color='k',fontsize=14)
                    txt.set_path_effects([PathEffects.withStroke(linewidth=5, foreground='w')])
                    
                    subplot(122)
                    ax2=gca()
                    txt=ax2.text(0.03,0.8, r'$k +1 = %i$'%(k+1), transform=ax2.transAxes, color='k',fontsize=14)
                    txt.set_path_effects([PathEffects.withStroke(linewidth=5, foreground='w')])
                    if component=='d':
                        sph.image(h2.d,qty="rho",units="g cm^-3",width=80,subplot=ax2,show_cbar=False,vmin=vmin,vmax=vmax)#,cmap="Greys")
                        txt=ax2.text(0.03,0.9, r'Dark matter', transform=ax2.transAxes, color='k',fontsize=14)
                    elif component=='s':
                        sph.image(h2.s,qty="rho",units="g cm^-3",width=80,subplot=ax2,show_cbar=False,vmin=vmin,vmax=vmax)#,cmap="Greys")
                        txt=ax2.text(0.03,0.9, r'Stars', transform=ax2.transAxes, color='k',fontsize=14)
                    elif component=='g':
                        sph.image(h2.g,qty="rho",units="g cm^-3",width=80,subplot=ax2,show_cbar=False,vmin=vmin,vmax=vmax)#,cmap="Greys")
                        txt=ax2.text(0.03,0.9, r'Gas', transform=ax2.transAxes, color='k',fontsize=14)
                    else:
                        print 'Components not valid'
                    txt.set_path_effects([PathEffects.withStroke(linewidth=5, foreground='w')])
                    
                    if angle=='faceon':
                        txt=ax2.text(0.75,0.9, r'Face on', transform=ax2.transAxes, color='k',fontsize=14)
                    elif angle=='sideon':
                        txt=ax2.text(0.75,0.9, r'Side on', transform=ax2.transAxes, color='k',fontsize=14)
                    txt.set_path_effects([PathEffects.withStroke(linewidth=5, foreground='w')])
       
                    ax=gca()
                    circle=Circle((0,0),xcircle*Rvir1,color='k',fill=False)
                    ax.add_artist(circle)
                    
                    show()
                    
            elif style=='HST':
                    figure(figsize=figsize)

                    subplot(121)
                    ax1=gca()
                    ax1.text(0.03,0.9, r'HST image', transform=ax1.transAxes, color='w',fontsize=14)
                    ax1.text(0.03,0.8, r'$k = %i$'%(k), transform=ax1.transAxes, color='w',fontsize=14)
                    pynbody.plot.stars.render(s1,width='20 kpc',axes=ax1,plot=True,clear=False)
                    
                    u_st = s1['pos'].units.latex()
                    xlabel("$x/%s$" % u_st)
                    ylabel("$y/%s$" % u_st)       

                    ax=gca()
                    circle=Circle((0,0),xcircle*Rvir1,color='w',fill=False)
                    ax.add_artist(circle)
                    
                    if angle=='faceon':
                        ax1.text(0.75,0.9, r'Face on', transform=ax1.transAxes, color='w',fontsize=14)
                    elif angle=='sideon':
                        ax1.text(0.75,0.9, r'Side on', transform=ax1.transAxes, color='w',fontsize=14)
                        
                    subplot(122)
                    ax2=gca()
                    ax2.text(0.03,0.9, r'HST image', transform=ax2.transAxes, color='w',fontsize=14)
                    ax2.text(0.03,0.8, r'$k +1 = %i$'%(k+1), transform=ax2.transAxes, color='w',fontsize=14)
                    pynbody.plot.stars.render(s2,width='20 kpc',axes=ax2,plot=True,clear=False)
                    u_st = s1['pos'].units.latex()
                    xlabel("$x/%s$" % u_st)
                    ylabel("$y/%s$" % u_st)

                    ax=gca()
                    circle=Circle((0,0),xcircle*Rvir1,color='w',fill=False)
                    ax.add_artist(circle)
                    
                    if angle=='faceon':
                        ax2.text(0.75,0.9, r'Face on', transform=ax2.transAxes, color='w',fontsize=14)
                    elif angle=='sideon':
                        ax2.text(0.75,0.9, r'Side on', transform=ax2.transAxes, color='w',fontsize=14)
                        
                    show()
                  
            else:
                print 'Style not compatible'

def view_profile(sim,k,components=['d'],figsize=(6,6),textfont=16,plot_prediction=False,rmin_evolve=0.01,rmax_evolve=1,Ttype='Tmulti',use_RMvir=False):
    for comp in components:
        figure(figsize=figsize)
        j=1
        for (i,col,label) in zip((k,k+1),('gray','black'),('before','after')):
            ss = gl[i]
            t = ss['t']
            r = ss[comp]['r']
            brho = ss[comp]['brho']
            brhofit = ss[comp][fitname]['brho']
            p = ss[comp][fitname]['p']

            Rvir=rvir_fangzhou[i]
            Mvir=mvir_fangzhou[i]

            if j==1: 
                brhov = brho[-1]
                pi=p
                ssi=ss
                gammai=ss[component]['gamma']
                betai =ss[component]['beta_smooth']
                alphai =ss[component]['alpha']

                plot(log10(r/Rvir),log10(brho/brhov),color=col,label=r'$k=%i$'%k)
                plot(log10(r/Rvir),log10(brhofit/brhov),'--',color=col,lw=2)  
            
            if j == 2:
                pf=p
                rf=r
                t1 = ss['flowdata']['t1']
                t2 = ss['flowdata']['t2']
                
                plot(log10(r/Rvir),log10(brho/brhov),color=col,label=r'$k=%i$'%(k+1))
                plot(log10(r/Rvir),log10(brhofit/brhov),'--',color=col,lw=2)  
            
                if plot_prediction:
                        # Define m
                        m=[]
                        for i in range(size(rf)):
                            mval=get_m(ss,rf[i])
                            m.append(mval)
                        m=array(m)
                        Mi = prf.M(rf, pi)
                        
                        # EVOLVE PROFILE
                        wcore = 0*(r<=rmin_evolve*Rvir)+1*((r>rmin_evolve*Rvir)&(r<=rmax_evolve*Rvir))+0*(r>rmax_evolve*Rvir)
                        nancore=nan*ones_like(wcore)
                        nancore[where(wcore==1)]=1
                        
                        if Ttype=='Tmulti':
                            M_rmin=max(ss[comp]['r'][0],ss['all']['r'][0])
                            M_d = ss[comp]['M'][where(ss[comp]['r']>M_rmin)]
                            M_a = ss['all']['M'][where(ss[comp]['r']>M_rmin)]
                            M_rr=r[where(ss[comp]['r']>M_rmin)]
                            M_ratio=M_a/M_d

                            slope, intercept,_,_,_ = stats.linregress(log10(M_rr/Rviri),log10(M_ratio))
                            Mratio=10**intercept
                            Mn=-slope
                            add_params=[Mratio,Mn]
                        
                        if use_RMvir:
                                res = evolve(rf, rf, Mi, pi, m, alphai,alphai,gammai,gammai,betai,betai, Ttype=Ttype, w=wcore, method='halo',add_params=add_params,Rvirf=Rvir,Mvirf=Mvir)
                        else:
                                res = evolve(rf, rf, Mi, pi, m, alphai,alphai,gammai,gammai,betai,betai, Ttype=Ttype, w=wcore, method='halo',add_params=add_params)
                                                              
                        p_model = res['pf']
                        brhofit_model = prf.brho(rf,res['pf'])
                        plot(log10(rf/Rvir),log10(brhofit_model/brhov),'-',color='red',label=r'$\rm model$',lw=2)
                            

            j+=1
        
        ax = gca()
        ax.axis(limits)
        ax.text(0.95,0.92, r'${:.01f}$'.format(t1) + r'$\rm \rightarrow$' + r'${:.01f}$'.format(t2) + r' $\rmGyr$', transform=ax.transAxes, color=col,fontsize=textfont,ha='right')
        
        maxr = rmax_fit*Rvirf
        rr=rf[rf<=maxr]
        Darea=fit.Delta_area(log10(rr/Rvirf),log10(prf.brho(rr, pf)/brhov),log10(prf.brho(rr, p_model)/brhov),xlimits=delta_xlim,ymin=delta_ymin)
        Darea_i=fit.Delta_area(log10(rr/Rvir),log10(prf.brho(rr, pi)/brhov),log10(prf.brho(rr, p_model)/brhov),xlimits=delta_xlim,ymin=delta_ymin)
        Darea_if=fit.Delta_area(log10(rr/Rvir),log10(prf.brho(rr, pi)/brhov),log10(prf.brho(rr, pf)/brhov),xlimits=delta_xlim,ymin=delta_ymin)

        maxr = rmax_evolve*Rvirf
        fmeani= mean(abs(array_nonan(fprofile[rf<=maxr])))
        fstdi = std(abs(array_nonan(fprofile[rf<=maxr])))
        
        print 'delta    = ', Darea
        print 'delta_i  = ', Darea_i
        print 'delta_if = ', Darea_if
        print 'fmean    = ', fmeani
        print 'fstd     = ', fstdi
 
   
        logr=log10(r/Rvir)
        ifill=where((logr>delta_xlim[0])&(logr<delta_xlim[1]))
        fill_between(logr[ifill],delta_ymin*ones_like(logr[ifill]),log10(prf.brho(rr, pi)/brhov)[ifill],facecolor='lightcoral',alpha=0.2,lw=0)
        fill_between(logr[ifill],delta_ymin*ones_like(logr[ifill]),log10(prf.brho(rr, pf)/brhov)[ifill],facecolor='lightcoral',alpha=0.2,lw=0)
        if plot_prediction:
            fill_between(logr[ifill],log10(prf.brho(rr, pf)/brhov)[ifill],log10(brhofit_model/brhov)[ifill],facecolor='blue',alpha=0.4,lw=0) #hatch="X",edgecolor="k"
        
        legend(loc='lower left',frameon=False)
        
        xlabel(r'$\log(r/R_{\rmvir})$',fontsize=textfont+2)
        ylabel(r'$\log(\bar{\rho}/\bar{\rho}_{\rmvir})$',fontsize=textfont+2)

def view_fprofile(sim,k,components=['d'],figsize=(6,6),textfont=16,limits=[-2, 0.0, 0, 0.5]):
    for comp in components:
        figure(figsize=figsize)
        j=1
        
        plot(log10(rf/Rvir),nancore*fprofile,color='b',lw=2)
        plot(log10(rf/Rvir),fprofile,color='b',linestyle=':',lw=2)
        plot(log10(rf/Rvir),-nancore*fprofile,color='r',lw=2)
        plot(log10(rf/Rvir),-fprofile,color='r',linestyle=':',lw=2)
        plot(log10(rf/Rvir),ones_like(rf)*fmean[k+1],color='k',lw=1)
        
        ax = gca()
        ax.axis(limits)
        ax.text(0.95,0.92, r'${:.01f}$'.format(t1) + r'$\rm \rightarrow$' + r'${:.01f}$'.format(t2) + r' $\rmGyr$', transform=ax.transAxes, color=col,fontsize=textfont,ha='right')
        ax.text(-0.7,fmean[k+1]+0.01, r'$|f|_{\rm RMS} = %.2f$'%fmean[k+1], fontsize = textfont,color='k')
        
        xlabel(r'$\log(r/R_{\rmvir})$',fontsize=textfont+2)
        ylabel(r'$|f|$',fontsize=textfont+2)

def view_kprofile(sim,k,components=['d'],figsize=(6,6),textfont=16,limits=[-2, 0, -0.5, 0.5]):
    for comp in components:
        figure(figsize=figsize)
        j=1
        
        Tviri=0.5*prf.G*Mviri/Rviri
        Tvirf=0.5*prf.G*Mvirf/Rvirf
        
        # Treal
        for (i,linestyle) in zip((k,k+1),('-','--')):
            r_real=treal[i][0]
            Rvir_real=treal[i][1]
            T_real=treal[i][2]
            Tline,=plot(log10(r_real/Rvir_real),log10(T_real/Tviri),'r',linestyle=linestyle,lw=2)
            if i==k:
                Tline.set_label(r'$K$')

        plot(log10(rf/Rviri),log10(nancore*Ti_evol/Tviri),'b',linestyle='-',label=r'$K_{\rm multi}$',lw=2)
        plot(log10(rf/Rviri),log10(Ti_evol/Tviri),'b',linestyle=':',lw=2)
        
        plot(log10(rf_evol/Rvirf),log10(nancore*Tf_evol/Tviri),'b',linestyle='--',lw=2)
        plot(log10(rf_evol/Rvirf),log10(Tf_evol/Tviri),'b',linestyle=':',lw=2)
       

        ax=gca()
        ax.axis(limits)
        ax.text(0.95,0.92, r'${:.01f}$'.format(t1) + r'$\rm \rightarrow$' + r'${:.01f}$'.format(t2) + r' $\rmGyr$', transform=ax.transAxes, color=col,fontsize=textfont,ha='right')
                
        legend(loc='lower left',frameon=False)

        xlabel(r'$\log(r/R_{\rmvir})$',fontsize=textfont+2)
        ylabel(r'$\log(K/K_{\rm vir})$',fontsize=textfont+2)
        
def view_fkprofile(sim,k,components=['d'],figsize=(6,12),textfont=16,limits1=[-2, 0.0, 0, 0.4],limits2=[-2, 0, -0.5, 0.5]):
    for comp in components:
        figure(figsize=figsize)
        
        ax1=subplot(211)
        
        plot(log10(rf/Rvir),nancore*fprofile,color='b',lw=2)
        plot(log10(rf/Rvir),fprofile,color='b',linestyle=':',lw=2)
        plot(log10(rf/Rvir),-nancore*fprofile,color='r',lw=2)
        plot(log10(rf/Rvir),-fprofile,color='r',linestyle=':',lw=2)
        plot(log10(rf/Rvir),ones_like(rf)*fmean[k+1],color='k',lw=1)
        
        ax1.axis(limits1)
        ax1.text(0.95,0.92, r'${:.01f}$'.format(t1) + r'$\rm \rightarrow$' + r'${:.01f}$'.format(t2) + r' $\rmGyr$', transform=ax1.transAxes, color=col,fontsize=textfont,ha='right')
        ax1.text(-0.7,fmean[k+1]+0.01, r'$|f|_{\rm RMS} = %.2f$'%fmean[k+1], fontsize = textfont,color='k')
        
        setp(ax1.get_xticklabels(), visible=False)
        ylabel(r'$|f|$',fontsize=textfont+2)
        
        ax2=subplot(212, sharex=ax1)
        Tviri=0.5*prf.G*Mviri/Rviri
        Tvirf=0.5*prf.G*Mvirf/Rvirf
        
        # Treal
        for (i,linestyle) in zip((k,k+1),('-','--')):
            r_real=treal[i][0]
            Rvir_real=treal[i][1]
            T_real=treal[i][2]
            Tline,=plot(log10(r_real/Rvir_real),log10(T_real/Tviri),'r',linestyle=linestyle,lw=2)
            if i==k:
                Tline.set_label(r'$K$')

        plot(log10(rf/Rviri),log10(nancore*Ti_evol/Tviri),'b',linestyle='-',label=r'$K_{\rm multi}$',lw=2)
        plot(log10(rf/Rviri),log10(Ti_evol/Tviri),'b',linestyle=':',lw=2)
        
        plot(log10(rf_evol/Rvirf),log10(nancore*Tf_evol/Tviri),'b',linestyle='--',lw=2)
        plot(log10(rf_evol/Rvirf),log10(Tf_evol/Tviri),'b',linestyle=':',lw=2)
       
        ax2.axis(limits2)        
        legend(loc='lower left',frameon=False)

        xlabel(r'$\log(r/R_{\rmvir})$',fontsize=textfont+2)
        ylabel(r'$\log(K/K_{\rm vir})$',fontsize=textfont+2)
        
        subplots_adjust(wspace=0, hspace=0.1)
