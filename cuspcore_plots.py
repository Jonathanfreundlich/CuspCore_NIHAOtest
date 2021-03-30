###############################################################################
#                  PLOTS WHEN RUNNING THE CUSP-CORE MODEL                     #
###############################################################################

from numpy import *
from matplotlib.pylab import *
from decimal import Decimal
from scipy import stats

sys.path.insert(0, '/cs/sci/freundlich/CUSPCORE/Analysis/General')
import general_functions
from general_functions import *
reload(general_functions)

sys.path.insert(0, '/cs/sci/freundlich/CUSPCORE/Analysis/Model')
import evolving
import cuspcore_aux
reload(evolving)
reload(cuspcore_aux)
from evolving import *
from cuspcore_aux import *

###############################################################################

def icol(ncols,nrows):
    position=zeros((nrows,ncols))
    for j in range(ncols):
        position[:,j]=j
    return position.flatten()

def irow(ncols,nrows):
    position=zeros((nrows,ncols))
    for i in range(nrows):
        position[i,:]=i
    return position.flatten()

###############################################################################

# COLOR SUBPLOT AND COUNT SUCCESSES
def do_color_subplot(ax,counts,plot_criterion,criteria,criteria_threshold,bgcolor_big_failure='moccasin',bgcolor_big_success ='lightgreen',bgcolor_small_failure='lightgoldenrodyellow',bgcolor_small_success='white'):
    t1,Dsnapshot,fmerger,delta,delta_i,delta_if,fmean = criteria
    t_min,Dsnapshot_threshold,merger_thr,Dfit_threshold,fmean_min = criteria_threshold
    [nearly,nochange,nmergers,nsuccess]=counts
    
    # COLOR TO REMOVE BELOW SNAPSHOT THRESHOLD
    if plot_criterion=='Dsnapshot':
        # Flag first snapshots
        if (t1<t_min):
            bgcolor='lightgreen'
            ax.set_axis_bgcolor('lightgreen')
            nearly=nearly+1
        # Flag snapshots where nothing happens
        elif (Dsnapshot<Dsnapshot_threshold):
            bgcolor='lightblue'
            ax.set_axis_bgcolor('lightblue')
            nochange=nochange+1
        # Flag mergers
        elif (abs(fmerger) > merger_thr):
            col = 'red'
            bgcolor='pink'
            ax.set_axis_bgcolor('pink')
            nmergers=nmergers+1
        # Flag when the model is wrong
        elif (delta>Dfit_threshold):
            col = 'black'
            bgcolor='lightcoral'
            ax.set_axis_bgcolor('lightcoral')            
        # Flag when the model is closer to the initial profile 
        # (while the final and initial profiles can be distinguished)
        elif (delta_if>Dsnapshot_threshold) and (delta>delta_i):
            col = 'black'
            bgcolor='lightgoldenrodyellow'
            ax.set_axis_bgcolor('lightgoldenrodyellow') 
        else:
            col = 'black'
            bgcolor='white'
            nsuccess=nsuccess+1
    # COLOR TO REMOVE BELOW F THRESHOLD
    elif plot_criterion=='f':
        # Flag first snapshotss
        if (t1<t_min):
            col='black'
            bgcolor='lightgreen'
            ax.set_axis_bgcolor('lightgreen')
            nearly=nearly+1
        # Flag snapshots where f below f_min
        #elif (abs(f)<f_min):
        #    bgcolor='lightblue'
        #    ax.set_axis_bgcolor('lightblue')
        #    nochange=nochange+1
        # Flag mergers
        elif (abs(fmerger) > merger_thr):
            col = 'red'
            bgcolor='pink'
            ax.set_axis_bgcolor('pink')  
            nmergers=nmergers+1  
        # Flag when the model is wrong
        elif (fmean>fmean_min) and (delta>Dfit_threshold):
            col = 'black'
            bgcolor=bgcolor_big_failure #'moccasin'
            ax.set_axis_bgcolor(bgcolor)  
        elif (fmean>fmean_min) and (delta>delta_i) and (delta_if>0.03):
            col = 'black'
            bgcolor=bgcolor_big_failure #'moccasin'
            ax.set_axis_bgcolor(bgcolor) 
        elif (fmean<fmean_min) and (delta>Dfit_threshold):
            col = 'black'
            bgcolor=bgcolor_small_failure #'lightgoldenrodyellow'
            ax.set_axis_bgcolor(bgcolor) 
        # Flag when the model is closer to the initial profile 
        # (while the final and initial profiles can be distinguished)
        elif (fmean<fmean_min) and (delta>delta_i) and (delta_if>0.03):
            col = 'black'
            bgcolor=bgcolor_small_failure #'lightgoldenrodyellow' 
            ax.set_axis_bgcolor(bgcolor)
        else: # SUCCESSES
            col = 'black' 
            bgcolor=bgcolor_small_success #'white'
            if (fmean>fmean_min):
                bgcolor=bgcolor_big_success #'lightgreen'
                ax.set_axis_bgcolor(bgcolor) 
            nsuccess=nsuccess+1 
    else:
        print 'Warning: no plot criterion was specified'
           
    countsf=[nearly,nochange,nmergers,nsuccess]
    return countsf,col,bgcolor
          
###############################################################################

# PLOT F PROFILE
    
def do_plot_fprofile(relevent,xprofile_relevent,fprofile_relevent,plot_criterion,limits_fprofile,fitname,Dsnapshots,fmerger,delta,delta_i,delta_if,fmean,fstd,criteria_threshold,counts,savefigure,params,t1_array,t2_array,nrows=8,ncols=8,figsize=(20,20),textfont=12,title=[],savedir='.',do_color=True,bgcolor_big_failure='moccasin',bgcolor_big_success ='lightgreen',bgcolor_small_failure='lightgoldenrodyellow',bgcolor_small_success='white'):
    [nearly,nochange,nmergers,nsuccess]=counts
    sim,xi,xi_merger,rmax_fit,rmax_evolve,merger_thr,t_min,Dsnapshot_threshold, limits,constrain_fit,constrain_evolution,nsuccess_all,nsample_all = params
    t_min,Dsnapshot_threshold,merger_thr,Dfit_threshold,fmean_min = criteria_threshold
    
    fig = figure(figsize=figsize)
    for index in range(size(relevent[0])):
        
        subplot(nrows,ncols,index+1)
        k=relevent[0][index]-1
        t1=t1_array[k+1]
        t2=t2_array[k+1]
               
        plot(log10(xprofile_relevent[index]),fprofile_relevent[index],color='b',lw=2)
        plot(log10(xprofile_relevent[index]),-fprofile_relevent[index],color='r',lw=2)
        axhline(y=fmean[k+1],color='k',linestyle='-')

        ax = gca()
        criteria = t1,    Dsnapshots[k+1],     fmerger[k+1], delta[k+1],delta_i[k+1],delta_if[k+1],fmean[k+1]
        criteria_threshold = t_min, Dsnapshot_threshold, merger_thr,   Dfit_threshold,fmean_min
        counts             = [nearly,nochange,nmergers,nsuccess]
        if do_color:
            do_color_subplot(ax,counts,plot_criterion,criteria,criteria_threshold,bgcolor_big_failure=bgcolor_big_failure,bgcolor_big_success =bgcolor_big_success,bgcolor_small_failure=bgcolor_small_failure,bgcolor_small_success=bgcolor_small_success)                                        

        
        ############################################################
        # INDICATE SOME INFORMATION ON THE PLOTS

        ax.text(0.95,0.9, r'$\rm %i \rightarrow %i$ / '%(k,k+1) +r'${:.01f}$'.format(t1) + r'$\rm \rightarrow$' + r'${:.01f}$'.format(t2) + r'$\rm Gyr$', transform=ax.transAxes, color='k',fontsize=textfont,ha='right')
        ax.text(-0.78,fmean[k+1]+0.015, r'$|f|_{\rm RMS} = %.2f$'%fmean[k+1], fontsize = textfont,color='k')
        
        if irow(ncols,nrows)[index]==nrows-1:
            if icol(ncols,nrows)[index]==0:
                    ax.set_xticks([-2,-1.5,-1,-0.5])
                    ax.set_xticklabels([r'$-2$',r'$-1.5$',r'$-1$',r'$-0.5$'], fontdict={'fontsize':textfont})
            elif icol(ncols,nrows)[index]==ncols-1:
                    ax.set_xticks([-2,-1.5,-1,-0.5,0])
                    ax.set_xticklabels([r'$0/-2$',r'$-1.5$',r'$-1$',r'$-0.5$',r'$0$'], fontdict={'fontsize':textfont})
            else:
                    ax.set_xticks([-2,-1.5,-1,-0.5])
                    ax.set_xticklabels([r'$0/-2$',r'$-1.5$',r'$-1$',r'$-0.5$'], fontdict={'fontsize':textfont})
        
        if icol(ncols,nrows)[index]==0:
            if irow(ncols,nrows)[index]==0:
                    ax.set_yticks([-0.1,0.,0.1,0.2,0.3,0.4,0.5])
                    ax.set_yticklabels([r'',r'$0/0.5$',r'$0.1$',r'$0.2$',r'$0.3$',r'$0.4$',r'$0.5$'], fontdict={'fontsize':textfont})
            elif irow(ncols,nrows)[index]==nrows-1:
                    ax.set_yticks([-0.1,0.,0.1,0.2,0.3,0.4])
                    ax.set_yticklabels([r'',r'$0$',r'$0.1$',r'$0.2$',r'$0.3$',r'$0.4$'], fontdict={'fontsize':textfont})
            else:
                    ax.set_yticks([-0.1,0.,0.1,0.2,0.3,0.4])
                    ax.set_yticklabels([r'',r'$0/0.5$',r'$0.1$',r'$0.2$',r'$0.3$',r'$0.4$'], fontdict={'fontsize':textfont})
          
       
        ax.axis(limits_fprofile)
        lastrow = size(relevent[0])-ncols
        if not(index >= lastrow): ax.set_xticklabels([])
        if index >= lastrow: xlabel(r'$\log(r/R_{\rmvir})$',fontsize=textfont+2)
        if not(Decimal(str(index)) % Decimal(str(ncols)) == 0): ax.set_yticklabels([])
        if Decimal(str(index)) % Decimal(str(ncols)) == 0: ylabel(r'$|f|$',fontsize=textfont+2)
        setp(ax.get_yticklabels()[0], visible=False)

    fig.subplots_adjust(hspace=0)
    fig.subplots_adjust(wspace=0)
    if size(title)>0:
        suptitle(title)
    else:
        do_plot_title(sim,relevent,fmean,delta,counts,criteria_threshold)

    show()
    if savefigure:
        fig.savefig(savedir+'/%s_fprofile.pdf'%sim.replace('.',''),bbox_inches='tight')
        
###############################################################################

def do_plot_T(sim, relevent,gl,Treal,plot_criterion,limits_T,fitname,Dsnapshots,fmerger,delta,delta_i,delta_if,fmean,criteria_threshold,counts,savefigure,nrows=8,ncols=8,figsize=(20,20),textfont=12,title=[],Ttype='jeans-alpha',savedir='.',do_color=True,component='d',bgcolor_big_failure='moccasin',bgcolor_big_success ='lightgreen',bgcolor_small_failure='lightgoldenrodyellow',bgcolor_small_success='white'):
    # Requires relevent, gl,plot_criterion,limits_T,Dsnapshots,fmerger,delta,f,criteria_threshold,savefigure
    [nearly,nochange,nmergers,nsuccess]=counts
    
    fig = figure(figsize=figsize)
    for index in range(size(relevent[0])):    
                    subplot(nrows,ncols,index+1)
                    k=relevent[0][index]-1
                    
                    ssi = gl[k]
                    ssf = gl[k+1]
                    ri = ssi[component]['r']
                    rf = ssf[component]['r']
                    pi = ssi[component][fitname]['p']
                    pf = ssf[component][fitname]['p']
                    (ci, ai, bi, gi, Rviri, Mviri)=pi
                    (cf, af, bf, gf, Rvirf, Mvirf)=pf
                    t1 = ssf['flowdata']['t1']
                    t2 = ssf['flowdata']['t2']
               
                    Tviri=0.5*prf.G*Mviri/Rviri
                    Tvirf=0.5*prf.G*Mvirf/Rvirf
                    
                    Mi_fit = prf.M(rf, pi)
                    Mi_real=[]
                    for i in range(size(rf)):
                        R=rf[i]
                        Mi_real.append(ssi[component]['M'][abs(ssi[component]['r'] - R) == min(abs(ssi[component]['r'] - R))][0])
                    Mi_real=array(Mi_real)
                       
                    add_params=[]
                    if Ttype=='jeans-Mreal' or Ttype=='alpha-Mreal':
                        Mi = Mi_real
                        add_params=Mi_real
                    else: 
                        Mi=Mi_fit
                    
                    rf_model = prf.inv_M(Mi, pf)
                    
                    Ui = prf.U(rf, pi)
                    Uf = prf.U(rf_model, pf)
                    DU=Uf-Ui
                    
                    # Treal
                    r_real=Treal[k][0]
                    Rvir_real=Treal[k][1]
                    T_real=Treal[k][2]
                    plot(log10(r_real/Rvir_real),T_real/Tviri,'k',linestyle='-',label=r'$T_{real}$')
                    
                    r_real=Treal[k+1][0]
                    Rvir_real=Treal[k+1][1]
                    T_real=Treal[k+1][2]
                    plot(log10(r_real/Rvir_real),T_real/Tviri,'k',linestyle='--',label=r'$T_{real}$')
                    
                    alphaf =ssf[component]['alpha']
                    betaf =ssf[component]['beta_smooth']
                    gammaf=ssf[component]['gamma']
                    m=[]
                    for i in range(size(rf)):
                        mval=get_m(ssf,rf[i])
                        m.append(mval)
                    m=array(m)  
                    
                    if Ttype=='gamma-Treal':
                        Tr=Tr_relevent[index][2]
                        add_params=Tr
                        gammai=1.5*fe.G*Mi/rf/Tr-prf.alpha_Dekel(rf,pi)
                        gammaf=gammai
                        plot(log10(rf/Rviri),Tr/Tviri,'r',linestyle=':',label=r'$T_r$')
                    
                    if Ttype=='Tmulti':
                        M_rmin=max(ssi[component]['r'][0],ssi['all']['r'][0])
                        M_d = ssi[component]['M'][where(ssi[component]['r']>M_rmin)]
                        M_a = ssi['all']['M'][where(ssi[component]['r']>M_rmin)]
                        M_rr=ri[where(ssi[component]['r']>M_rmin)]
                        M_ratio=M_a/M_d

                        slope, intercept,_,_,_ = stats.linregress(log10(M_rr/Rviri),log10(M_ratio))
                        Mratio=10**intercept
                        Mn=-slope
                        add_params=[Mratio,Mn]
                        
                    Ti=get_T(rf,[alphaf,betaf,gammaf,pi],m=0.,Ttype=Ttype,add_params=add_params)
                    Tf=get_T(rf_model,[alphaf,betaf,gammaf,pf],m=m,Ttype=Ttype,add_params=add_params)
                    
                    Tfidi=prf.T_fid(rf,pi)
                    Tfidf=prf.T_fid(rf_model,pf)
                    
                    plot(log10(rf/Rviri),Ti/Tviri,'r',linestyle='-',label=r'$T_i$')
                    plot(log10(rf_model/Rviri),Tf/Tviri,'r',linestyle='--',label=r'$T_f$')

                    ############################################################
                    # COLOR THE PLOTS
                    
                    ax = gca()
                    criteria = t1,    Dsnapshots[k+1],     fmerger[k+1], delta[k+1],delta_i[k+1],delta_if[k+1],fmean[k+1]
                    counts             = [0,0,0,0]
                    if do_color:
                        do_color_subplot(ax,counts,plot_criterion,criteria,criteria_threshold,bgcolor_big_failure=bgcolor_big_failure,bgcolor_big_success =bgcolor_big_success,bgcolor_small_failure=bgcolor_small_failure,bgcolor_small_success=bgcolor_small_success)

                    ############################################################
                    # INDICATE SOME INFORMATION ON THE PLOTS
                    
                    
                    ax.text(0.95,0.9, r'$\rm %i \rightarrow %i$ / '%(k,k+1) +r'${:.01f}$'.format(t1) + r'$\rm \rightarrow$' + r'${:.01f}$'.format(t2) + r'$\rm Gyr$', transform=ax.transAxes, color='k',fontsize=textfont,ha='right')

                    if irow(ncols,nrows)[index]==nrows-1:
                        if icol(ncols,nrows)[index]==0:
                                ax.set_xticks([-2,-1.5,-1,-0.5])
                                ax.set_xticklabels(['-2','-1.5','-1','-0.5'])
                        elif icol(ncols,nrows)[index]==ncols-1:
                                ax.set_xticks([-2,-1.5,-1,-0.5,0])
                                ax.set_xticklabels(['0/-2','-1.5','-1','-0.5','0'])
                        else:
                                ax.set_xticks([-2,-1.5,-1,-0.5])
                                ax.set_xticklabels(['0/-2','-1.5','-1','-0.5'])
                    
                    if icol(ncols,nrows)[index]==0:
                        if irow(ncols,nrows)[index]==0:
                            ax.set_yticks([0,0.5,1.,1.5,2.,2.5,3.])
                            ax.set_yticklabels(['0','0.5','1.0','1.5','2.0','2.5','3.0'])
                        elif irow(ncols,nrows)[index]==nrows-1:
                            ax.set_yticks([-0.5,0.,0.5,1.,1.5,2,2.5,3.])
                            ax.set_yticklabels(['','0','0.5','1.0','1.5','2.0','2.5','0/3.0'])
                        else:
                            ax.set_yticks([0.,0.5,1.,1.5,2,2.5,3.])
                            ax.set_yticklabels(['0','0.5','1.0','1.5','2.0','2.5','0/3.0'])
                    
                    ax.axis(limits_T)
                    lastrow = size(relevent[0])-ncols                    
                    if not(index >= lastrow): ax.set_xticklabels([])
                    if index >= lastrow: xlabel(r'$\log(r/R_{\rm vir})$',fontsize=textfont+2)
                    if not(Decimal(str(index)) % Decimal(str(ncols)) == 0): ax.set_yticklabels([])
                    if Decimal(str(index)) % Decimal(str(ncols)) == 0: ylabel(r'$U,T,T_{fid}$ $[T_{vir}]$')
                    setp(ax.get_yticklabels()[0], visible=False)

    fig.subplots_adjust(hspace=0)
    fig.subplots_adjust(wspace=0) 
    
    if size(title)>0:
        suptitle(title)
    else:
        do_plot_title(sim,relevent,fmean,delta,counts,criteria_threshold)

    show()
    if savefigure:
        fig.savefig(savedir+'/%s_kinetic.pdf'%sim.replace('.',''),bbox_inches='tight')
    
def do_plot_alphabetagamma(sim, relevent,gl,plot_criterion,limits_abg,fitname,Dsnapshots,fmerger,delta,delta_i,delta_if,fmean,criteria_threshold,counts,savefigure,nrows=8,ncols=8,figsize=(20,20),textfont=12,title=[],linear_slopes=False,rlim=[-2,0],savedir='.',do_color=True,component='d',bgcolor_big_failure='moccasin',bgcolor_big_success ='lightgreen',bgcolor_small_failure='lightgoldenrodyellow',bgcolor_small_success='white'):
    [nearly,nochange,nmergers,nsuccess]=counts
    fig = figure(figsize=figsize)
    for index in range(size(relevent[0])): 
        subplot(nrows,ncols,index+1)
        k=relevent[0][index]-1
        
        ssi = gl[k]
        ri = ssi[component]['r']
        pi = ssi[component][fitname]['p']
        (ci, ai, bi, gi, Rviri, Mviri)=pi   
        alphai = ssi[component]['alpha']
        alphai_model = prf.alpha_Dekel(ri,pi)
        betai = ssi[component]['beta_smooth']
        gammai = ssi[component]['gamma']
        
        ssf = gl[k+1]
        rf = ssf[component]['r']
        pf = ssf[component][fitname]['p']
        (cf, af, bf, gf, Rvirf, Mvirf)=pf  
        alphaf = ssf[component]['alpha']
        alphaf_model = prf.alpha_Dekel(rf,pf)
        betaf = ssf[component]['beta_smooth']
        gammaf = ssf[component]['gamma']
        
        if linear_slopes:
            alphai=linearize(alphai,log10(ri/Rviri),rlim)
            betai=linearize(betai,log10(ri/Rviri),rlim)
            gammai=linearize(gammai,log10(ri/Rviri),rlim)
            alphaf=linearize(alphaf,log10(rf/Rvirf),rlim)
            betaf=linearize(betaf,log10(rf/Rvirf),rlim)
            gammaf=linearize(gammaf,log10(rf/Rvirf),rlim)
            
        t1 = ssf['flowdata']['t1']
        t2 = ssf['flowdata']['t2']
            
        plot(log10(ri/Rviri),alphai/4.,'r',linestyle='-',label=r'$\alpha/4$ $(r)$')
        plot(log10(ri/Rviri),alphai_model/4.,'r',linestyle=':')
        plot(log10(rf/Rvirf),alphaf/4.,'r',linestyle='--')
        
        plot(log10(ri/Rviri),betai,'b',linestyle='-',label=r'$\beta$ $(b)$')
        plot(log10(rf/Rvirf),betaf,'b',linestyle='--')
        
        plot(log10(ri/Rviri),gammai,'g',linestyle='-',label=r'$\gamma$ $(g)$')
        plot(log10(rf/Rvirf),gammaf,'g',linestyle='--')
        axhline(0)
        ############################################################
        # COLOR THE PLOTS

        ax = gca()
        criteria = t1,    Dsnapshots[k+1],     fmerger[k+1], delta[k+1],delta_i[k+1],delta_if[k+1],fmean[k+1]
        counts             = [0,0,0,0]
        if do_color:
            do_color_subplot(ax,counts,plot_criterion,criteria,criteria_threshold,bgcolor_big_failure=bgcolor_big_failure,bgcolor_big_success =bgcolor_big_success,bgcolor_small_failure=bgcolor_small_failure,bgcolor_small_success=bgcolor_small_success)

        ############################################################
        # INDICATE SOME INFORMATION ON THE PLOTS

        ax.axis(limits_abg)
        ax.text(0.03,0.9, r'$\rm %i \rightarrow %i$: '%(k,k+1), transform=ax.transAxes, color='k',fontsize=textfont)
        ax.text(0.03,0.83, r'${:.01f}$'.format(t1) + r'$\rm \rightarrow$' + r'${:.01f}$'.format(t2) + r'$\rmGyr$', transform=ax.transAxes, color='k',fontsize=textfont)
        
        lastrow = size(relevent[0])-ncols
        if not(index >= lastrow): ax.set_xticklabels([])
        if index >= lastrow: xlabel(r'$\log(r/R_{\rm vir})$',fontsize=textfont+2)
        if not(Decimal(str(index)) % Decimal(str(ncols)) == 0): ax.set_yticklabels([])
        if Decimal(str(index)) % Decimal(str(ncols)) == 0: ylabel(r'$\alpha/4$ $(r)$, $\beta$ $(b)$, $\gamma$ $(g)$')
        setp(ax.get_yticklabels()[0], visible=False)
        
    fig.subplots_adjust(hspace=0)
    fig.subplots_adjust(wspace=0)

    show()
    if savefigure:
        fig.savefig(savedir+'/%s_slopes.pdf'%sim.replace('.',''),bbox_inches='tight')
        
def do_plot_jeans( relevent,gl,plot_criterion,limits_abg,fitname,Dsnapshots,fmerger,delta,delta_i,delta_if,fmean,criteria_threshold,counts,savefigure,nrows=8,ncols=8,figsize=(20,20),textfont=12,title=[],Ttype='jeans-alpha',linear_slopes=False,rlim=[-2,0],savedir='.',do_color=True,component='d',bgcolor_big_failure='moccasin',bgcolor_big_success ='lightgreen',bgcolor_small_failure='lightgoldenrodyellow',bgcolor_small_success='white'):
    [nearly,nochange,nmergers,nsuccess]=counts
    fig = figure(figsize=figsize)
    for index in range(size(relevent[0])): 
        subplot(nrows,ncols,index+1)
        k=relevent[0][index]-1
        
        ssi = gl[k]
        ri = ssi[component]['r']
        pi = ssi[component][fitname]['p']
        (ci, ai, bi, gi, Rviri, Mviri)=pi   
        alphai = ssi[component]['alpha']
        alphai_model = prf.alpha_Dekel(ri,pi)
        betai = ssi[component]['beta_smooth']
        gammai = ssi[component]['gamma']
             
        ssf = gl[k+1]
        rf = ssf[component]['r']
        pf = ssf[component][fitname]['p']
        (cf, af, bf, gf, Rvirf, Mvirf)=pf  
        alphaf = ssf[component]['alpha']
        alphaf_model = prf.alpha_Dekel(rf,pf)
        betaf = ssf[component]['beta_smooth']
        gammaf = ssf[component]['gamma']

        if linear_slopes:
            alphai=linearize(alphai,log10(ri/Rviri),rlim)
            betai=linearize(betai,log10(ri/Rviri),rlim)
            gammai=linearize(gammai,log10(ri/Rviri),rlim)
            alphaf=linearize(alphaf,log10(rf/Rvirf),rlim)
            betaf=linearize(betaf,log10(rf/Rvirf),rlim)
            gammaf=linearize(gammaf,log10(rf/Rvirf),rlim)
        
        t1 = ssf['flowdata']['t1']
        t2 = ssf['flowdata']['t2']
            
        numerator_i=3.-2.*betai
        denominator_i=get_denominator(ri,[alphai,betai,gammai,pi],Ttype='jeans')
        denominator_model_i=get_denominator(ri,[alphai,betai,gammai,pi],Ttype='jeans-alpha')
        
        plot(log10(ri/Rviri),numerator_i,'r',linestyle='-')
        plot(log10(ri/Rviri),denominator_i,'b',linestyle='-')
        plot(log10(ri/Rviri),denominator_model_i,'b',linestyle=':')
        axhline(0)
        
        ############################################################
        # COLOR THE PLOTS

        ax = gca()
        criteria = t1,    Dsnapshots[k+1],     fmerger[k+1], delta[k+1],delta_i[k+1],delta_if[k+1],fmean[k+1]
        counts             = [0,0,0,0]
        if do_color:
            do_color_subplot(ax,counts,plot_criterion,criteria,criteria_threshold,bgcolor_big_failure=bgcolor_big_failure,bgcolor_big_success =bgcolor_big_success,bgcolor_small_failure=bgcolor_small_failure,bgcolor_small_success=bgcolor_small_success)

        ############################################################
        # INDICATE SOME INFORMATION ON THE PLOTS

        ax.axis([-2,0,-1,6])
        ax.set_xlim([limits_abg[0], limits_abg[1]])
        ax.text(0.03,0.9, r'$\rm %i \rightarrow %i$: '%(k,k+1), transform=ax.transAxes, color='k',fontsize=textfont)
        ax.text(0.03,0.83, r'${:.01f}$'.format(t1) + r'$\rm \rightarrow$' + r'${:.01f}$'.format(t2) + r'$\rmGyr$', transform=ax.transAxes, color='k',fontsize=textfont)
        
        
        lastrow = size(relevent[0])-ncols
        if not(index >= lastrow): ax.set_xticklabels([])
        if index >= lastrow: xlabel(r'$\log(r/R_{\rm vir})$',fontsize=textfont+2)
        if not(Decimal(str(index)) % Decimal(str(ncols)) == 0): ax.set_yticklabels([])
        if Decimal(str(index)) % Decimal(str(ncols)) == 0: ylabel(r'$3-2\beta$ $(r)$, $\alpha+\gamma-2\beta$ $(b)$')
        setp(ax.get_yticklabels()[0], visible=False)
                    
    fig.subplots_adjust(hspace=0)
    fig.subplots_adjust(wspace=0)

    show()
    if savefigure:
        fig.savefig(savedir+'/%s_jeans.pdf'%sim.replace('.',''),bbox_inches='tight')    
    
###############################################################################   
    
# AUXILIARY FUNCTIONS

def do_plot_title(sim,relevent,delta,fmean,counts,criteria_threshold):
    [nearly,nochange,nmergers,nsuccess]=counts
    t_min,Dsnapshot_threshold,merger_thr,Dfit_threshold,fmean_min = criteria_threshold
    suptitle(sim+': '+r'$\rm median$ $\widetilde{|f|}= %.2f$, '%median(abs(array_nonan(fmean[relevent])))
                     +r'$\langle \Delta fit \rangle=$' + r'${:.03f}$, '.format(mean(delta[relevent]))
                     +r'$\rm %i/%i$ '%(nsuccess,size(relevent[0])-nearly-nmergers-nochange)
                     +r'$\rm successes$ '
                     +r'$\rm with$ '
                     +r'$\Delta < %.3f$ '%Dfit_threshold
                     +r'$(%.f \rm percent)$'%(nsuccess/float(size(relevent[0])-nearly-nmergers-nochange)*100.)
             ,fontsize=22,y=0.93)

def do_print_summary(params,print_summary):
    sims,xi,xi_merger,rmax_fit,rmax_evolve,merger_thr,t_min,Dsnapshot_threshold, limits,constrain_fit,constrain_evolution,nsuccess_all,nsample_all = params
    
    import glob,re
    filenumbers=[]
    for filename in glob.glob("/cs/sci/freundlich/CUSPCORE/logs/*.txt"):
        filenumbers.append(int(re.search(r"(?<=log-).*?(?=.txt)", filename).group(0)))
    filename='/cs/sci/freundlich/CUSPCORE/logs/log-%i.txt'%(max(filenumbers)+1)

    file = open(filename,'w') 
    file.write(' \n')
    file.write('SUMMARY\n')
    file.write(' \n')
    file.write('Parameters\n')
    file.write('xi                  = %.2f Rvir sphere size for evaluating f\n'%xi)
    file.write('xi_merger           = %.2f Rvir sphere size for merger identification\n'%xi_merger)
    file.write('rmax_fit            = %.2f Rvir sphere size for fitting\n'%rmax_fit)
    file.write('rmax_evolve         = %.2f Rvir sphere size when evolving\n'%rmax_evolve)
    file.write('merger_thr          = %.2f merger threshold\n'%merger_thr)
    file.write('t_min               = %.2f Gyr minimum time\n'%t_min)
    file.write('Dnsapshot_threshold = %.2f minimum Darea between snapshots\n'%Dsnapshot_threshold)
    file.write('Dfit_threshold      = %.2f minimum Darea for successfull fits\n'%Dsnapshot_threshold)
    file.write('limits              = %.1f %.1f %.1f %.1f\n'%(limits[0],limits[1],limits[2],limits[3]))
    file.write('constraint_fit      = '+str(constrain_fit)+'\n')
    file.write('constraint_evolution= '+str(constrain_evolution)+'\n')

    file.write(' \n')
    file.write('Results\n')
    for i in range(size(sims)):
        percent=nsuccess_all[i]/float(nsample_all[i])*100.
        file.write('%s: %i/%i successes (%.f percent)\n'%(sims[i],nsuccess_all[i],nsample_all[i],percent))

    print ' \n'
    percent_tot=sum(nsuccess_all)/float(sum(nsample_all))*100.
    file.write('TOTAL   : %i/%i successes (%.f percent)\n'%(sum(nsuccess_all),sum(nsample_all),percent_tot))
    file.close()

    print 'A log file has been saved as %s:'%filename
    if print_summary:
        file = open(filename, 'r') 
        print file.read()
    
    return percent_tot
       