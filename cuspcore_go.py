###############################################################################
#            TESTING THE MODEL FOR THE CUSP-CORE TRANSFORMATION               #
###############################################################################

import os
import sys
from numpy import *
from matplotlib.pylab import *
from decimal import Decimal
import pickle
import inspect
from scipy import stats

import treal_functions
import general_functions
import slopes_functions
import prepare_functions
import fitting as fit
import cuspcore_aux
import cuspcore_plots
import evolving

from general_functions import *
from slopes_functions import *
from prepare_functions import *
from cuspcore_aux import *
from cuspcore_plots import *
from evolving import *

###############################################################################

def go(sims,
       directory='/cs/sci/freundlich/CUSPCORE/NIHAO_data/',
       # PARAMETERS
       xi=0.01,
       xi_merger=0.1,
       rmax_fit=1.,
       rmin_fit=0.01,
       rmax_evolve=1.,
       rmin_evolve=0.01,
       merger_thr=0.15,
       f_min=0.10,
       fmean_min=0.06,
       t_min=2.,
       Dsnapshot_threshold=0.07,
       Dfit_threshold=0.07,
       m_constraint=2.,
       component='d',
       delta_xlim=[-2,-1.5],
       delta_ymin=2.,
       # OPTIONS
       plotting_option='tmin-crit-nomergers',
       plot_criterion='f',
       constrain_fit=False,
       constrain_evolution=False,
       Ttype='jeans',
       # DISPLAY/PLOTTING OPTIONS
       selection=[],
       nrows=8,
       ncols=8,
       figsize=(20,20),
       textfont=12,
       title=[],
       plot_fprofile=True,
       plot_T=False,
       plot_alphabetagamma=False,
       test_best=False,
       print_summary=True,
       limits=[-2,-1,2,3],
       limits_fprofile=[-2,0.,0,0.5],
       limits_T=(-2,0,0,4),
       limits_abg=(-2,0,-1,1),
       xticks_val=(-2,-1.5,-1),
       xticks_fprofile=(-2,-1.5,-1,-0.5,0),
       savefigure=False,
       savedir='.',
       do_color=True,
       # SMOOTHING OPTIONS FOR ALPHABETAGAMMA
       polyorder=3,
       sigma = 21,
       mode= 'interp',
       double_smooth=False,
       linear_slopes=False,
       multi=False,
       use_RMvir=False,
       bgcolor_big_failure='moccasin',
       bgcolor_big_success ='lightgreen',
       bgcolor_small_failure='lightgoldenrodyellow',
       bgcolor_small_success='white',
       imin=0,
       return_p=False):
    
    # PRINT ALL PARAMETERS OF THE FUNCTION
    frame = inspect.currentframe()
    args, _, _, values = inspect.getargvalues(frame)
    print ' '
    print 'Parameters of function %s:' % inspect.getframeinfo(frame)[2]
    for i in args:
        print "    %s = %s" % (i, values[i])
    print ' '
    
    # CHECK FUNCTION PARAMETERS
    if constrain_fit:
        fitname = 'lsfit_brho_b2_g3_constrained'
    else: 
        fitname = 'lsfit_brho_b2_g3_unconstrained'

    if constrain_fit and constrain_evolution:
        filestart='constrained'
    elif not constrain_fit and not constrain_evolution:
        filestart='unconstrained'
    else:
        filestart='mixed'
        
    if test_best and not plot_fprofile:
        plot_fprofile = True
        print '    plot_fprofile was changed to True to allow test_best'
        print ' ' 
    
    rlim=[log10(rmin_fit),log10(rmax_fit)]
    
    # INITIALIZATION 
    delta_all=[]
    delta_i_all=[]
    Dsnapshot_all=[]
    fmean_all=[]
    fstd_all=[]
    fmerger_all=[]
    
    #f_all=[]
    f_tot=[]
    delta_if_all=[]
    delta_tot=[]
    nsuccess_all=[]
    nrelevent_all=[]
    nmergers_all=[]
    nochange_all=[]
    nsample_all=[]
    nearly_all=[]
    fprofile_max=[]
    #xi_all=[]
    
    relevent_fangzhou_all=[]
    relevent_all=[]
    t1_array_all=[]
    t2_array_all=[]
    fmerger_all=[]

    fbest_all=[]
    Rbest_all=[]
    xbest_all=[]
    rprofile_all=[]
    xprofile_all=[]
    
    fprofile_all=[]
    fprofile_simname=[]
    mprofile_all=[]
    
    R12_all=[]
    Rvir_all=[]
    Mvir_all=[]
    xi_all=[]
    ok_fangzhou_all=[]
    Rvir_halo_all=[]
    Mvir_halo_all=[]
    Mstar_all=[]
    

    R12_coldgas_all=[]
    R12_coldbar_all=[]
    R12_SFR_all=[]
    Rs_NFW_all=[]
    R2_Einasto_all=[]
    Rc_Dekel_all=[]
    Re_Sersic_star_all=[]
    
    pi_all=[]
    pf_all=[]
    pmodel_all=[]

    for sim in sims:
        
        ############################################################
        # LOAD DATA 
        
        print 'Simulation %s'%sim
        with open(directory+'NIHAO-%s.pickle'%sim[1:]) as f:
            gl = pickle.load(f)
        
        gl=gl[:]
        gl[imin:] = slopes_functions.derive_slopes( gl[imin:],polyorder=polyorder,sigma=sigma,mode=mode,double_smooth=double_smooth,rlim=rlim,use_fangzhou_Rvir=True,linearize=True,betanull=False)
        gl[imin:] = prepare_functions.define_brho( gl[imin:],polyorder=polyorder,sigma=sigma,mode=mode,double_smooth=double_smooth,rlim=rlim,use_fangzhou_Rvir=True)
        
        treal= treal_functions.load_or_create_gl(sim,directory=directory,use_fangzhou_Rvir=True,c=component)
       
        
        fitrange=prepare_functions.get_fitrange(gl, use_fangzhou_Rvir=True,component=component)
        fitrange={'all':fitrange['all'][imin:],'d':fitrange['d'][imin:],'s':fitrange['s'][imin:],'g':fitrange['g'][imin:]}
        gl[imin:]=prepare_functions.reduce_range_gl(gl[imin:],fitrange)
        treal[imin:]=prepare_functions.reduce_range_Treal(treal[imin:],fitrange)
        
        ############################################################
        # GET TIME, RVIR AND MVIR FOR EACH OUTPUT

        a_array=[]
        t1_array=[nan]
        t2_array=[nan]
        rvir_halo=[]
        mvir_halo=[]
        
        for (i,ss) in zip(range(size(gl)),gl):
            a_array.append(ss['a'])
            rvir_halo.append(ss['all']['Rvir'])
            mvir_halo.append(ss['all']['Mvir'])
            if i>0:  
                try:
                    t1_array.append(ss['flowdata']['t1'])
                    t2_array.append(ss['flowdata']['t2'])       
                except:
                    t1_array.append(nan)
                    t2_array.append(nan)            
        
        a_array=array(a_array)
        t1_array=array(t1_array)
        t2_array=array(t2_array)
        rvir_halo=array(rvir_halo)
        mvir_halo=array(mvir_halo)
        
        ############################################################
        # GET RADII FROM FANGZHOU'S FILES
          ok_fangzhou,r12_fangzhou,rvir_fangzhou,mvir_fangzhou,mstar_fangzhou=get_fangzhou_radii(sim,a_array,get_stars=True)
        
        Mstar_all.append(mstar_fangzhou)
        ok_fangzhou_all.append(ok_fangzhou)
        R12_all.append(r12_fangzhou)
        Rvir_all.append(rvir_fangzhou)
        Mvir_all.append(mvir_fangzhou)
        Rvir_halo_all.append(rvir_halo)
        Mvir_halo_all.append(mvir_halo)
        xi_all.append(r12_fangzhou/rvir_fangzhou)
        
        ############################################################
        # CARRY OUT FITS FOR EACH SNAPSHOTS
        
        if component=='s':
            gl[imin:]=fit.do_fits(gl[imin:], r12_fangzhou[imin:], mstar_fangzhou[imin:]/2., rmax_fit,rmin_fit,m_constraint,components=[component])
        else:
            gl[imin:]=fit.do_fits(gl[imin:], rvir_fangzhou[imin:], mvir_fangzhou[imin:], rmax_fit, rmin_fit, m_constraint,components=[component])
        print 'Fitting done'
           
        ############################################################
        # DEFINE DIFFERENT QUANTITIES PER OUTPUT
        
        glsim = [ss for ss in gl[1:]]
        fmerger = array([nan]+[get_f(gl,ss, xi_merger,Rtype='ratio',com='all') for ss in glsim])
        previous_fmerger=array([nan]+fmerger[:-1].tolist())
        
        rms = [nan]
        for ss in glsim:
            if 'rms' in ss['all'].keys(): rms.append(ss['all']['rms'])
            else: rms.append(nan)
        rms=array(rms)

        Dsnapshots=ones(size(gl))
        Dsnapshots[0]=nan
        for i in range(size(glsim)):
            try:
                rfit = gl[i+1][component]['r']
                Rvir = rvir_fangzhou[i+1]
                rr=rfit[(rfit>=rmin_fit*Rvir)&(rfit<=rmax_fit*Rvir)]
                if size(gl[i][component]['brho'])>0:
                    brhov = gl[i][component]['brho'][0]
                else:
                    brhov = nan
                p1=gl[i][component][fitname]['p']
                p2=gl[i+1][component][fitname]['p']
                fit1=prf.brho(rr, p1, model='an')
                fit2=prf.brho(rr, p2, model='an')

                Dsnapshots[i+1]= fit.Delta_area(log10(rr/Rvir),log10(fit1/brhov), log10(fit2/brhov), xlimits=delta_xlim,ymin=delta_ymin)
            except:
                Dsnapshots[i+1]=nan

        relevent_fangzhou=ones(size(rvir_fangzhou),dtype=bool)
        for i in range(size(rvir_fangzhou)):
            if isnan(rvir_fangzhou[i]):
                relevent_fangzhou[i]=False
                if i<>size(rvir_fangzhou)-1:relevent_fangzhou[i+1]=False

        delta=nan*ones(size(gl))
        delta_i=nan*ones(size(gl))
        delta_if=nan*ones(size(gl))
        fmean=nan*ones(size(gl))
        fstd=nan*ones(size(gl))
        
        pi_array=[(nan, nan, nan, nan, nan, nan)]*size(gl)
        pf_array=[(nan, nan, nan, nan, nan, nan)]*size(gl)
        pmodel_array=[(nan, nan, nan, nan, nan, nan)]*size(gl)

        ############################################################
        # CHOOSE PLOTTING OPTION
        
        relevent=where((relevent_fangzhou) & (t1_array>t_min) & (abs(fmerger) < merger_thr) & (abs(previous_fmerger) < merger_thr))
        print 'Number of relevent snapshots: ',size(relevent)
        
        if size(selection)>0:
            relevent=(array(selection),)
            
        ############################################################
        # APPLY MODEL AND PLOT PROFILES
        
        rmss=[]
        fs=[]
        pmodel=[None] * size(gl)
        pmodel_T=[None] * size(gl)
        
        nsuccess=0
        nmergers=0
        nochange=0
        nearly=0
        
        fig = figure(figsize=figsize)
        
        fprofile_relevent=[]
        mprofile_relevent=[]
        Mi_relevent=[]
        rprofile_relevent=[]
        xprofile_relevent=[]
        Treal_relevent=[]
        
        for index in range(size(relevent[0])):
                subplot(nrows,ncols,index+1)
                k=relevent[0][index]-1
                j=1 #1 for 'before', 2 for 'after'
                for (i,col,label) in zip((k,k+1),('gray','black'),('before','after')):

                    ############################################################
                    # RETRIEVE PROFILE FOR EACH OUTPUT
                    ss = gl[i]
                    a = array(ss['a'])
                    t = ss['t']
                    r = ss[component]['r']
    
                    brho = ss[component]['brho']
                    brhofit = ss[component][fitname]['brho']
                    rho = ss[component]['rho']
                    p = ss[component][fitname]['p']
                    
                    Rvir=rvir_fangzhou[i]
                    Mvir=mvir_fangzhou[i]
                    R12=r12_fangzhou[i]
                    Mstar=mstar_fangzhou[i]

                    if j==1: 
                        brhov = brho[0]
                        Rviri = Rvir
                        Mviri = Mvir
                        Rstari= R12
                        Mstari= Mstar
                        ssi = ss
                    
                    if j==2:
                        Rvirf = Rvir
                        Mvirf = Mvir
                        Rstarf= R12
                        Mstarf= Mstar

                    ############################################################
                    # PLOT SIMULATED PROFILES
                    
                    plot(log10(r/Rvir),log10(brho/brhov),color=col,label=label)
                    plot(log10(r/Rvir),log10(brhofit/brhov),'--',color=col,lw=2)  

                    if j == 1:
                        ri = r
                        pi = p
                        brhofit_i = brhofit
                        gammai=ss[component]['gamma']
                        betai =ss[component]['beta_smooth']
                        alphai =ss[component]['alpha']
                        
                        
                    if j == 2:
                        rf = r
                        pf = p
                        brhofit_f = brhofit
                        gammaf=ss[component]['gamma']
                        betaf =ss[component]['beta_smooth'] 
                        alphaf =ss[component]['alpha']
                        t1 = ss['flowdata']['t1']
                        t2 = ss['flowdata']['t2']
                        
                        ############################################################
                        # DEFINE m and f
                        
                        m=[]
                        for i in range(size(rf)):
                            mval=get_m(ss,rf[i])
                            m.append(mval)
                        m=array(m)

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
                        
                        fprofile=m/Mi 
                        
                        #print 'Mi_real=', Mi_real
                        ############################################################

                        mprofile_relevent.append(m)
                        rprofile_relevent.append(rf)
                        xprofile_relevent.append(rf/Rvir)
                        fprofile_relevent.append(fprofile)
                        Mi_relevent.append(Mi)
                        
                        Treal_relevent.append(treal[k][:3])
                        
                        #####################################################
                        # EVOLVE PROFILE
                        
                        if component=='s':
                            wcore = 0*(r<=rmin_evolve*Rvirf)+1*((r>rmin_evolve*Rvirf)&(r<=rmax_evolve*Rstarf))+0*(r>rmax_evolve*Rstarf)
                        else:
                            wcore = 0*(r<=rmin_evolve*Rvirf)+1*((r>rmin_evolve*Rvirf)&(r<=rmax_evolve*Rvirf))+0*(r>rmax_evolve*Rvirf)
                        nancore=nan*ones_like(wcore)
                        nancore[where(wcore==1)]=1
                        
                        if Ttype=='gamma-Treal':
                            add_params=Treal_relevent[index][2]
                           
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
                           
                        if linear_slopes:
                            alphai_evol=linearized_function(log10(rf/Rvir),log10(ri/Rviri),alphai,rlim)
                            betai_evol=linearized_function(log10(rf/Rvir),log10(ri/Rviri),betai,rlim)
                            gammai_evol=linearized_function(log10(rf/Rvir),log10(ri/Rviri),gammai,rlim)

                            alphaf_evol=linearized_function(log10(rf/Rvir),log10(rf/Rvir),alphaf,rlim)
                            betaf_evol=linearized_function(log10(rf/Rvir),log10(rf/Rvir),betaf,rlim)
                            gammaf_evol=linearized_function(log10(rf/Rvir),log10(rf/Rvir),gammaf,rlim)                   
                        else:
                            alphai_evol=alphaf
                            betai_evol=betaf
                            gammai_evol=gammaf

                            alphaf_evol=alphaf
                            betaf_evol=betaf
                            gammaf_evol=gammaf
                            
                        if constrain_evolution:
                            exec("constraint=({'type': 'ineq', 'fun': lambda x:  x[1]+3*sqrt(x[0]*%.12f/%.12f)})"%(rmin,float(Rvir)))
                            res = evolve_constrained(rf, rf, Mi, pi, m, constraint,alphai_evol,alphaf_evol,gammai_evol, gammaf_evol,betai_evol,betaf_evol, Ttype=Ttype,w=wcore, method='halo',add_params=add_params)

                        else:
                            if use_RMvir:
                                res = evolve(rf, rf, Mi, pi, m, alphai_evol,alphaf_evol,gammai_evol, gammaf_evol,betai_evol,betaf_evol, Ttype=Ttype, w=wcore, method='halo',add_params=add_params,Rvirf=Rvirf,Mvirf=Mvirf)
                            else:
                                res = evolve(rf, rf, Mi, pi, m, alphai_evol,alphaf_evol,gammai_evol, gammaf_evol,betai_evol,betaf_evol, Ttype=Ttype, w=wcore, method='halo',add_params=add_params)

                        if Ttype=='jeans-gamma':
                            print 'gamma = ', res['gammaf']
                        
                        p_model = res['pf']
                        brhofit_model = prf.brho(rf,res['pf'])
                        rhofit_model = prf.rho(rf,res['pf'])
                        ss.update({'evolution_rms':fit.rmslog(brhofit_f,brhofit_model)})
                        rmss.append(fit.rmslog(brhofit_f,brhofit_model))

                        ############################################################
                        # QUANTITIES USED IN THE EVOLUTION

                        Mi_evol = prf.M(rf, pi, "an")
                        rf_evol = prf.inv_M(Mi_evol, p_model, "an") #the new radii that enclose the same masses
                        Ui_evol = prf.U(rf, pi, "an")
                        Uf_evol = prf.U(rf, p_model, "an")
                        Ti_evol=get_T(rf,[alphai_evol,betai_evol,gammai_evol,pi],m=0.,Ttype=Ttype,do_smooth=False,add_params=add_params)
                        Tf_evol=get_T(rf,[alphaf_evol,betaf_evol,gammaf_evol,pf],m=m ,Ttype=Ttype,do_smooth=False,add_params=add_params)
                        Ei_evol = Ui_evol - G*m/rf + Ti_evol #0.5*G*Mi/ri
                        Ef_evol = Uf_evol - G*m/rf_evol + Tf_evol #0.5*G*Mi/rf 

                        ############################################################
                        # PLOT PREDICTED PROFILE
                        
                        plot(log10(rf/Rvir),log10(brhofit_model/brhov),'-',color='red',label='toy model',lw=2)

                    j+=1
                    
                ############################################################
                # ASSESS THE PREDICTION SUCCESS
                
                maxr = rmax_fit*Rvirf
                rr=rf[rf<=maxr]
                Darea=fit.Delta_area(log10(rr/Rvirf),log10(prf.brho(rr, pf)/brhov),log10(prf.brho(rr, p_model)/brhov),xlimits=delta_xlim,ymin=delta_ymin)
                delta[k+1]=Darea
                
                Darea_i=fit.Delta_area(log10(rr/Rvir),log10(prf.brho(rr, pi)/brhov),log10(prf.brho(rr, p_model)/brhov),xlimits=delta_xlim,ymin=delta_ymin)
                delta_i[k+1]=Darea_i
                
                Darea_if=fit.Delta_area(log10(rr/Rvir),log10(prf.brho(rr, pi)/brhov),log10(prf.brho(rr, pf)/brhov),xlimits=delta_xlim,ymin=delta_ymin)
                delta_if[k+1]=Darea_if     
                
                maxr = rmax_evolve*Rvirf
                fmean[k+1] = mean(abs(array_nonan(fprofile[rf<=maxr])))
                fstd[k+1] = std(abs(array_nonan(fprofile[rf<=maxr])))

                pi_array[k+1]=pi
                pf_array[k+1]=pf
                pmodel_array[k+1]=p_model
                
                ############################################################
                # COLOR THE PLOTS AND COUNT THE NUMBER OF SUCCESSES

                ax = gca()
                criteria           = t1,    Dsnapshots[k+1],     fmerger[k+1], Darea,Darea_i,Darea_if, fmean[k+1]
                criteria_threshold = t_min, Dsnapshot_threshold, merger_thr,   Dfit_threshold,fmean_min
                counts             = [nearly,nochange,nmergers,nsuccess]
                if do_color:
                    counts,col,bgcolor = cuspcore_plots.do_color_subplot(ax,counts,plot_criterion,criteria,criteria_threshold,bgcolor_big_failure=bgcolor_big_failure,bgcolor_big_success =bgcolor_big_success,bgcolor_small_failure=bgcolor_small_failure,bgcolor_small_success=bgcolor_small_success)
                    [nearly,nochange,nmergers,nsuccess] = counts

                ############################################################
                # INDICATE SOME INFORMATION ON THE PLOTS

                ax.text(0.95,0.9, r'$\rm %i \rightarrow %i$ / '%(k,k+1) +r'${:.01f}$'.format(t1) + r'$\rm \rightarrow$' + r'${:.01f}$'.format(t2) + r'$\rmGyr$', transform=ax.transAxes, color=col,fontsize=textfont,ha='right')
                ax.text(0.03,0.27, r'$f_{\rm merger}=%.2f$'%fmerger[k+1], fontsize = textfont, transform=ax.transAxes,color=col)
                ax.text(0.03,0.20, r'$\delta_{\rm sim}=%.2f$'%Darea_if, fontsize = textfont, transform=ax.transAxes,color=col)
                ax.text(0.03,0.13, r'$\delta_0=%.2f$'%Darea_i, fontsize = textfont, transform=ax.transAxes,color=col)
                ax.text(0.03,0.06, r'$\delta =%.2f$'%Darea, fontsize = textfont, transform=ax.transAxes,color=col)

                if icol(ncols,nrows)[index]==0:
                    ax.set_xticks([-2,-1.5])
                    ax.set_xticklabels([r'$-2$',r'$-1.5$'], fontdict={'fontsize':textfont})
                elif icol(ncols,nrows)[index]==ncols-1:
                    ax.set_xticks([-2,-1.5,-1])
                    ax.set_xticklabels([r'$-1/-2$',r'$-1.5$',r'$-1$'], fontdict={'fontsize':textfont})
                else:
                    ax.set_xticks([-2,-1.5])
                    ax.set_xticklabels([r'$-1/-2$',r'$-1.5$'], fontdict={'fontsize':textfont})
                 
                if irow(ncols,nrows)[index]==0:
                    ax.set_yticks([1.8,2.,2.2,2.4,2.6,2.8,3,3.2])
                    ax.set_yticklabels([r'',r'$2.0/3.2$',r'$2.2$',r'$2.4$',r'$2.6$',r'$2.8$',r'$3.0$',r'$3.2$'],
                                       fontdict={'fontsize':textfont})
                elif irow(ncols,nrows)[index]==nrows-1:
                    ax.set_yticks([1.8,2.,2.2,2.4,2.6,2.8,3])
                    ax.set_yticklabels([r'',r'$2.0$',r'$2.2$',r'$2.4$',r'$2.6$',r'$2.8$',r'$3.0$'],
                                       fontdict={'fontsize':textfont})
                else:
                    ax.set_yticks([1.8,2.,2.2,2.4,2.6,2.8,3])
                    ax.set_yticklabels([r'',r'$2.0/3.2$',r'$2.2$',r'$2.4$',r'$2.6$',r'$2.8$',r'$3.0$'], 
                                       fontdict={'fontsize':textfont})
                
                ax.axis(limits)
                lastrow = size(relevent[0])-ncols 
                if not(index >= lastrow): ax.set_xticklabels([])
                if index >= lastrow: xlabel(r'$\log(r/R_{\rmvir})$',fontsize=textfont+2)
                if not(Decimal(str(index)) % Decimal(str(ncols)) == 0): ax.set_yticklabels([])
                if Decimal(str(index)) % Decimal(str(ncols)) == 0: ylabel(r'$\log(\bar{\rho}/\bar{\rho}_{\rmvir})$',fontsize=textfont+2)
                setp(ax.get_yticklabels()[0], visible=False)
                xticks(xticks_val)

        fig.subplots_adjust(hspace=0)
        fig.subplots_adjust(wspace=0) 
        
        if size(title)>0:
            suptitle(title)
        else:
            cuspcore_plots.do_plot_title(sim,relevent,delta,fmean,counts,criteria_threshold)
        
        show()
        
        if savefigure:
            fig.savefig(savedir+'/%s_density.pdf'%sim.replace('.',''),bbox_inches='tight')

        ############################################################
        # DEFINE SOME QUANTITIES FOR ALL SIMULATIONS
        
        relevent_fangzhou_all.append(relevent_fangzhou)
        relevent_all.append(relevent)
        t1_array_all.append(t1_array)
        t2_array_all.append(t2_array)
        fmerger_all.append(fmerger)
        Dsnapshot_all.append(Dsnapshots)
        delta_all.append(array(delta))
        delta_i_all.append(array(delta_i))
        delta_if_all.append(array(delta_if))
        nearly_all.append(nearly)
        nsuccess_all.append(nsuccess)
        nrelevent_all.append(size(relevent[0]))
        nmergers_all.append(nmergers)
        nochange_all.append(nochange)
        nsample_all.append(size(relevent[0])-nearly-nmergers-nochange)
        fmean_all.append(array(fmean))
        fstd_all.append(array(fstd))
        
        pi_all.append(pi_array)
        pf_all.append(pf_array)
        pmodel_all.append(pmodel_array)
        
        ##################################################################################################################
        # PLOT F PROFILE
        
        if plot_fprofile:
            criteria_threshold = t_min, Dsnapshot_threshold, merger_thr,   Dfit_threshold, fmean_min
            counts=[nearly,0,nmergers,nsuccess]
            params=sim,xi,xi_merger,rmax_fit,rmax_evolve,merger_thr,t_min,Dsnapshot_threshold,\
    limits,constrain_fit,constrain_evolution,0,0
    
            cuspcore_plots.do_plot_fprofile( relevent,xprofile_relevent,fprofile_relevent,plot_criterion,limits_fprofile,fitname,Dsnapshots,fmerger,delta,delta_i,delta_if,fmean,fstd,criteria_threshold,counts,savefigure,params,t1_array,t2_array,nrows=nrows,ncols=ncols,figsize=figsize,textfont=textfont,title=title,savedir=savedir,do_color=do_color,bgcolor_big_failure=bgcolor_big_failure,bgcolor_big_success =bgcolor_big_success,bgcolor_small_failure=bgcolor_small_failure,bgcolor_small_success=bgcolor_small_success)
                
        ############################################################
        # PLOT T PROFILE
        
        if plot_T:
            criteria_threshold = t_min, Dsnapshot_threshold, merger_thr,   Dfit_threshold, fmean_min
            counts=[nearly,0,nmergers,nsuccess]
            cuspcore_plots.do_plot_T(        sim,relevent,gl,treal,plot_criterion,limits_T,fitname,Dsnapshots,fmerger,delta,delta_i,delta_if,fmean,criteria_threshold,counts,savefigure,nrows=nrows,ncols=ncols,figsize=figsize,textfont=textfont,title=title,Ttype=Ttype,savedir=savedir,do_color=do_color,component=component,bgcolor_big_failure=bgcolor_big_failure,bgcolor_big_success =bgcolor_big_success,bgcolor_small_failure=bgcolor_small_failure,bgcolor_small_success=bgcolor_small_success)
            
        ############################################################
        # PLOT alpha, beta, gamma PROFILE
        
        if plot_alphabetagamma:
            criteria_threshold = t_min, Dsnapshot_threshold, merger_thr,   Dfit_threshold, fmean_min
            counts=[nearly,0,nmergers,nsuccess]
            cuspcore_plots.do_plot_alphabetagamma( sim,relevent,gl,plot_criterion,limits_abg,fitname,Dsnapshots,fmerger,delta,delta_i,delta_if,fmean,criteria_threshold,counts,savefigure,nrows=nrows,ncols=ncols,figsize=figsize,textfont=textfont,title=title,linear_slopes=linear_slopes,rlim=rlim,savedir=savedir,do_color=do_color,component=component)
            cuspcore_plots.do_plot_jeans( relevent,gl,plot_criterion,limits_abg,fitname,Dsnapshots,fmerger,delta,delta_i,delta_if,fmean,criteria_threshold,counts,savefigure,nrows=nrows,ncols=ncols,figsize=figsize,textfont=textfont,title=title,linear_slopes=linear_slopes,rlim=rlim,savedir=savedir,do_color=do_color,component=component,bgcolor_big_failure=bgcolor_big_failure,bgcolor_big_success =bgcolor_big_success,bgcolor_small_failure=bgcolor_small_failure,bgcolor_small_success=bgcolor_small_success)
                
    ############################################################
    # PRINT SUMMARY
    params=sims,xi,xi_merger,rmax_fit,rmax_evolve,merger_thr,t_min,Dsnapshot_threshold,\
    limits,constrain_fit,constrain_evolution,nsuccess_all,nsample_all
    percent_tot=cuspcore_plots.do_print_summary(params,print_summary)
    
    ############################################################
    ### BEEP WHEN DONE
    do_beep()
    
    if return_p:
        return [t1_array_all, t2_array_all, Dsnapshot_all, delta_all, delta_i_all, fmean_all, fmerger_all, pi_all, pf_all, pmodel_all]
    else:
        return [percent_tot, t1_array_all, t2_array_all, Dsnapshot_all, delta_all, delta_i_all, fmean_all, fstd_all, fmerger_all, fprofile_all, fprofile_simname, criteria_threshold]



