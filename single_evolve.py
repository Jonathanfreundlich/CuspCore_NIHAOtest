# PLOT THE PROFILES AND THE MODEL PREDICTION

global ri, pi, brhofit_i, gammai, betai, alphai
global rf, pf, brhofit_f, gammaf, betaf, alphaf
global t1, t2
global Mviri, Mvirf
global m, Mi_fit, Mi_real, Miall_real, Mall, fprofile, wcore,nancore
global Mi_evol,rf_evol, Ui_evol,Uf_evol,Ti_evol,Tf_evol,Ei_evol,Ef_evol
global Darea,Darea_i,Darea_if     

if not multiple_snapshots:
    figure()
else: 
    subplot(nrows,ncols,index+1)

j=1 #1 for 'before', 2 for 'after'
for (i,col,label) in zip((k,k+1),('gray','black'),('before','after')):

    ############################################################
    # RETRIEVE PROFILE FOR EACH OUTPUT
    ss = gl[i]
    a = array(ss['a'])
    t = ss['t']
    r = ss['d']['r']

    brho = ss['d']['brho']
    brhofit = ss['d'][fitname]['brho']
    rho = ss['d']['rho']
    p = ss['d'][fitname]['p']

    Rvir=rvir_fangzhou[i]
    Mvir=mvir_fangzhou[i]
    R12=r12_fangzhou[i]

    if j==1: 
        brhov = brho[-1]
        Rviri = Rvir
        Mviri = Mvir
        ssi = ss

    if j==2:
        Rvirf = Rvir
        Mvirf = Mvir

    ############################################################
    # PLOT SIMULATED PROFILES

    plot(log10(r/Rvir),log10(brho/brhov),color=col,label=label)
    plot(log10(r/Rvir),log10(brhofit/brhov),'--',color=col)  

    if j == 1:
        ri = r
        pi = p
        brhofit_i = brhofit
        gammai=ss['d']['gamma']
        betai =ss['d']['beta_smooth']
        alphai =ss['d']['alpha']


    if j == 2:
        rf = r
        pf = p
        brhofit_f = brhofit
        gammaf=ss['d']['gamma']
        betaf =ss['d']['beta_smooth'] 
        alphaf =ss['d']['alpha']
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
        Miall_real=[]
        for i in range(size(rf)):
            R=rf[i]
            Mi_real.append(ssi['d']['M'][abs(ssi['d']['r'] - R) == min(abs(ssi['d']['r'] - R))][0])
            Miall_real.append(ssi['d']['Mall'][abs(ssi['d']['r'] - R) == min(abs(ssi['d']['r'] - R))][0])
        Mi_real=array(Mi_real)
        Miall_real=array(Miall_real)
        Mall = ss['d']['Mall']
        
        add_params=[]
        if Ttype[-5:]=='Mreal':
            Mi = Miall_real
            add_params=Miall_real
        else: 
            Mi=Mi_fit
                            
        fprofile=m/Mi

        #####################################################
        # EVOLVE PROFILE

        wcore = 0*(r<=rmin_evolve*Rvirf)+1*((r>rmin_evolve*Rvirf)&(r<=rmax_evolve*Rvirf))+0*(r>rmax_evolve*Rvirf)
        nancore=nan*ones_like(wcore)
        nancore[where(wcore==1)]=1

        if Ttype=='gamma-Treal':
            add_params=Treal_relevent[index][2]

        if Ttype=='Tmulti':
            component='d'
            M_rmin=max(ss[component]['r'][0],ss['all']['r'][0])
            M_d = ss[component]['M'][where(ss[component]['r']>M_rmin)]
            M_a = ss['all']['M'][where(ss[component]['r']>M_rmin)]
            M_rr=r[where(ss[component]['r']>M_rmin)]
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
        
        if Ptype=='f=i':
            print 'Ptype = ', Ptype
            alphaf_evol=alphai_evol
            betaf_evol=betai_evol
            gammaf_evol=gammai_evol
        elif Ptype=='b=0':
            betai_evol=0.
            betaf_evol=0.

        if constrain_evolution:
            exec("constraint=({'type': 'ineq', 'fun': lambda x:  x[1]+3*sqrt(x[0]*%.12f)})"%(rmin_evolve))
            res = evolve_constrained(rf, rf, Mi, pi, m, constraint,alphai_evol,alphaf_evol,gammai_evol, gammaf_evol,betai_evol,betaf_evol, Ttype=Ttype,w=wcore, method='halo',add_params=add_params)
        else:
            res = evolve(rf, rf, Mi, pi, m, alphai_evol,alphaf_evol,gammai_evol, gammaf_evol,betai_evol,betaf_evol, Ttype=Ttype, w=wcore, method='halo',add_params=add_params)

        if Ttype=='jeans-gamma':
            print 'gamma = ', res['gammaf']

        p_model = res['pf']
        brhofit_model = prf.brho(rf,res['pf'])
        rhofit_model = prf.rho(rf,res['pf'])
        ss.update({'evolution_rms':fit.rmslog(brhofit_f,brhofit_model)})

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

        plot(log10(rf/Rvir),log10(brhofit_model/brhov),'--',color='red',label='toy model')

    j+=1

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
                
ax = gca()
ax.axis(limits)
ax.text(0.03,0.9, r'$\rm %i \rightarrow %i$: '%(k,k+1), transform=ax.transAxes, color=col,fontsize=textfont)
ax.text(0.03,0.83, r'${:.01f}$'.format(t1) + r'$\rm \rightarrow$' + r'${:.01f}$'.format(t2) + r'$\rmGyr$', transform=ax.transAxes, color=col,fontsize=textfont)
ax.text(0.03,0.20, r'$\Delta snap={:.2f}$'.format(Dsnapshots[k+1]), fontsize = textfont, transform=ax.transAxes,color=col)
ax.text(0.03,0.13, r'$|f|={:.2f} \pm {:.2f} $'.format(fmean[k+1],fstd[k+1]), fontsize = textfont, transform=ax.transAxes,color=col)
ax.text(0.03,0.06, r'$\Delta fit={:.2f}$'.format(Darea), fontsize = textfont, transform=ax.transAxes,color=col)

xlabel(r'$\log(r/R_{\rmvir})$',fontsize=textfont+2)
ylabel(r'$\log(\bar{\rho}/\bar{\rho}_{\rmvir})$',fontsize=textfont+2)

############################################################
# COLOR THE PLOTS AND COUNT THE NUMBER OF SUCCESSES

criteria           = t1,    Dsnapshots[k+1],     fmerger[k+1], Darea,Darea_i,Darea_if, fmean[k+1]
criteria_threshold = t_min, Dsnapshot_threshold, merger_thr,   Dfit_threshold,fmean_min
counts             = [nearly,nochange,nmergers,nsuccess]
counts,col,bgcolor = cuspcore_plots.do_color_subplot(ax,counts,plot_criterion,criteria,criteria_threshold)
[nearly,nochange,nmergers,nsuccess] = counts

############################################################
# FORMAT IF multiple_snapshots

if multiple_snapshots:# and size(relevent[0])<>0:
    lastrow = size(relevent[0])-ncols #round(size(relevent[0])/ncols-1)*nrows
    if not(index >= lastrow): ax.set_xticklabels([])
    if index >= lastrow: xlabel(r'$\log(r/R_{\rmvir})$',fontsize=textfont+2)
    if not(Decimal(str(index)) % Decimal(str(ncols)) == 0): ax.set_yticklabels([])
    if Decimal(str(index)) % Decimal(str(ncols)) == 0: ylabel(r'$\log(\bar{\rho}/\bar{\rho}_{\rmvir})$',fontsize=textfont+2)
    setp(ax.get_yticklabels()[0], visible=False)
    xticks(xticks_val)#,rotation=90)

############################################################

print 'k = %i - %i: delta = %.3f'%(k,k+1,Darea)
print 'p_init : c=%.2f, a=%.2f, 1/b=%.2f, g=%.2f, Rvir=%.2f kpc, Mvir=%.2e Msun'%pi
print 'p_final: c=%.2f, a=%.2f, 1/b=%.2f, g=%.2f, Rvir=%.2f kpc, Mvir=%.2e Msun'%pf
print 'p_model: c=%.2f, a=%.2f, 1/b=%.2f, g=%.2f, Rvir=%.2f kpc, Mvir=%.2e Msun'%p_model
print ' '
    
    