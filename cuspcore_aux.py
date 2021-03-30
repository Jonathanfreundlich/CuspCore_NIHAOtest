###############################################################################
#               AUXILIARY FUNCTIONS FOR CUSP-CORE MODELLING                   #
###############################################################################

import sys
from numpy import *
from matplotlib.pylab import *

import format_functions
import treal_functions
import general_functions
from general_functions import *
import fitting as fit

from scipy.signal import savgol_filter
from lmfit import *
from IPython.display import Audio,display

components=['all', 'd', 's', 'g']

###############################################################################

# GET THE MASS FLOWING INSIDE A SHELL AND THE CORRESPONDING FRACTION

def get_m(ss,R,dr='net',Rtype='radius',com='all'):
    if dr == 'net':
        return get_m(ss,R,'in',Rtype,com)-get_m(ss,R,'out',Rtype,com)
    else:
        if size(ss['d']['r'])>0:
            if Rtype=='ratio':
                Rvir=ss[com]['Rvir']
                R=R*Rvir
            mval=ss['flowdata'][com][dr]['m'][abs(ss[com]['r'] - R) == min(abs(ss[com]['r'] - R))][0]
            return mval
        else:
            return nan

def get_f(gl,ss, R, dr='net',Rtype='radius',com='all'):
    if dr == 'net':
        return get_f(gl,ss,R,'in',Rtype,com)-get_f(gl,ss,R,'out',Rtype,com)
    else:
        prev_ss = gl[gl.index(ss)-1]
        if (size(prev_ss['d']['r'])>1) & (size(ss['d']['r'])>1):
            if Rtype=='ratio':
                Rvir=ss[com]['Rvir']
                R=R*Rvir
            M = prev_ss['d']['M'][abs(prev_ss['d']['r'] - R) == min(abs(prev_ss['d']['r'] - R))][0]
            return ss['flowdata'][com][dr]['m'][abs(ss[com]['r'] - R) == min(abs(ss[com]['r'] - R))][0]/M   
        else:
            return nan

###############################################################################

def do_beep():
    # BEEP WHEN DONE
    sound_file = '/cs/sci/freundlich/CUSPCORE/beep-01a.mp3'
    display(Audio(url='//www.soundjay.com/button/beep-09.mp3', autoplay=True))  
    
###############################################################################    

def count_successes_all(sims,output,merger_threshold=[],fmean_threshold=[],fit_threshold=[],t_min=[],print_line=True,table_line=False): # WHEN INCLUDING ALL GALAXIES
    [percent_tot,t1_array_all,t2_array_all,Dsnapshot_all,delta_all,delta_i_all,fmean_all,fstd_all,fmerger_all,fprofile_all, fprofile_simname,criteria_threshold]=output
    tmin,Dsnapshot_threshold,merger_thr,Dfit_threshold,fmean_min = criteria_threshold
    ngalaxies = shape(t1_array_all)[0]
    
    if merger_threshold==[]: merger_threshold=merger_thr
    if fmean_threshold==[]:  fmean_threshold=fmean_min
    if fit_threshold==[]: fit_threshold=Dfit_threshold
    if t_min==[]: t_min=tmin
    
    count_output=[]
    n_all_tot=0
    n_nomergers_tot=0
    n_fmean_tot=0
    n_success_all_tot=0
    n_success_nomergers_tot=0
    n_success_fmean_tot=0
    
    for igal in range(ngalaxies):
        galaxy=sims[igal]
        
        fmerger = fmerger_all[igal].copy()
        fmean=fmean_all[igal].copy()
        previous_fmerger=array([nan]+fmerger[:-1].tolist())
        delta=delta_all[igal].copy()
        
        subsample_nomergers=where((t1_array_all[igal]>=t_min) & (abs(fmerger) < merger_threshold)& (abs(previous_fmerger) < merger_threshold))
        subsample_fmean=where((t1_array_all[igal]>=t_min) & (abs(fmerger) < merger_threshold)& (abs(previous_fmerger) < merger_threshold) & (fmean>=fmean_threshold))
        
        n_all=size(fmerger[where(t1_array_all[igal]>=t_min)])
        n_nomergers=size(subsample_nomergers)
        n_fmean=size(subsample_fmean)
        n_others=n_nomergers-n_fmean
        
        # COUNT SUCCESSES TOTAL
        success_cond1 = where((t1_array_all[igal]>=t_min) & (delta<=fit_threshold))
        success_cond2 = where((t1_array_all[igal]>=t_min) & (delta<=delta_i_all[igal]) | (Dsnapshot_all[igal]<=0.03 ))
        n_success_all = size(intersect1d(success_cond1,success_cond2))
        try:
            percent_all = float(n_success_all)/float(n_all)
        except:
            percent_all=nan
        
        # COUNT SUCCESSES NO MERGERS
        success_cond1 = where((delta[subsample_nomergers]<=fit_threshold))
        success_cond2 = where((delta[subsample_nomergers]<=delta_i_all[igal][subsample_nomergers]) | (Dsnapshot_all[igal][subsample_nomergers]<=0.03 ))
        n_success_nomergers = size(intersect1d(success_cond1,success_cond2))
        try:
            percent_nomergers = float(n_success_nomergers)/float(n_nomergers)
        except:
            percent_nomergers=nan
        
        # COUNT SUCCESSES NO MERGERS FMEAN
        success_cond1 = where((delta[subsample_fmean]<=fit_threshold))
        success_cond2 = where((delta[subsample_fmean]<=delta_i_all[igal][subsample_fmean]) | (Dsnapshot_all[igal][subsample_fmean]<=0.03 ))
        n_success_fmean = size(intersect1d(success_cond1,success_cond2))
        try:
            percent_fmean = float(n_success_fmean)/float(n_fmean)
        except:
            percent_fmean=nan
        
        n_success_others=n_success_nomergers-n_success_fmean
        try:
            percent_others= float(n_success_others)/float(n_others)
        except:
            percent_others=nan

        if print_line:
            print '{:8s}: '.format(galaxy)+r'%i/%i = %.0f percent success'%(n_success_all, n_all,percent_all*100.)
            print '{:8s}: '.format('')    +r'%i/%i = %.0f percent success without mergers'%(n_success_nomergers, n_nomergers,percent_nomergers*100)
            print '{:8s}: '.format('')    +r'%i/%i = %.0f percent success when fmean > %.2f'%(n_success_fmean, n_fmean,percent_fmean*100,fmean_threshold)
            print '{:8s}: '.format('')    +r'%i/%i = %.0f percent success when fmean < %.2f'%(n_success_nomergers-n_success_fmean, n_nomergers-n_fmean,percent_others*100,fmean_threshold)
        
        if table_line:
            print '{:8s}: '.format(galaxy) +r'& %i/%i = %.0f percent'%(n_success_all, n_all,percent_all*100.)+r'& %i/%i = %.0f percent'%(n_success_nomergers, n_nomergers,percent_nomergers*100)+r'& %i/%i = %.0f percent'%(n_success_fmean, n_fmean,percent_fmean*100)+r'& %i/%i = %.0f percent'%(n_success_nomergers-n_success_fmean, n_nomergers-n_fmean,percent_others*100)
        
        n_all_tot=n_all_tot+n_all
        n_nomergers_tot=n_nomergers_tot+n_nomergers
        n_fmean_tot=n_fmean_tot+n_fmean
        
        n_success_all_tot=n_success_all_tot+n_success_all
        n_success_nomergers_tot=n_success_nomergers_tot+n_success_nomergers
        n_success_fmean_tot=n_success_fmean_tot+n_success_fmean   
        
        count_output.append([galaxy,n_all,n_success_all,percent_all, n_nomergers,n_success_nomergers,percent_nomergers,n_fmean,n_success_fmean,percent_fmean,n_others,n_success_others,percent_others])

            
    n_others_tot=n_nomergers_tot-n_fmean_tot
    n_success_others_tot=n_success_nomergers_tot-n_success_fmean_tot
    
    percent_all_tot = float(n_success_all_tot)/float(n_all_tot)
    percent_nomergers_tot = float(n_success_nomergers_tot)/float(n_nomergers_tot)
    percent_fmean_tot = float(n_success_fmean_tot)/float(n_fmean_tot)
    percent_others_tot=float(n_success_others_tot)/float(n_others_tot)
    
    if print_line:
        print ' '
        print '{:8s}: '.format('TOTAL   ')+r'%i/%i = %.0f percent success'%(n_success_all_tot, n_all_tot,percent_all_tot*100.)
        print '{:8s}: '.format('')    +r'%i/%i = %.0f percent success without mergers'%(n_success_nomergers_tot, n_nomergers_tot,percent_nomergers_tot*100)
        print '{:8s}: '.format('')    +r'%i/%i = %.0f percent success when fmean > %.2f'%(n_success_fmean_tot, n_fmean_tot,percent_fmean_tot*100,fmean_threshold)
        print '{:8s}: '.format('')    +r'%i/%i = %.0f percent success when fmean < %.2f'%(n_success_others_tot, n_others_tot,percent_others_tot*100,fmean_threshold)     
    
    if table_line:
        print ' '
        print '{:8s}: '.format('TOTAL   ')+r'& %i/%i = %.0f percent'%(n_success_all_tot, n_all_tot,percent_all_tot*100.)+r'& %i/%i = %.0f percent'%(n_success_nomergers_tot, n_nomergers_tot,percent_nomergers_tot*100)+r'& %i/%i = %.0f percent'%(n_success_fmean_tot, n_fmean_tot,percent_fmean_tot*100)+r'& %i/%i = %.0f percent'%(n_success_others_tot, n_others_tot,percent_others_tot*100)  
    
    count_output.append(['TOTAL',n_all_tot,n_success_all_tot,percent_all_tot, n_nomergers_tot,n_success_nomergers_tot,percent_nomergers_tot, n_fmean_tot,n_success_fmean_tot,percent_fmean_tot, n_others_tot,n_success_others_tot,percent_others_tot])
    return count_output
         
    
def count_successes(sims,output):
    [percent_tot,t1_array_all,t2_array_all,Dsnapshot_all,delta_all,delta_i_all,fmean_all,fstd_all,fmerger_all,fprofile_all, fprofile_simname,criteria_threshold]=output
    t_min,Dsnapshot_threshold,merger_thr,Dfit_threshold,fmean_min = criteria_threshold
    ngalaxies = shape(t1_array_all)[0]
    output=[]
    n_all_tot=0
    n_fmean_tot=0
    n_success_all_tot=0
    n_success_fmean_tot=0
    for igal in range(ngalaxies):
        galaxy=sims[igal]
        subsample=where((t1_array_all[igal]>=t_min) & (abs(fmerger_all[igal]) <= merger_thr))
        fmean_subsample=where((t1_array_all[igal]>=t_min) & (abs(fmerger_all[igal]) <= merger_thr) & (fmean_all[igal]>=fmean_min))
        n_all = size(subsample)
        n_fmean = size(fmean_subsample)
        
        # Count successes
        success_cond1 = where((delta_all[igal][subsample]<=Dfit_threshold))
        success_cond2 = where((delta_all[igal][subsample]<=delta_i_all[igal][subsample]) | (Dsnapshot_all[igal][subsample]<=0.03 ))
        n_success_all = size(intersect1d(success_cond1,success_cond2))
        success_cond1 = where((delta_all[igal][fmean_subsample]<=Dfit_threshold))
        success_cond2 = where((delta_all[igal][fmean_subsample]<=delta_i_all[igal][fmean_subsample]) | (Dsnapshot_all[igal][fmean_subsample]<=0.03 ))
        n_success_fmean = size(intersect1d(success_cond1,success_cond2))
        
        # Percents
        percent_all = float(n_success_all)/float(n_all)
        percent_fmean = float(n_success_fmean)/float(n_fmean)
        percent_others= float(n_success_all-n_success_fmean)/float(n_all-n_fmean)
        output.append([galaxy,n_all,n_success_all,percent_all,n_fmean,n_success_fmean,percent_fmean])
        #print output[igal]
        print '{:8s}: '.format(galaxy)+r'%i/%i = %.0f percent success'%(n_success_all, n_all,percent_all*100.)
        print '{:8s}: '.format('')    +r'%i/%i = %.0f percent success when fmean > %.2f'%(n_success_fmean, n_fmean,percent_fmean*100,fmean_min)
        print '{:8s}: '.format('')    +r'%i/%i = %.0f percent success when fmean < %.2f'%(n_success_all-n_success_fmean, n_all-n_fmean,percent_others*100,fmean_min)
        
        n_all_tot=n_all_tot+n_all
        n_fmean_tot=n_fmean_tot+n_fmean
        n_success_all_tot=n_success_all_tot+n_success_all
        n_success_fmean_tot=n_success_fmean_tot+n_success_fmean       
    
    percent_all_tot = float(n_success_all_tot)/float(n_all_tot)
    percent_fmean_tot = float(n_success_fmean_tot)/float(n_fmean_tot)
    percent_others_tot=float(n_success_all_tot-n_success_fmean_tot)/float(n_all_tot-n_fmean_tot)
    
    print ' '
    print '{:8s}: '.format('TOTAL   ')+r'%i/%i = %.0f percent success'%(n_success_all_tot, n_all_tot,percent_all_tot*100.)
    print '{:8s}: '.format('')    +r'%i/%i = %.0f percent success when fmean > %.2f'%(n_success_fmean_tot, n_fmean_tot,percent_fmean_tot*100,fmean_min)
    print '{:8s}: '.format('')    +r'%i/%i = %.0f percent success when fmean < %.2f'%(n_success_all_tot-n_success_fmean_tot, n_all_tot-n_fmean_tot,percent_others_tot*100,fmean_min)
            output.append(['TOTAL',n_all_tot,n_success_all_tot,percent_all_tot,n_fmean_tot,n_success_fmean_tot,percent_fmean_tot])
    return output      
        
def plot_successes(output):
    [percent_tot,t1_array_all,t2_array_all,Dsnapshot_all,delta_all,delta_i_all,fmean_all,fstd_all,fmerger_all,fprofile_all, fprofile_simname,criteria_threshold] =output
    t_min,Dsnapshot_threshold,merger_thr,Dfit_threshold,fmean_min = criteria_threshold
    ngalaxies = shape(t1_array_all)[0]
    output=[]
    for igal in range(ngalaxies):
        subsample=where((t1_array_all[igal]>=t_min) & (abs(fmerger_all[igal]) <= merger_thr))
        fmean_subsample=where((t1_array_all[igal]>=t_min) & (abs(fmerger_all[igal]) <= merger_thr) & (fmean_all[igal]>=fmean_min))
        success_cond1 = where((delta_all[igal][subsample]<=Dfit_threshold))
        success_cond2 = where((delta_all[igal][subsample]<=delta_i_all[igal][subsample]) | (Dsnapshot_all[igal][subsample]<=0.03 ))
        successful = intersect1d(success_cond1,success_cond2)
        #
        # FIGURE 1
        figure()
        scatter(fmean_all[igal][subsample],delta_all[igal][subsample],color='b',label='fail')
        scatter(fmean_all[igal][subsample][successful],delta_all[igal][subsample][successful],color='r',label='success')
        axhline(y=Dfit_threshold,color='k')
        axvline(x=fmean_min,color='k')
        xlabel(r'$|f| mean$')
        ylabel(r'$\Delta fit$')
        legend(loc='upper left')
        #
        # FIGURE 2
        figure()
        scatter(fstd_all[igal][subsample],delta_all[igal][subsample],color='b',label='fail')
        scatter(fstd_all[igal][subsample][successful],delta_all[igal][subsample][successful],color='r',label='success')
        axhline(y=Dfit_threshold,color='k')
        axvline(x=fmean_min,color='k')
        xlabel(r'$|f| std$')
        ylabel(r'$\Delta fit$')    
        legend(loc='upper left')
        #
        # FIGURE 3
        figure()
        scatter(fmean_all[igal][subsample],Dsnapshot_all[igal][subsample],color='b',label='fail')
        scatter(fmean_all[igal][subsample][successful],Dsnapshot_all[igal][subsample][successful],color='r',label='success')
        axhline(y=Dfit_threshold,color='k')
        axvline(x=fmean_min,color='k')
        xlabel(r'$|f| mean$')
        ylabel(r'$\Delta snap$')    
        legend(loc='upper left')  
        
def plot_successes_all(output):
    [percent_tot,t1_array_all,t2_array_all,Dsnapshot_all,delta_all,delta_i_all,fmean_all,fstd_all,fmerger_all,fprofile_all, fprofile_simname,criteria_threshold] =output
    t_min,Dsnapshot_threshold,merger_thr,Dfit_threshold,fmean_min = criteria_threshold
    ngalaxies = shape(t1_array_all)[0]
    #
    # FIGURE 1
    figure()
    axhline(y=Dfit_threshold,color='k')
    axvline(x=fmean_min,color='k')
    xlabel(r'$|f| mean$')
    ylabel(r'$\Delta fit$')
    for igal in range(ngalaxies):
        subsample=where((t1_array_all[igal]>=t_min) & (abs(fmerger_all[igal]) <= merger_thr))
        fmean_subsample=where((t1_array_all[igal]>=t_min) & (abs(fmerger_all[igal]) <= merger_thr) & (fmean_all[igal]>=fmean_min))
        success_cond1 = where((delta_all[igal][subsample]<=Dfit_threshold))
        success_cond2 = where((delta_all[igal][subsample]<=delta_i_all[igal][subsample]) | (Dsnapshot_all[igal][subsample]<=0.03 ))
        successful = intersect1d(success_cond1,success_cond2)
        #
        scatter(fmean_all[igal][subsample],delta_all[igal][subsample],color='b',label='fail')
        scatter(fmean_all[igal][subsample][successful],delta_all[igal][subsample][successful],color='r',label='success')
    #
    # FIGURE 2
    figure()
    axhline(y=Dfit_threshold,color='k')
    xlabel(r'$|f| std$')
    ylabel(r'$\Delta fit$')    
    for igal in range(ngalaxies):
        subsample=where((t1_array_all[igal]>=t_min) & (abs(fmerger_all[igal]) <= merger_thr))
        fmean_subsample=where((t1_array_all[igal]>=t_min) & (abs(fmerger_all[igal]) <= merger_thr) & (fmean_all[igal]>=fmean_min))
        success_cond1 = where((delta_all[igal][subsample]<=Dfit_threshold))
        success_cond2 = where((delta_all[igal][subsample]<=delta_i_all[igal][subsample]) | (Dsnapshot_all[igal][subsample]<=0.03 ))
        successful = intersect1d(success_cond1,success_cond2)
        #
        scatter(fstd_all[igal][subsample],delta_all[igal][subsample],color='b',label='fail')
        scatter(fstd_all[igal][subsample][successful],delta_all[igal][subsample][successful],color='r',label='success')
    #
    # FIGURE 3
    figure()
    axhline(y=Dsnapshot_threshold,color='k')
    axvline(x=fmean_min,color='k')
    xlabel(r'$|f| mean$')
    ylabel(r'$\Delta snap$')    
    for igal in range(ngalaxies):
        subsample=where((t1_array_all[igal]>=t_min) & (abs(fmerger_all[igal]) <= merger_thr))
        fmean_subsample=where((t1_array_all[igal]>=t_min) & (abs(fmerger_all[igal]) <= merger_thr) & (fmean_all[igal]>=fmean_min))
        success_cond1 = where((delta_all[igal][subsample]<=Dfit_threshold))
        success_cond2 = where((delta_all[igal][subsample]<=delta_i_all[igal][subsample]) | (Dsnapshot_all[igal][subsample]<=0.03 ))
        successful = intersect1d(success_cond1,success_cond2)
        #
        scatter(fmean_all[igal][subsample],Dsnapshot_all[igal][subsample],color='b',label='fail')
        scatter(fmean_all[igal][subsample][successful],Dsnapshot_all[igal][subsample][successful],color='r',label='success')
        
def latex_float(f):
    float_str = "{0:.3g}".format(f)
    if "e" in float_str:
        base, exponent = float_str.split("e")
        return r"{0} \times 10^{{{1}}}".format(base, int(exponent))
    else:
        return float_str
 
def success_table(sims,output):
    [percent_tot,t1_array_all,t2_array_all,Dsnapshot_all,delta_all,delta_i_all,fmean_all,fstd_all,fmerger_all,fprofile_all, fprofile_simname,criteria_threshold]=output
    t_min,Dsnapshot_threshold,merger_thr,Dfit_threshold,fmean_min = criteria_threshold
    ngalaxies = shape(t1_array_all)[0]
    output=[]
    n_all_tot=0
    n_fmean_tot=0
    n_success_all_tot=0
    n_success_fmean_tot=0
    for igal in range(ngalaxies):
        galaxy=sims[igal]
        subsample=where((t1_array_all[igal]>=t_min) & (abs(fmerger_all[igal]) <= merger_thr))
        fmean_subsample=where((t1_array_all[igal]>=t_min) & (abs(fmerger_all[igal]) <= merger_thr) & (fmean_all[igal]>=fmean_min))
        n_all = size(subsample)
        n_fmean = size(fmean_subsample)
        
        # Mstar, Mvir
        ok_fangzhou,r12_fangzhou,rvir_fangzhou,mvir_fangzhou,mstar_fangzhou=get_fangzhou_radii(galaxy,[1.],get_stars=True)
        
        # Count successes
        success_cond1 = where((delta_all[igal][subsample]<=Dfit_threshold))
        success_cond2 = where((delta_all[igal][subsample]<=delta_i_all[igal][subsample]) | (Dsnapshot_all[igal][subsample]<=0.03 ))
        n_success_all = size(intersect1d(success_cond1,success_cond2))
        success_cond1 = where((delta_all[igal][fmean_subsample]<=Dfit_threshold))
        success_cond2 = where((delta_all[igal][fmean_subsample]<=delta_i_all[igal][fmean_subsample]) | (Dsnapshot_all[igal][fmean_subsample]<=0.03 ))
        n_success_fmean = size(intersect1d(success_cond1,success_cond2))
        
        # Percents
        percent_all = float(n_success_all)/float(n_all)
        percent_fmean = float(n_success_fmean)/float(n_fmean)
        percent_others= float(n_success_all-n_success_fmean)/float(n_all-n_fmean)
        output.append([galaxy,n_all,n_success_all,percent_all,n_fmean,n_success_fmean,percent_fmean])
        print "\ texttt{%8s}"%galaxy+ "& $%s$"%latex_float(mvir_fangzhou[0])+ "& $%s$"%latex_float(mstar_fangzhou[0])+" & $%.0f \percent$ $(%i/%i)$"%(percent_fmean*100,n_success_fmean, n_fmean)+" & $%.0f \percent$ $(%i/%i)$"%(percent_others*100,n_success_all-n_success_fmean, n_all-n_fmean)+" & $%.0f \percent$ $(%i/%i)$ \\"%(percent_all*100.,n_success_all, n_all)
        
        n_all_tot=n_all_tot+n_all
        n_fmean_tot=n_fmean_tot+n_fmean
        n_success_all_tot=n_success_all_tot+n_success_all
        n_success_fmean_tot=n_success_fmean_tot+n_success_fmean       
    
    percent_all_tot = float(n_success_all_tot)/float(n_all_tot)
    percent_fmean_tot = float(n_success_fmean_tot)/float(n_fmean_tot)
    percent_others_tot=float(n_success_all_tot-n_success_fmean_tot)/float(n_all_tot-n_fmean_tot)
    
    print ' '
    print 'All'+'& & & $%.0f \percent$ $(%i/%i)$'%(percent_fmean_tot*100,n_success_fmean_tot, n_fmean_tot)+' & $%.0f \percent$ $(%i/%i)$'%(percent_others_tot*100,n_success_all_tot-n_success_fmean_tot, n_all_tot-n_fmean_tot)+' & $%.0f \percent$ $(%i/%i)$ \\'%(percent_all_tot*100.,n_success_all_tot, n_all_tot)
        output.append(['TOTAL',n_all_tot,n_success_all_tot,percent_all_tot,n_fmean_tot,n_success_fmean_tot,percent_fmean_tot])
    return output #print 'mean success rate: ', mean(output[:][])        

def successes_table_all(sims,output,merger_threshold=[],fmean_threshold=[],fit_threshold=[],t_min=[],print_line=True,table_line=False,use_previous_mergers=0,use_cond2=True): # WHEN INCLUDING ALL GALAXIES
    [percent_tot,t1_array_all,t2_array_all,Dsnapshot_all,delta_all,delta_i_all,fmean_all,fstd_all,fmerger_all,fprofile_all, fprofile_simname,criteria_threshold]=output
    tmin,Dsnapshot_threshold,merger_thr,Dfit_threshold,fmean_min = criteria_threshold
    ngalaxies = shape(t1_array_all)[0]
    print 'Number of galaxies: ', ngalaxies
    
    if merger_threshold==[]: merger_threshold=merger_thr
    if fmean_threshold==[]:  fmean_threshold=fmean_min
    if fit_threshold==[]: fit_threshold=Dfit_threshold
    if t_min==[]: t_min=tmin
    if use_previous_mergers<>0:
        previous_merger_threshold=merger_threshold
    else:
        previous_merger_threshold=10.
    
    # FANGZHOU'S CATALOG
    data=genfromtxt('/cs/sci/freundlich/CUSPCORE/catalogs/NIHAO_a1.0000.txt',skip_header =1)
    
    count_output=[]
    n_all_tot=0
    n_nomergers_tot=0
    n_fmean_tot=0
    n_success_all_tot=0
    n_success_nomergers_tot=0
    n_success_fmean_tot=0
    
    for igal in range(ngalaxies):
        galaxy=sims[igal]
        
        fmerger = fmerger_all[igal].copy()
        fmean=fmean_all[igal].copy()
        delta=delta_all[igal].copy()
        
        #print igal
        previous_fmerger=array([nan]+fmerger[:-1].copy().tolist())
        if use_previous_mergers<>0:
            for i in range(use_previous_mergers):
                previous_fmergeri=array([nan]*(i+1)+fmerger[:-(i+1)].copy().tolist())
                previous_fmerger=maximum(abs(previous_fmerger),abs(previous_fmergeri))
    
        # Mstar, Mvir
        ok_fangzhou,r12_fangzhou,rvir_fangzhou,mvir_fangzhou,mstar_fangzhou=get_fangzhou_radii(galaxy,[1.],get_stars=True)
        
        subsample_nomergers=where((t1_array_all[igal]>=t_min) & (abs(fmerger) < merger_threshold)& (abs(previous_fmerger) < previous_merger_threshold) & (~isnan(delta)))
        subsample_fmean=where((t1_array_all[igal]>=t_min) & (abs(fmerger) < merger_threshold)& (abs(previous_fmerger) < previous_merger_threshold) & (fmean>=fmean_threshold) & (~isnan(delta)))
        
        n_all=size(fmerger[where((t1_array_all[igal]>=t_min) & (~isnan(delta)))])
        n_nomergers=size(subsample_nomergers)
        n_fmean=size(subsample_fmean)
        n_others=n_nomergers-n_fmean
        #print n_all, n_nomergers, n_fmean, n_others
        
        # COUNT SUCCESSES TOTAL
        success_cond1 = where((t1_array_all[igal]>=t_min) & (delta<=fit_threshold) & (~isnan(delta)))
        success_cond2 = where((t1_array_all[igal]>=t_min) & (~isnan(delta)) & (delta<=delta_i_all[igal]) | (Dsnapshot_all[igal]<=0.03 ))
        #print size(success_cond1), size(success_cond2)
        if use_cond2:
            n_success_all = size(intersect1d(success_cond1,success_cond2))
        else:
            n_success_all = size(success_cond1)
        try:
            percent_all = float(n_success_all)/float(n_all)
        except:
            percent_all=nan
        
        # COUNT SUCCESSES NO MERGERS
        success_cond1 = where((delta[subsample_nomergers]<=fit_threshold))
        success_cond2 = where((delta[subsample_nomergers]<=delta_i_all[igal][subsample_nomergers]) | (Dsnapshot_all[igal][subsample_nomergers]<=0.03 ))
        if use_cond2:
            n_success_nomergers = size(intersect1d(success_cond1,success_cond2))
        else:
            n_success_nomergers = size(success_cond1)
        try:
            percent_nomergers = float(n_success_nomergers)/float(n_nomergers)
        except:
            percent_nomergers=nan
        
        # COUNT SUCCESSES NO MERGERS FMEAN
        success_cond1 = where((delta[subsample_fmean]<=fit_threshold))
        success_cond2 = where((delta[subsample_fmean]<=delta_i_all[igal][subsample_fmean]) | (Dsnapshot_all[igal][subsample_fmean]<=0.03 ))
        if use_cond2:
            n_success_fmean = size(intersect1d(success_cond1,success_cond2))
        else:
            n_success_fmean = size(success_cond1)
        try:
            percent_fmean = float(n_success_fmean)/float(n_fmean)
        except:
            percent_fmean=nan
        
        n_success_others=n_success_nomergers-n_success_fmean
        try:
            percent_others= float(n_success_others)/float(n_others)
        except:
            percent_others=nan

        # GET RE, UDG
        dagger=''
        try:
            I=float(galaxy.replace('g',''))
            i_ID=where(data[:,0]==I)[0][0]
            mstar=data[i_ID,17]
            msub=data[i_ID,7]
            re_star=data[i_ID,21]
            muV_central=data[i_ID,28]
            UDG=False
            if muV_central>=24 and re_star>=1.5:
                UDG=True
            if UDG==True:
                dagger='$^\dagger$'
        except:
            msub=nan
            mstar=nan
            re_star=nan
            muV_central=nan
            UDG=False
            print 'Galaxy %s not in Fangzhou catalog'%galaxy
        
        
        
        if galaxy=='g1.08e11':
           print "textbf{texttt{%8s}}"%galaxy+dagger+ "& $mathbf{%s}$"%latex_float(mvir_fangzhou[0])+ "& $%.1f$"%rvir_fangzhou[0]+ "& $mathbf{%s}$"%latex_float(mstar)+'& $mathbf{%.2f}$ '%re_star+'& $mathbf{%.2f}$'%muV_central+'& textbf{%02i/%02i = %.0f\percent}'%(n_success_all, n_all,percent_all*100.)+'& textbf{%02i/%02i = %.0f\percent}'%(n_success_nomergers, n_nomergers,percent_nomergers*100)+'& textbf{%02i/%02i = %.0f\percent}'%(n_success_fmean, n_fmean,percent_fmean*100)+'& textbf{%02i/%02i = %.0f\percent}'%(n_success_nomergers-n_success_fmean, n_nomergers-n_fmean,percent_others*100)+' endline' 
        else:
            print "texttt{%8s}"%galaxy+dagger+ "& $%s$"%latex_float(mvir_fangzhou[0])+ "& $%.1f$"%rvir_fangzhou[0]+ "& $%s$"%latex_float(mstar)+'& $%.2f$ '%re_star+'& $%.2f$'%muV_central+'& %02i/%02i = %.0f\percent'%(n_success_all, n_all,percent_all*100.)+'& %02i/%02i = %.0f\percent'%(n_success_nomergers, n_nomergers,percent_nomergers*100)+'& %02i/%02i = %.0f\percent'%(n_success_fmean, n_fmean,percent_fmean*100)+'& %02i/%02i = %.0f\percent'%(n_success_nomergers-n_success_fmean, n_nomergers-n_fmean,percent_others*100)+' endline' 
        
        n_all_tot=n_all_tot+n_all
        n_nomergers_tot=n_nomergers_tot+n_nomergers
        n_fmean_tot=n_fmean_tot+n_fmean
        
        n_success_all_tot=n_success_all_tot+n_success_all
        n_success_nomergers_tot=n_success_nomergers_tot+n_success_nomergers
        n_success_fmean_tot=n_success_fmean_tot+n_success_fmean   
        
        count_output.append([galaxy,n_all,n_success_all,percent_all, n_nomergers,n_success_nomergers,percent_nomergers,n_fmean,n_success_fmean,percent_fmean,n_others,n_success_others,percent_others])

            
    n_others_tot=n_nomergers_tot-n_fmean_tot
    n_success_others_tot=n_success_nomergers_tot-n_success_fmean_tot
    
    percent_all_tot = float(n_success_all_tot)/float(n_all_tot)
    percent_nomergers_tot = float(n_success_nomergers_tot)/float(n_nomergers_tot)
    percent_fmean_tot = float(n_success_fmean_tot)/float(n_fmean_tot)
    percent_others_tot=float(n_success_others_tot)/float(n_others_tot)
    
    print 'textbf{All}'+'& & & & & '+r'& textbf{%i/%i = %.0f\percent}'%(n_success_all_tot, n_all_tot,percent_all_tot*100.)+r'& textbf{%i/%i = %.0f\percent}'%(n_success_nomergers_tot, n_nomergers_tot,percent_nomergers_tot*100)+r'& textbf{%i/%i = %.0f\percent}'%(n_success_fmean_tot, n_fmean_tot,percent_fmean_tot*100)+r'& textbf{%i/%i = %.0f\percent}'%(n_success_others_tot, n_others_tot,percent_others_tot*100)+' endline'
    return [percent_all_tot*100.,percent_nomergers_tot*100,percent_fmean_tot*100,percent_others_tot*100]
    
        