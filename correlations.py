import scipy
import numpy
import sys
sys.path.insert(0, '/cs/sci/freundlich/CUSPCORE/Analysis/Model')
sys.path.insert(0, '/cs/sci/freundlich/CUSPCORE/Analysis/General')
import cuspcore_go
reload(cuspcore_go)
from cuspcore_aux import *
import general_functions
from general_functions import *
reload(general_functions)
from scipy.special import erf
import matplotlib.patches as patches

def m(x, w):
    """Weighted Mean"""
    return np.sum(x * w) / np.sum(w)

def cov(x, y, w):
    """Weighted Covariance"""
    return np.sum(w * (x - m(x, w)) * (y - m(y, w))) / np.sum(w)

def corr(x, y, w):
    """Weighted Correlation"""
    return cov(x, y, w) / np.sqrt(cov(x, x, w) * cov(y, y, w))

def plot_correlation(x,y,qname='q',w=[],success=[],c=[],figsize=(7,5),xname=r'x',yname=r'y',fontsize=24,Dfit_threshold=0.07,xlim=[],ylim=[],clabel=[],savefigure=False,x_ticks=[],x_labels=[],y_ticks=[],y_labels=[],plot_line=True,standalone=True,ax=[],return_p=False,do_annotate=False,plot_errorbar=True,plot_colorbar=False,cname='t',c_ticks=[-2.5,-2.0,-1.5,-1.,-0.5,0],c_labels=[r'$-2.5$',r'$-2.0$',r'$-1.5$',r'$-1.0$',r'$-0.5$',r'$-0.0$'],rasterized= False):   
    
    #print qname
    space_down=0.18
    space_left=0.07
    width=0.21
    height=width*18./6.*3.5/3.
    
    x=array_nonan(x,y)
    if success<>[]: s=array_nonan(success,y)==1
    if c<>[]: c=array_nonan(c,y)
    if w<>[]: w=array_nonan(w,y)
    y=array_nonan(y,y)
    
    y=array_nonan(y,x)
    if success<>[]: s=array_nonan(s,x)==1
    if c<>[]: c=array_nonan(c,x)
    if w<>[]: w=array_nonan(w,x)
    x=array_nonan(x,x)

    if w==[]: w=ones_like(x)
    
    if c<>[]:
        cmin=c.min()
        cmax=c.max()
    else:
        c='k'
        cmin=nan
        cmax=nan
       
    p = numpy.polyfit(x, y, 1)
    chi_squared = numpy.sum((numpy.polyval(p, x) - y) ** 2)
    
    rpearson=scipy.stats.pearsonr(x,y)[0]
    stdev=abs(std(y-p[1]-p[0]*x))
    rcor=corr(x,y,w)
    
    if success<>[]:
        p_success = numpy.polyfit(x[s], y[s], 1)
        p_failure = numpy.polyfit(x[invert(s)], y[invert(s)], 1)
        r_success = scipy.stats.pearsonr(x[s], y[s])[0]
        r_failure = scipy.stats.pearsonr(x[invert(s)], y[invert(s)])[0]
        print "p_success = ", p_success
        print "p_failure = ", p_failure
        print "r_success = ", r_success
        print "r_failure = ", r_failure

    # Bin averages
    bins = numpy.linspace(x.min(), x.max(), 6)
    digitized_x = numpy.digitize(x, bins)
    x_means = [x[digitized_x == i].mean() for i in range(1, len(bins))]
    y_means = [y[digitized_x == i].mean() for i in range(1, len(bins))]
    x_std   = [x[digitized_x == i].std() for i in range(1, len(bins))]
    y_std   = [y[digitized_x == i].std() for i in range(1, len(bins))]
    
    y_nsuccess = asarray([float(size(where(y[digitized_x == i]<Dfit_threshold))) for i in range(1, len(bins))])
    y_ntot     = asarray([float(size(y[digitized_x == i])) for i in range(1, len(bins))])
    psuccess   = y_nsuccess/y_ntot
    dpsuccess  = psuccess*sqrt(1./y_nsuccess+1./y_ntot) # ASSUMING POISSON NOISE
    dpsuccess[where(psuccess==0)]=0
        
    # SCATTER PLOT 
    if standalone:
        fig, ax = subplots(nrows=1, ncols=1,figsize=figsize)
        clf()
        ax=gca()
    else:
        if ax==[]:
            ax=gca()
            
    if plot_line: 
        ax.axhline(y=Dfit_threshold,color='gray')
    
    if c<>[]: 
        #scatter_plot =
        mappable= ax.scatter(x,y,c=c,marker='o',s=20,vmin=cmin,vmax=cmax,edgecolor='None',alpha=0.7,rasterized= rasterized)#,linewidth=0.5)#,facecolors='r')#c=z
    else:
        mappable= ax.scatter(x,y,c='r',marker='o',s=20,edgecolor=None,rasterized= rasterized)#,alpha=0.5)

    #ax=gca()
    xmin, xmax = ax.get_xlim()
    if xlim<>[]:
        xmin=xlim[0]
        xmax=xlim[1]
    
    if standalone:
        ax.set_position([space_left+0.13,space_down, height*5./7., height])
    
    if plot_errorbar:
        ax.errorbar(x_means,y_means,xerr=x_std,yerr=y_std,ecolor='k',fmt=None,linewidth=2)

    abscisse=linspace(xmin,xmax,100)
    ax.plot(abscisse,p[1]+p[0]*abscisse,'k')
    print qname, ': y = %.2f + %.4f x '%(p[1],p[0])
    
    if success<>[]:
        ax.scatter(x[s],y[s],marker='o',s=40,c='b',rasterized= rasterized)
        ax.scatter(x[invert(s)], y[invert(s)],marker='o',s=40,c='r',rasterized= rasterized)
        ax.plot(abscisse,p_success[1]+p_success[0]*abscisse,'b')
        ax.plot(abscisse,p_failure[1]+p_failure[0]*abscisse,'r')
    
    ax.set_xlabel(xname,fontsize=fontsize)
    ax.set_ylabel(yname,fontsize=fontsize)    
    rect=patches.Rectangle((0.01,0.01),0.3,0.15,transform=ax.transAxes,facecolor='white',color='white',alpha=0.6,zorder=11)
    ax.add_patch(rect)
    ax.annotate(r'$%.2f$'%abs(rpearson), xy=(0.03, 0.03),xycoords='axes fraction',fontsize=fontsize-4,zorder=12)
    if do_annotate:
        annotate(r'$%.2f$'%abs(rpearson), xy=(0.03, 0.03),xycoords='axes fraction',fontsize=fontsize-4)

    if x_ticks<>[]:
        ax.set_xticks(x_ticks)
        if x_labels<>[]:
            ax.set_xticklabels(x_labels,fontsize=fontsize-8)
    
    if y_ticks<>[]:
        ax.set_yticks(y_ticks)
        if y_labels<>[]:
            ax.set_yticklabels(y_labels,fontsize=fontsize-8)
    
    if xlim<>[]:
        ax.set_xlim(xlim)
    else:
        ax.set_xlim([xmin,xmax])
        
    if ylim<>[]:
        ax.set_ylim(ylim)
    
    if plot_colorbar:
        fig=gcf()
        cbar_ax = fig.add_axes([space_left+0.13+0.06+height*5./7.,space_down, 0.02*18./7., height])
        cbar=colorbar(mappable,cax=cbar_ax)
        cbar.set_ticks(c_ticks)
        cbar.set_ticklabels(c_labels)
        cbar.set_label(cname,fontsize=fontsize)
        tick_params(axis='both', which='major', labelsize=fontsize-4,direction='in')
            
    if standalone and savefigure:
        savefig('/cs/sci/freundlich/CUSPCORE/ARTICLE/images/correlations/fig_success_%s.pdf'%qname)
    
    if return_p:
        return p
    if not standalone:
        return ax
    else: 
        show()

def plot_success(x,y,qname='q',success=[],c=[],figsize=(7,5),xname=r'x',yname=r'y',fontsize=24,Dfit_threshold=0.07,xlim=[],ylim=[],clabel=[],savefigure=False,x_ticks=[],x_labels=[],y_ticks=[],y_labels=[],standalone=True,ax=[],rasterized= False):   

    space_down=0.18
    space_left=0.07
    width=0.21
    height=width*18./6.*3.5/3.#0.78

    x=array_nonan(x,y)
    if success<>[]: s=array_nonan(success,y)==1
    if c<>[]: c=array_nonan(c,y)
    y=array_nonan(y,y)
    
    y=array_nonan(y,x)
    if success<>[]: s=array_nonan(s,x)==1
    if c<>[]: c=array_nonan(c,x)
    x=array_nonan(x,x)

    if c<>[]:
        cmin=c.min()
        cmax=c.max()
    else:
        c='k'
        cmin=nan
        cmax=nan
    
    p = numpy.polyfit(x, y, 1)
    stdev=abs(std(y-p[1]-p[0]*x))
    
    # Bin averages
    bins = numpy.linspace(x.min(), x.max(), 6)
    digitized_x = numpy.digitize(x, bins)
    x_means = [x[digitized_x == i].mean() for i in range(1, len(bins))]
    y_means = [y[digitized_x == i].mean() for i in range(1, len(bins))]
    if c<>'k':
        c_means = [c[digitized_x == i].mean() for i in range(1, len(bins))]
    x_std   = [x[digitized_x == i].std() for i in range(1, len(bins))]
    y_std   = [y[digitized_x == i].std() for i in range(1, len(bins))]
    
    y_nsuccess = asarray([float(size(where(y[digitized_x == i]<Dfit_threshold))) for i in range(1, len(bins))])
    y_ntot     = asarray([float(size(y[digitized_x == i])) for i in range(1, len(bins))])
    psuccess   = y_nsuccess/y_ntot
    dpsuccess  = psuccess*sqrt(1./y_nsuccess+1./y_ntot) # ASSUMING POISSON NOISE
    dpsuccess[where(psuccess==0)]=sqrt(1./y_ntot[where(psuccess==0)])
    
    rpearson=scipy.stats.pearsonr(x_means,psuccess)[0]

    # SUCCESS PROBABILITY PLOT
    #if plot_proba:
    if standalone:
        fig, ax = subplots(nrows=1, ncols=1,figsize=figsize)
        clf()
        ax=gca()
    else:
        if ax==[]:
            ax=gca()
    
    if c<>'k':
        ax.scatter(x_means,psuccess,c=c_means,marker='o',s=70,vmin=cmin,vmax=cmax,zorder=10,rasterized= rasterized)
    
    ax.errorbar(x_means,psuccess,xerr=x_std,yerr=dpsuccess,ecolor='k',fmt=None,linewidth=2,zorder=0)
    ax.axhline(y=0.5,color='gray')

    xmin, xmax = ax.get_xlim()
    if xlim<>[]:
        xmin=xlim[0]
        xmax=xlim[1]
    
    if xlim==[]:
        ylim=[0,1]
    
    abscisse=linspace(xmin,xmax,100)
    ax.plot(abscisse,0.5*(1+erf((Dfit_threshold-p[1]-p[0]*abscisse)/sqrt(2.*stdev**2))),'k',lw=2)

    if x_ticks<>[]:
        ax.set_xticks(x_ticks)
        if x_labels<>[]:
            ax.set_xticklabels(x_labels,fontsize=fontsize-8)
    
    if y_ticks<>[]:
        ax.set_yticks(y_ticks)
        if y_labels<>[]:
            ax.set_xticklabels(y_labels,fontsize=fontsize-8)
    
    ax.set_xlim([xmin,xmax])
    ax.set_ylim(ylim)

    if standalone:
        ax.set_position([space_left+0.13,space_down, height*5./7., height])

    ax.set_xlabel(xname,fontsize=fontsize)
    ax.set_ylabel(r'$p$',fontsize=fontsize)
     rect=patches.Rectangle((0.01,0.01),0.3,0.15,transform=ax.transAxes,facecolor='white',color='white',alpha=0.6,zorder=11)
    ax.add_patch(rect)
    
    if standalone and savefigure:
        savefig('/cs/sci/freundlich/CUSPCORE/ARTICLE/images/correlations/fig_success_p_%s.pdf'%qname)
     
    if not standalone:
        return ax
    else:
        show()
    

def plot_success2D(x,y,delta,nbins=5,figsize=(7,5),xname=r'x',yname=r'y',fontsize=24,Dfit_threshold=0.07,xlim=[],ylim=[],xticks=[],yticks=[],clabel=[],savefigure=False): 
    space_down=0.18
    space_left=0.07
    width=0.21
    height=width*18./6.*3.5/3.#0.78

    x=array_nonan(x,y)
    delta=array_nonan(delta,y)
    y=array_nonan(y,y)
    
    y=array_nonan(y,x)
    delta=array_nonan(delta,x)
    x=array_nonan(x,x) 
    
    weights=1*(delta<Dfit_threshold)

    fig, axes = subplots(nrows=1, ncols=1,figsize=figsize)
    clf()
    
    H, xedges, yedges = np.histogram2d(x, y, bins=nbins, weights=weights)
    H2, _, _ = np.histogram2d(x,y, bins=nbins)
    extent = [0,1, xedges[0],xedges[-1]]
    plt.imshow(H/H2, extent=extent,interpolation='nearest',aspect='auto')
    
    ax=gca()
    ax.set_position([space_left+0.13,space_down, height*5./7., height])

    xlabel(yname,fontsize=fontsize)
    ylabel(xname,fontsize=fontsize)

    if xticks<>[]:
        plt.xticks(xticks)
    
    if yticks<>[]:
        plt.yticks(yticks)
    
    clabel=r'$p$'
    cbar_ax = fig.add_axes([space_left+0.13+0.06+height*5./7.,space_down, 0.02*18./7., height])
    cbar=colorbar(label=clabel,cax=cbar_ax)

    if savefigure:
        savefig('/cs/sci/freundlich/CUSPCORE/ARTICLE/images/correlations/fig_success_p2d_%s.pdf'%qname)
    
    show()