global sim, directory,components,k
global rmin_evolve,rmax_evolve,rmin_fit,rmax_fit,linear_slopes,Ttype,Ptype
global Dsnapshot_threshold,Dfit_threshold,constrain_fit,constrain_evolution
global xi,xi_merger,merger_thr,f_min,fmean_min,t_min,m_constraint,com,delta_xlim,delta_ymin
global limits,limits_fprofile,limits_T,limits_abg,xticks_val,xticks_fprofile
global polyorder,sigma,mode,double_smooth,linear_slopes
global textfont
global selection,multiple_snapshots,plot_criterion
global figsize,nrows,n

rmin_evolve=0.01
rmax_evolve=1.
rmin_fit=0.01
rmax_fit=1.
Ttype='jeans'
Ptype=None
Dsnapshot_threshold=0.05
Dfit_threshold=0.07
constrain_fit=False
constrain_evolution=False

xi=0.01
xi_merger=0.1
merger_thr=0.15
f_min=0.10
fmean_min=0.06
t_min=2.
m_constraint=2.
com='all'
delta_xlim=[-2,-1.5]
delta_ymin=2.

limits=[-2,-1,2,3]
limits_fprofile=[-2,0.,0,0.5]
limits_T=(-2,0,0,4)
limits_abg=(-2,0,-1,1)
xticks_val=(-2,-1.5,-1)
xticks_fprofile=(-2,-1.5,-1,-0.5,0)

polyorder=3
sigma = 21
mode= 'interp'
double_smooth=False
linear_slopes=False
#rlim=[-2,0]

textfont=12

directory='/cs/sci/freundlich/CUSPCORE/Michael/'
sims = ['g1.08e11']
sim = sims[0]
components=['d']
k=32
selection=[]

multiple_snapshots=False
plot_criterion='f'
figsize=(20,20)
nrows=8
ncols=8

D200=False