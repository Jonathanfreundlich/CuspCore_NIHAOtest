# LOAD DATA 

global fitname, filestart,rlim
global gl,treal,glsim,times
global a_array,t1_array,t2_array,rvir_halo,mvir_halo
global ok_fangzhou,r12_fangzhou,rvir_fangzhou,mvir_fangzhou,relevent_fangzhou
global fmerger,rms,Dsnapshots,delta,delta_i,delta_if,fmean,fstd
global components

#############################################################
# DERIVED PARAMETERS
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

rlim=[log10(rmin_fit),log10(rmax_fit)] 

#############################################################

print 'Simulation %s'%sim
with open(directory+'NIHAO-%s.pickle'%sim[1:]) as f:
    gl = pickle.load(f)   

gl = slopes_functions.derive_slopes(gl,polyorder=polyorder,sigma=sigma,mode=mode,double_smooth=double_smooth,rlim=rlim,use_fangzhou_Rvir=True,linearize=linear_slopes,betanull=False,D200=D200)
gl = prepare_functions.define_brho(gl,polyorder=polyorder,sigma=sigma,mode=mode,double_smooth=double_smooth,rlim=rlim,use_fangzhou_Rvir=True,D200=D200)
treal=treal_functions.load_or_create_gl(sim,use_fangzhou_Rvir=True,D200=D200)
fitrange=prepare_functions.get_fitrange(gl,use_fangzhou_Rvir=True,D200=D200)
gl=prepare_functions.reduce_range_gl(gl,fitrange)
treal=prepare_functions.reduce_range_Treal(treal,fitrange)
print ' '

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
        t1_array.append(ss['flowdata']['t1'])
        t2_array.append(ss['flowdata']['t2']) 

a_array=array(a_array)
t1_array=array(t1_array)
t2_array=array(t2_array)
rvir_halo=array(rvir_halo)
mvir_halo=array(mvir_halo)

# GET RADII FROM FANGZHOU'S FILES
ok_fangzhou,r12_fangzhou,rvir_fangzhou,mvir_fangzhou=get_fangzhou_radii(sim,a_array,get_all=False,D200=D200)

ok_fangzhou_all.append(ok_fangzhou)
R12_all.append(r12_fangzhou)
Rvir_all.append(rvir_fangzhou)
Mvir_all.append(mvir_fangzhou)
Rvir_halo_all.append(rvir_halo)
Mvir_halo_all.append(mvir_halo)
xi_all.append(r12_fangzhou/rvir_fangzhou)

# CARRY OUT FITS FOR EACH SNAPSHOTS
gl=fit.do_fits(gl,rvir_fangzhou,mvir_fangzhou,rmax_fit,rmin_fit,m_constraint,components=components)
          
# DEFINE DIFFERENT QUANTITIES PER OUTPUT
glsim = [ss for ss in gl[1:]]
times = [float(ss['t'].in_units('Gyr')) for ss in gl]
fmerger = array([nan]+[get_f(gl,ss, xi_merger,Rtype='ratio',com=com) for ss in glsim])
rms = [nan]
for ss in glsim:
    if 'rms' in ss['all'].keys(): rms.append(ss['all']['rms'])
    else: rms.append(nan)
rms=array(rms)

Dsnapshots=ones(size(gl))
Dsnapshots[0]=nan
for i in range(size(glsim)):
    rfit = gl[i+1]['d']['r']
    Rvir = rvir_fangzhou[i+1]
    rr=rfit[(rfit>=rmin_fit*Rvir)&(rfit<=rmax_fit*Rvir)]
    if size(gl[i]['d']['brho'])>0:
        brhov = gl[i]['d']['brho'][-1]
    else:
        brhov = nan
    p1=gl[i]['d'][fitname]['p']
    p2=gl[i+1]['d'][fitname]['p']
    fit1=prf.brho(rr, p1, model='an')
    fit2=prf.brho(rr, p2, model='an')
    try:
        Dsnapshots[i+1]=fit.Delta_area(log10(rr/Rvir),log10(fit1/brhov),log10(fit2/brhov),xlimits=delta_xlim,ymin=delta_ymin)
    except:
        Dsnapshots[i+1]=nan

relevent_fangzhou=ones(size(rvir_fangzhou),dtype=bool)
for i in range(size(rvir_fangzhou)):
    if isnan(rvir_fangzhou[i]):
        relevent_fangzhou[i]=False
        relevent_fangzhou[i+1]=False

delta=nan*ones(size(gl))
delta_i=nan*ones(size(gl))
delta_if=nan*ones(size(gl))
fmean=nan*ones(size(gl))
fstd=nan*ones(size(gl))