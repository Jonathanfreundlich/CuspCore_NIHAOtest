<div>
<div style="text-align:center"><b><span style="color:rgb(0,0,0)">CREATION OF A DATA FILE CONTAINING DIFFERENT PROFILES WITH format_functions.load_or_create_gl()</span><br>
<br>
</b>######################################<br>
</div>
<br>
<b><u>Structure of the data files: </u><br>
</b><br>
<b># PROFILE DATA WITH getdata()</b><br>
<br>
<table border="1" bordercolor="#888" cellspacing="0" style="border-collapse:collapse;border-color:rgb(136,136,136);border-width:1px">
<tbody>
<tr>
<td style="width:372px;height:384px">&nbsp;data = { 'sim':sim,<br>
<span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp;&nbsp; </span>'a':a,&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <br>
&nbsp;<span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; </span>'z':z,&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <br>
<span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span>'t':t,&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <br>
<span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp;&nbsp; </span>'Rvir':Rvir,&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <br>
<span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp;&nbsp; </span>'rho_crit':rho_crit,&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <br>
<span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp;&nbsp; </span>'fsv':fsv,&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <br>
<span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp;&nbsp; </span>'er':er,&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <br>
<span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp;&nbsp; </span>'Mstar':Mstar_tot,&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <br>
<span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp;&nbsp; </span>'SFR':SFR_tot,&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <br>
<span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp;&nbsp; </span>'SFR_Ha':SFR_Ha_tot,&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <br>
<span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp;&nbsp; </span>'SFR_UV':SFR_UV_tot,&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <br>
<span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp;&nbsp; </span>['all','d','g','s']:<br>
<span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; </span>{'r':p['rbins'].in_units('kpc'),<br>
<span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; </span>'M':p['mass_enc'].in_units('Msol'),<br>
<span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; </span>'n':cumsum(p['n']),&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; &nbsp; &nbsp; &nbsp;&nbsp; <br>
<span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; </span>'Rvir':Rvir,<br>
<span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; </span>'Mvir':float(sum(h1['mass'].in_units('Msol'))), <br>
<span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; </span>'eps':min(s['eps']),<br>
<span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; </span>'rho':p['density'].in_units('Msol kpc^-3'), <span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; </span>'vr_disp':p['vr_disp'].in_units('km s^-1'), <br>
<span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; </span>'vt_disp':p['vt_disp'].in_units('km s^-1'),<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; &nbsp; 'beta':p['beta']}} <br>
</td>
<td style="width:245px;height:384px"># simulation name <br>
# scale factor<br>
# redshift<br>
# time<br>
# virial radius (imprecise)<br>
# critical density<br>
<font size="1"><font size="2"># </font>sum(h1.s['mass'])/(sum(h1['mass'])*fbar)</font><br>
# h1.properties['Rmax']/max(h1['r'])<br>
# stellar mass within 0.15 Rvir<br>
# SFR [40-80 Myr]<br>
# SFR [30-60 Myr]<br>
# SFR [80-120 Myr]<br>
<br>
# radii<br>
# enclosed mass<br>
# number of particles enclosed<br>
# Virial radius (imprecise)<br>
# Virial mass (imprecise)<br>
# gravitational softening length<br>
# density<br>
# radial velocity dispersion<br>
# tangential velocity dispersion<br>
# anisotropy parameter<br>
<br>
</td>
</tr>
</tbody>
</table>
<br>
<b># STAR FORMATION BETWEEN TWO SNAPSHOTS WITH derive_SFRsnapshot()</b><br>
<br>
<table border="1" bordercolor="#888" cellspacing="0" style="border-collapse:collapse;border-color:rgb(136,136,136);border-width:1px">
<tbody>
<tr>
<td style="width:371px;height:16px">data.update({'DM_snapshot':massform,<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;&nbsp; 'SFR_snapshot':SFR_snapshot,<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;&nbsp; 'SSFR_snapshot':SSFR_snapshot})</td>
<td style="width:246px;height:16px"># Stellar mass formed since last snapshot<br>
# Associated SFR<br>
# Associated SSFR<br>
</td>
</tr>
</tbody>
</table>
<br>
<br>
<b># FLOW DATA WITH derive_flowdata()</b><br>
<br>
<table border="1" bordercolor="#888" cellspacing="0" style="border-collapse:collapse;border-color:rgb(136,136,136);border-width:1px">
<tbody>
<tr>
<td style="width:373px;height:16px">data.update({'flowdata':flowdata})<br>
<br>
flowdata = {&nbsp; 'r':r,&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <br>
<span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span>'z1':z1,&nbsp;&nbsp;&nbsp; <br>
<span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp;&nbsp; </span>'a1':a1,&nbsp; <br>
<span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp;&nbsp; </span>'t1':t1,&nbsp;&nbsp;&nbsp; <br>
<span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp;&nbsp; </span>'z2':z2,&nbsp; <br>
<span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp;&nbsp; </span>'a2':a2, <br>
<span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp;&nbsp; </span>'t2':t2, <br>
 <span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp;&nbsp; </span>'dt':dt,&nbsp;&nbsp;&nbsp;&nbsp; <br>
<span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp;&nbsp; </span>['all','d','g','s']:{['in','out']:<br>
<span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp;&nbsp; <span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span></span>{'vr' = mean(vr), <br>
<span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp;&nbsp; <span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; </span></span>'vr_ad' = mean(vr[ad]),&nbsp;&nbsp;&nbsp; <br>
 <span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp;&nbsp; <span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; </span></span>'vr_imp' = mean(vr[imp]),<br>
 <span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp;&nbsp; <span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; </span></span>'vcirc' = vcirc,&nbsp;&nbsp;&nbsp; <br>
 <span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp;&nbsp; <span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; </span></span>'m' = sum(pclsc['mass']), <br>
<span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp;&nbsp; <span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; </span></span>'m_ad' = sum(pclsc[ad]['mass']), <br>
<span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp;&nbsp; <span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; </span></span>'m_imp' = sum(pclsc[imp]['mass']) }}}<br>
<br>
</td>
<td style="width:244px;height:16px">&nbsp;<br>
<br>
# radii, same as data['all','d','g','s']['r']<br>
# redshift of the previous timestep<br>
# scale factor of the previous timestep<br>
# time of the previous timestep<br>
# redshift of the current timestep<br>
# scale factor of the current timestep<br>
# time of the current timestep<br>
# time difference<br>
<br>
# flow velocity (at each r -- it is a list)<br>
# adiabatic part's velocity (at each r)<br>
# impulsive part's velocity (at each r)<br>
# circular velocity (at each r)<br>
# flow mass (at each r)<br>
# adiabatic part's mass (at each r)<br>
# impulsive part's mass (at each r)<br>
<br>
</td>
</tr>
</tbody>
</table>
<b><br>
</b><br>
<div style="text-align:center"><b>######################################</b><br>
</div>
<br>
<b><u>Notes: </u><br>
<i><br>
</i></b><i>- r is a regularly log-spaced array of (default) nbins=150 bins between Rmin=s['eps'].min() and Rmax=10*Rvir<br>
<br>
<span style="font-family:arial,sans-serif">p = pynbody.analysis.profile.Profile(s, min=Rmin, max=Rmax, type='log', nbins=nbins, ndim=3)</span><br>
&nbsp; <br>
- data['all']['r'], data['d']['r'], data['g']['r'], data['s']['r'] and data['flowdata']['r'] are similar<br>
<br>
- quantities are defined for the following components: 'all','d' (dark matter),'g' (gas) and 's' (stars)<br>
<br>
- in flowdata, such quantities are separated between 'in' (incoming) and 'out' (outflowing)<br>
</i><br>
&nbsp; <br>
<div style="text-align:center"><b>######################################</b><br>
</div>
<br>
</div>
<u><b>Usage:</b></u><br>
<span><br>
import sys<br>
sys.path.insert(0, '/cs/sci/freundlich/CUSPCORE/Analysis/Formatting')<br>
import format_functions<br>
reload(format_functions)<br>
<br>
sims=['g1.37e11','g1.52e11','g1.57e11','g1.59e11','g2.19e11','g2.63e10','g8.06e11']<br>
<br>
for sim in sims:<br>
&nbsp;&nbsp;&nbsp; format_functions.load_or_create_gl(sim,directory='/cs/sci/freundlich/CUSPCORE/Analysis/DATA/',use_fangzhou_Rvir=True)<br>
<br>
</span>

