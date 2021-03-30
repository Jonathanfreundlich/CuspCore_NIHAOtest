<div style="text-align:center"><b><span style="color:rgb(0,0,0)">FITTING THE DENSITY PROFILES WITH fitting.do_fits()</span></b><br>
</div>
<div>
<div><b><br>
</b>
<div style="text-align:center">######################################<br>
<br>
<div style="text-align:left"><br>
</div>
<div style="text-align:left"><b><u>Fields <span style="color:rgb(0,0,0)">added:</span></u><br>
</b></div>
<div style="text-align:left"><b><br>
# PROFILE DATA WITH derive_slopes()</b><br>
<br>
<table border="1" bordercolor="#888" cellspacing="0" style="border-collapse:collapse;border-color:rgb(136,136,136);border-width:1px">
<tbody>
<tr>
<td style="width:336px;height:291px">&nbsp;data.update({['all','d','g','s']:<br>
<span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span>{['lsfit_brho_b2_g3_(un)constrained']:&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <br>
<span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span>{'p':p,&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <br>
<span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; </span>'brho':brho,&nbsp;&nbsp;&nbsp; <br>
<span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; '</span>rho':rho,&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <br>
 <span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; </span>'bs':bs,&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <br>
<span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; </span>'s':s,&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <br>
 <span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; </span>'M':M,&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <br>
<span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; </span>'V':V,&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <br>
<span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; </span>'s1':prf.s(xi1*Rvir, p, model),&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <br>
 <span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; </span>'s2':prf.s(xi2*Rvir, p, model),&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <br>
<span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; </span>'bs1':prf.bs(xi1*Rvir, p, model),&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <br>
<span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; </span>'bs2':prf.bs(xi2*Rvir, p, model),&nbsp;&nbsp; <br>
<span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; </span>'c2':c2,&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <br>
<span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; </span>'rms':rms,&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <br>
<span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; </span>'Merr':(prf.calc_total_mass(r, rho, Rvir)-<span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; </span>Mvir)/Mvir}}}&nbsp;&nbsp;&nbsp; <br>
<br>
</td>
<td style="width:281px;height:291px">&nbsp;<br>
<br>
# Fit parameters (c, a, b, g, Rvir, Mvir)<br>
# Modeled brho (prf.brho(r, p, model))<br>
 # Modeled rho (prf.rho(r, p, model))<br>
# Modeled bs (prf.bs(r, p, model))<br>
# Modeled s (prf.s(r, p, model))<br>
# Modeled M (prf.M(r, p, model))<br>
# Modeled V (V = prf.V(r, p, model))<br>
# Slope at (default) 0.015 Rvir<br>
# Slope at (default) Rvir<br>
# bs slope at (default) 0.015 Rvir<br>
# bs slope at (default) Rvir<br>
# Concentration (c2 = Rvir/r[bisect(bs,2)-1])<br>
# RMS error between modeled brho and data<br>
# Relative Mvir error <br>
<br>
<br>
</td>
</tr>
</tbody>
</table>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <b><br>
&nbsp;</b><br>
<div style="text-align:center">######################################<br>
</div>
<b><b><b><br>
</b></b></b></div>
<div style="text-align:left"><u><b>Notes: </b></u><b><br>
</b><i><br>
- The fields are based on the quantities defined in gl, ie, with reduced range if reduce_range_gl() was carried out<br>
<br>
- When the reduced range is [], the outputs are also []</i><b><br>
<br>
</b><br>
<div style="text-align:center"><b>######################################<br>
<br>
</b>
<div style="text-align:left"><b><u>Usage:<br>
</u></b><br>
<font size="1"><span>fitting.do_fits(gl,rvir_fangzhou,mvir_fangzhou,rmax_fit,rmin_fit,m_constraint)</span></font><br>
<br>
<br>
</div>
</div>
</div>
</div>
</div>
</div>
