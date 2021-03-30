<div style="text-align:center">
<div style="text-align:center"><b><span style="color:rgb(0,0,0)">DERIVATION OF DIFFERENT FIELDS WITH slopes_functions.derive_slopes() and prepare_functions.derive_brho() from the FORMATTED DATA FILE</span></b><br>
</div>
</div>
<div><b><br>
</b>
<div style="text-align:center">######################################<br>
<br>
<div style="text-align:left"><br>
</div>
<div style="text-align:left"><b><u>Fields added: </u><br>
</b></div>
<div style="text-align:left"><b><br>
# PROFILE DATA WITH derive_slopes()</b><b><br>
</b><br>
<table border="1" cellspacing="0" style="border-collapse:collapse;border-color:rgb(136,136,136);border-width:1px">
<tbody>
<tr>
<td style="width:372px;height:304px">data.update({['all','d','g','s']:{<br>
<span>&nbsp;&nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp;&nbsp; </span>'M':noinf(M), <br>
<span>&nbsp;&nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp;&nbsp; </span>'Mall':noinf(Mall),&nbsp;&nbsp;&nbsp; <br>
<span>&nbsp;&nbsp; &nbsp;</span><span> &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp;&nbsp; </span>'rho':10**noinf(logrho),&nbsp;&nbsp; <br>
<span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp;&nbsp; </span>'logrho':noinf(logrho), <br>
<span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span> &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp;&nbsp; </span>'rho_smooth':10**noinf(logrho_smooth), <br>
<span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp;&nbsp; </span>'logrho_smooth':noinf(logrho_smooth),&nbsp;&nbsp;&nbsp;&nbsp; <br>
<span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span> &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp;&nbsp; </span>'sigmar':noinf(sigmar),&nbsp;&nbsp;&nbsp; <br>
<span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp;&nbsp; </span>'logsigmar':noinf(logsigmar),&nbsp;&nbsp; <br>
 <span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span> &nbsp;</span><span>&nbsp;&nbsp;&nbsp; </span>'logsigmar2':noinf(logsigmar2),&nbsp;&nbsp;&nbsp;&nbsp; <br>
 <span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp;</span><span>&nbsp;&nbsp;&nbsp; </span>'logsigmar_smooth':noinf(logsigmar_smooth),&nbsp; <br>
 <span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp; </span>'logsigmar2_smooth':noinf(logsigmar2_smooth),&nbsp; <br>
 <span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp; </span>'sigmar_smooth':10**noinf(logsigmar_smooth),<br>
<span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp; </span>'sigmar2_smooth':10**noinf(logsigmar2_smooth),<br>
<span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp; </span>'alpha':noinf(alpha),&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <br>
<span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp; </span>'beta_smooth':noinf(beta_smooth),&nbsp;&nbsp;&nbsp;&nbsp; <br>
<span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp; </span>'gamma':noinf(gamma)}})&nbsp;&nbsp;&nbsp; <br>
<br>
<br>
</td>
<td style="width:245px;height:304px"><br>
# Enclosed mass (M=ss[c]['M'])<br>
# Total enclosed mass (ss['all']['M'])<br>
# Density profile<br>
# "<br>
# "<br>
# "<br>
# Velocity dispersion (ss[c]['vr_disp'])<br>
# "<br>
# "<br>
# "<br>
# "<br>
# "<br>
# "<br>
# Log slope of the density profile<br>
# Smoothed anisotropy parameter<br>
# Log slope of the sigmar^2 profile<br>
<br>
</td>
</tr>
</tbody>
</table>
<br>
<b># PROFILE DATA WITH derive_brho()</b><br>
<br>
<table border="1" bordercolor="#888" cellspacing="0" style="border-collapse:collapse;border-color:rgb(136,136,136);border-width:1px">
<tbody>
<tr>
<td style="width:214px;height:176px">&nbsp;data.update({['all','d','g','s']:{<br>
<span>&nbsp;&nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp;&nbsp; </span>'dM':dM,<br>
<span>&nbsp;&nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp;&nbsp; </span>'drho':drho,&nbsp; <br>
<span>&nbsp;&nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp;&nbsp; </span>'brho':brho,&nbsp;&nbsp;&nbsp; <br>
<span>&nbsp;&nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp;&nbsp; </span>'dbrho':dbrho,&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <br>
<span>&nbsp;&nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp;&nbsp; </span>'s':s,&nbsp;&nbsp;&nbsp; <br>
<span>&nbsp;&nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp;&nbsp; </span>'bs':bs,&nbsp;&nbsp;&nbsp;&nbsp; <br>
<span>&nbsp;&nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp;&nbsp; </span>'Rmax':Rmax,&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <br>
<span>&nbsp;&nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp;&nbsp; </span>'Mmax':Mmax,&nbsp;&nbsp; <br>
<span>&nbsp;&nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp; &nbsp;</span><span>&nbsp;&nbsp;&nbsp; </span>'Rthr':Rthr}})&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <br>
<br>
</td>
<td style="width:403px;height:176px">&nbsp;<br>
# Uncertainty on M (dM = M/sqrt(cumsum(n)))<br>
# Uncertainty on rho (drho=rho/sqrt(n))<br>
# Mean density<br>
# Uncertainty on brho (brho/sqrt(cumsum(n)))<br>
# Logarithmic slope of rho (=alpha)<br>
# Logarithmic slope of brho<br>
# Largest radius r[-1]<br>
# Largest mass M[-1]<br>
# Rthr = nan</td>
</tr>
</tbody>
</table>
<br>
<br>
</div>
######################################<br>
<br>
<div style="text-align:left"><b><u>Usage:<br>
<br>
<div></div>
</u></b><span><font size="1">import sys<br>
import pickle<br>
from numpy import *<br>
<br>
sys.path.insert(0, '/cs/sci/freundlich/CUSPCORE/Analysis/Preparing')<br>
sys.path.insert(0, '/cs/sci/freundlich/CUSPCORE/Analysis/Formatting')<br>
<br>
import treal_functions<br>
import slopes_functions<br>
import prepare_functions<br>
<br>
reload(treal_functions)<br>
reload(slopes_functions)<br>
reload(prepare_functions)<br>
<br>
sim = ['g1.08e11'][0]<br>
directory='/cs/sci/freundlich/CUSPCORE/Analysis/Formatting/Test/'<br>
directory='/cs/sci/freundlich/CUSPCORE/NIHAO_data/'<br>
<br>
with open(directory+'NIHAO-%s.pickle'%sim[1:]) as f:<br>
&nbsp;&nbsp;&nbsp; gl = pickle.load(f)<br>
&nbsp;<br>
gl = slopes_functions.derive_slopes(gl,polyorder=3,sigma = 11,mode= 'interp',double_smooth=False,rlim=[-2.,0.])<br>
gl = prepare_functions.define_brho(gl,polyorder=3,sigma = 11,mode= 'interp',double_smooth=False,rlim=[-2.,0.])<br>
treal=treal_functions.load_or_create_gl(sim)<br>
</font></span><b><u><br>
<br>
</u></b>
<div style="text-align:center">######################################<br>
</div>
</div>
<b><br>
</b></div>
<div style="text-align:center"><b>REDUCTION OF THE STUDY RANGE WITH prepare_functions.get_fitrange(), reduce_range_gl() and reduce_range_Treal</b><br>
<br>
######################################<br>
<br>
<div style="text-align:left">Cut the initial data files for the profiles and for Treal to a narrower radius range, for example for log(r/Rvir) between -2 and 0 (default range). <br>
<br>
</div>
<div style="text-align:left"><b><u>Usage:<br>
</u></b><span><font size="1"><br>
fitrange=prepare_functions.get_fitrange(gl)<br>
gl=prepare_functions.reduce_range_gl(gl,fitrange)<br>
treal=prepare_functions.reduce_range_Treal(treal,fitrange)<br>
</font></span><b><u><br>
</u></b></div>
<div style="text-align:left"><br>
<div style="text-align:center">######################################<br>
