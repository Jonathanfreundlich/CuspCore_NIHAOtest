<div style="text-align:center"><b><span style="color:rgb(0,0,0)">CREATION OF A DATA FILE CONTAINING THE REAL KINETIC ENERGY WITH treal_functions.load_or_create_gl()</span><br>
<br>
</b>######################################<br>
<br>
<div style="text-align:left">For each snapshot i, this function creates an array Treal such that Treal[i]=[r,Rvir,Tr]. <br>
<br>
</div>
<div style="text-align:left"><u><b>Notes: </b></u><br>
<i><br>
- The radii r are taken from the files NIHAO-xxx.pickle containing the profile and flow data<br>
<br>
- It is possible to choose Rvir from the simulation haloes or from Fangzhou's data<br>
</i><br>
<div style="text-align:center">######################################<br>
<br>
<div style="text-align:left"><b><u>Usage: </u></b><br>
<span><br>
import sys<br>
sys.path.insert(0, '/cs/sci/freundlich/CUSPCORE/Analysis/Formatting')<br>
import treal_functions<br>
reload(treal_functions)<br>
<br>
sims =['g1.37e11','g1.52e11','g1.57e11','g1.59e11','g2.19e11','g2.63e10','g8.06e11']<br>
<br>
for sim in sims:<br>
&nbsp;&nbsp;&nbsp; treal_functions.load_or_create_gl(sim,directory='/cs/sci/freundlich/CUSPCORE/Analysis/DATA/',use_fangzhou_Rvir=True)</span><br>
</div>
<br>
######################################<br>
<br>
<div style="text-align:left"><b><u>Treal profiles for galaxy g1.08e11 (all snapshots):</u></b><br>
</div>
<br>
<div style="display:block;text-align:left"><a href="https://sites.google.com/site/hujicosmo/cusp-core-transition/2-kinetic-energy/Treal.png?attredirects=0" imageanchor="1"><img border="0" src="Treal.png" style="width:50%"></a></div>
<br>
</div>
</div>
<div style="text-align:left"><br>
by Jonathan Freundlich (13/03/2018)</div>
<br>
</div>
