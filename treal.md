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
import treal_functions<br>
<br>
sims =['g1.37e11','g1.52e11','g1.57e11','g1.59e11','g2.19e11','g2.63e10','g8.06e11']<br>
<br>
for sim in sims:<br>
&nbsp;&nbsp;&nbsp; treal_functions.load_or_create_gl(sim,directory='/cs/sci/freundlich/CUSPCORE/Analysis/DATA/',use_fangzhou_Rvir=True)</span><br>
</div>

