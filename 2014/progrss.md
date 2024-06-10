### `7 June, Friday`

---
## RESULTS- 2014 paper

<span style="color:yellow">`MISTAKE IN INITIAL CONDITIONS, CHANGE 10e-2 TO 10e-3`</span>

<center><img src="/home/gayathri/MSc Thesis/Thesis-work-9th-sem/2014/data/realfresult.png" alt="Alt text" width="600" height="400"></center>
<center><img src="/home/gayathri/MSc Thesis/Thesis-work-9th-sem/2014/data/editresult.png" alt="Alt text" width="600" height="400"></center>
<center><img src="/home/gayathri/MSc Thesis/Thesis-work-9th-sem/2014/data/res1no.png" alt="Alt text" width="600" height="400"></center>
<center><img src="/home/gayathri/MSc Thesis/Thesis-work-9th-sem/2014/data/res1.png" alt="Alt text" width="600" height="400"></center>
<center><img src="/home/gayathri/MSc Thesis/Thesis-work-9th-sem/2014/data/res2no.png" alt="Alt text" width="600" height="400"></center>
<center><img src="/home/gayathri/MSc Thesis/Thesis-work-9th-sem/2014/data/res2.png" alt="Alt text" width="600" height="400"></center>

## WORK DONE
<ol>
<li>Implemented dynamic quenching.</li>
<li>Trials</li>
<ol>
<li>Other schemes for time stepping.
<li>Increased spatial resolution.
<li>Increased time steps.
<li>Slightly tuning parameters.
<li>Changing the equations- coefficient, ways to make it dimensionless etc.


next week-how to include halo, Bz/delZ,eta no more constant, start the der from scratch

</ol>
<li>Read 2016 paper.</li>
<li>Figuring out how to implement the paper.</li>
</ol>

## DOUBTS
`From the paper:`
<ul>
<li>The boundary conditions for thick disc case.
<li>Reproducing figures from 2016: thin disc? 2.5D?
<li>Dimensionless using u and td?
<li>Why is the rotation curve different?
</ul>

`Codes:`
<ul>
<li>How to improve upon 2014?
<ul>
<li>extend the grid in z in both ways, to incorporate halo + alpha fn of r + eta_t varying in z + h flared
<li>next step 2016 paper or extending 2014?
</ul>
<li>CHANG-ES and my work
<li>Rederiving the 2016 paper- a short summary if possible
<li>Incorporating halos - how to set the boundary
<li>Why am i not getting the correct plot for zero flux case, what happens to alpha_m there. The parameters are same as in paper:

| Parameter | Value in Dimensional Units | Value in Dimensionless Units |
|-----------|----------------------------|------------------------------|
| eta       | 0.3333333333333333 km kpc / s | 1.0                        |
| radius    | 4.0 kpc                    | 10.498687664041995           |
| r_d       | 4.0 kpc                    | 26.246719160104988           |
| h_d       | 0.35 kpc                   | 0.9186351706036745           |
| h         | 0.381 kpc                  | 1.0                          |
| td        | 0.42581189007421516 Gyr   | 0.977792221680789            |
| alpha_0   | 1.49071198499986 km / s    | 1.70388379885484             |
| G         | -45.6 km / (kpc s)         | -20.30904353673911           |
| omega_0   | 127.0 km / (kpc s)         | 148.45792059019814           |
| r_omega   | 2.0 kpc                    | 5.2493438320209975           |
| l         | 0.1 kpc                    | 0.26246719160104987          |
| B_0       | 8.2e-06 G                  | 1.0                          |
| R         | 20.0 kpc                   | 52.493438320209975           |
| U_0       | 1.0 km / s                 | 2.348454853681303            |
| k         | 0.1 km kpc / s             | 0.3068136495137085           |
| z_i       | -0.381 kpc                 | -1.0                         |
| z_f       | 0.381 kpc                  | 1.0                          |
| R_alpha   | 1.49071198499986 km / s    | 1.70388379885484             |
| R_omega   | 127.0 km / (kpc s)         | -20.30904353673911           |
| R_u       | 2.348454853681303          | 2.348454853681303            |
| R_k       | 0.3068136495137085         | 0.3068136495137085           |

</ul>












