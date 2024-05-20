PFIM 4.0  
 
Option: 1 

 
Project: OPT
 
Date: Thu Mar 28 00:10:48 2024
 

 
**************************** INPUT SUMMARY ********************************
 
Analytical function model:  
 
dose/V * (exp(-Cl/V * t)) 
 

 
Initial design: 

 
Sample times for response: A 
       Protocol subjects doses
1 c=(0.5, 1, 5)        3   0.9

 
Total number of samples: 9
 
Associated criterion value: 335140.5
 
Identical sampling times for each response: TRUE
 
Random effect model: Trand =  1
 
Variance error model response A : ( 0.0756271255506821 + 0 *f)^2

 

Optimization step:  
 
Sampling windows for the response: A 
Window 1 : t= 0.5 1 5 8 12 16 20 
    Nb of sampling points to be taken in this window, n[ 1 ]= 3 
Maximum total number of points in one elementary protocol : 3 
Minimum total number of points in one elementary protocol : 3 

 

Now evaluating the Fisher Information Matrix for the 35 protocols generated 

 
BEST ONE GROUP PROTOCOL: 
 
Sample times for response: A 
                           times freq
1 c(`1` = 0.5, `3` = 5, `4` = 8)    1
  Subjects doses
1        3   0.9

 
Associated criterion: 380836.6
 

 
**************************** OPTIMISED DESIGN *****************************
 

 
Optimised design: 
Sample times for response: A 
         times      freq Subjects doses
1  c(5, 8, 12) 0.3782433  1.13473   0.9
2 c(0.5, 1, 8) 0.6217567  1.86527   0.9

 

 
Associated optimised criterion: 447552.6
 

 
Computation of the Population Fisher information matrix: option =  1
 
******************* FISHER INFORMATION MATRIX ******************
 
         [,1]       [,2]       [,3]
[1,] 5620.834   5408.326       0.00
[2,] 5408.326 108715.838       0.00
[3,]    0.000      0.000 8012697.75
[4,]    0.000      0.000 5771113.32
[5,]    0.000      0.000   67071.97
           [,4]        [,5]
[1,]          0       0.000
[2,]          0       0.000
[3,]    5771113   67071.967
[4,] 2350654625 1394882.238
[5,]    1394882    2975.651

 

 
************************** EXPECTED STANDARD ERRORS ************************
 
------------------------ Fixed Effects Parameters -------------------------
 
         Beta    StdError      RSE  
V  0.49631050 0.013669440 2.754211 %
Cl 0.05937912 0.003108169 5.234448 %

 

 
------------------------- Variance of Inter-Subject Random Effects ----------------------
 
         omega2     StdError      RSE  
V  3.378473e-05 4.044730e-04 1197.207 %
Cl 6.030134e-07 2.503586e-05 4151.792 %

 
------------------------ Standard deviation of residual error ------------------------ 
 
                Sigma   StdError      RSE  
sig.interA 0.07562713 0.02468227 32.63679 %

 
******************************* DETERMINANT ********************************
 
1.795644e+28
 
******************************** CRITERION *********************************
 
447552.6
 

 

 
******************* EIGENVALUES OF THE FISHER INFORMATION MATRIX ******************
 
        FixedEffects VarianceComponents
min     7.998987e+06           1641.354
max     2.350670e+09        7998987.051
max/min 2.938709e+02           4873.408

 
******************* CORRELATION MATRIX ******************
 
           [,1]       [,2]       [,3]
[1,]  1.0000000 -0.2187843  0.0000000
[2,] -0.2187843  1.0000000  0.0000000
[3,]  0.0000000  0.0000000  1.0000000
[4,]  0.0000000  0.0000000  0.2444136
[5,]  0.0000000  0.0000000 -0.4855857
           [,4]       [,5]
[1,]  0.0000000  0.0000000
[2,]  0.0000000  0.0000000
[3,]  0.2444136 -0.4855857
[4,]  1.0000000 -0.5657604
[5,] -0.5657604  1.0000000


 
Time difference of 0.7022421 secs
sys.self 
   0.117 

 
