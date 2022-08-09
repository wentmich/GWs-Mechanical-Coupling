(*Limits Mathematica to requested resources*)
Unprotect[$ProcessorCount];$ProcessorCount = 20;

(*Prints the machine name that each kernel is running on*)
Print[ParallelEvaluate[$MachineName]];


argv = Rest @ $CommandLine;
argc = Length @ argv;

Print[argv]


(* CONVENTIONS
"parity" --- 1 means even parity and 0 means odd parity
"TMorTE" --- 1 means TM and 0 means TE
The indices here are rather weird. Need to be careful about the selection rules and all that.

SELECTION RULES FOR EM CAVITIES (EVEN)
0 <= m <= n azimuthal
0 <= n inclination
1 <= p radial

SELECTION RULES FOR MECHANICAL CAVITIES
m >= 0 inclination
n >= 1 radial
0 <= l <= m azimuthal *)







(* define all physical constants (not cavity parameters) *)
(* define all physical constants (not cavity parameters) *)
SPEEDOFLIGHT=1.0(*3.0*10.0^8.0;*);
shearModulus=38.0*10^9;
youngsModulus=105.0*10^9;
matieralDensity=8.57*10^(-3)/((10^(-2))^3);
(* define all derived constants *)
lambdaMaterial = -shearModulus*(youngsModulus-2*shearModulus)/(youngsModulus-3*shearModulus);
vt = (matieralDensity/shearModulus)^(-1/2);
vl = (matieralDensity/(lambdaMaterial+2*shearModulus))^(-1/2);
(* define all cavity parameters *)
innerRadius=ToExpression[argv[[3]]];
outerRadius=ToExpression[argv[[3]]] + ToExpression[argv[[4]]];
directory = argv[[5]];
filename = argv[[6]];
maxWaveNumber=10.0^3;

Print[innerRadius]
Print[outerRadius]





(* COORDINATE ROTATION *)
(* A is a vector in spherical coordinates A = {Ar, A\[Theta], A\[Phi]} defined at point p
p is a point in spherical coordinates p = {r, theta, phi} *)
spherevec2cartvec[A_,p_]:={A[[1]]*Sin[p[[2]]]*Cos[p[[3]]]+A[[2]]*Cos[p[[2]]]*Cos[p[[3]]]-A[[3]]*Sin[p[[3]]],A[[1]]*Sin[p[[2]]]*Sin[p[[3]]]+A[[2]]*Cos[p[[2]]]*Sin[p[[3]]]+A[[3]]*Cos[p[[3]]],A[[1]]*Cos[p[[2]]]-A[[2]]*Sin[p[[2]]]}

(* X is a vector in cartesian coordinates X = {X, Y, Z} defined at point p
p is a point in spherical coordinates p = {r, theta, phi} *)
cartvec2spherevec[X_,p_]:={X[[1]]*Sin[p[[2]]]*Cos[p[[3]]]+X[[2]]*Sin[p[[2]]]*Sin[p[[3]]]+X[[3]]*Cos[p[[2]]],X[[1]]*Cos[p[[2]]]*Cos[p[[3]]]+X[[2]]*Cos[p[[2]]]*Sin[p[[3]]]-X[[3]]*Sin[p[[2]]],-X[[1]]*Sin[p[[3]]]+X[[2]]*Cos[p[[3]]]}

(* write a cartesian rotation function *)
rotateCartesianVector[vector_,eulerangles_]:=RollPitchYawMatrix[eulerangles,{3,2,3}].vector

(* write a function to rotate a spherical vector *)
rotateSphericalVector[vector_,currentpoint_,eulerangles_]:=cartvec2spherevec[rotateCartesianVector[spherevec2cartvec[vector,currentpoint],eulerangles],currentpoint]








(* GET THE MECHANICAL MODE DATA FROM DIRECTORY *)
MODES =Quiet[ToExpression[Quiet[Import[StringJoin[directory, filename]]]]];
waveNumbersK=Table[Quiet["wavenumberk"/.(StringJoin["m",ToString[m],"n",ToString[n],"l",ToString[0]]/.MODES[[1]])],{m,0,3},{n, 0, 3}];
waveNumbersP=Table[Quiet["wavenumberp"/.(StringJoin["m",ToString[m],"n",ToString[n],"l",ToString[0]]/.MODES[[1]])],{m,0,3},{n, 0, 3}];

Print[MODES]

Print["loaded in mode data..."]


(* CALCULATE THE ZEROS OF ALL THE BESSEL FUNCTIONS *)
garbageVariableUsedToDefineFunctions=FullSimplify[D[x*SphericalBesselJ[n,x],x]];
sphericalBesselDerivative[n_,x_]:=Evaluate[garbageVariableUsedToDefineFunctions]
Clear[garbageVariableUsedToDefineFunctions]
garbageVariableUsedToDefineFunctions=x*SphericalBesselJ[n,x];
sphericalBessel[n_,x_]:=Evaluate[garbageVariableUsedToDefineFunctions]
Clear[garbageVariableUsedToDefineFunctions]
sphericalBesselJZero = Import["/home/wentzel4/GWs/GWs-Mechanical-Coupling/SphereModes/bessel-zeros/norm_spherical_bessel_zeros.txt", "CSV"];
sphericalBesselJDerivZero = Import["/home/wentzel4/GWs/GWs-Mechanical-Coupling/SphereModes/bessel-zeros/norm_spherical_bessel_derivatives_zeros.txt", "CSV"];
SphericalBesselJZero[n_,p_]:=sphericalBesselJZero[[n+1]][[p]]
SphericalBesselJDerivZero[n_,p_]:=sphericalBesselJDerivZero[[n+1]][[p]]



Print["got zeros of bessel functions loaded..."]



(* WRITE DOWN THE MECHANICAL MODES OF THE CAVITY *)
(* wave numbers *)
wavek[m_,n_,l_,a_]:=Inactivate[waveNumbersK[[m+1]][[n+1]]]
wavep[m_,n_,l_,a_]:=Inactivate[waveNumbersP[[m+1]][[n+1]]]
(*wavek[m_,n_,l_,a_]:= Sqrt[omegaMech[m,n,l,a]^2 / vt^2];
wavep[m_,n_,l_,a_]:= Sqrt[omegaMech[m,n,l,a]^2 / vl^2];*)
omegaMech[m_,n_,l_,a_]:=vt*Sqrt[wavek[m,n,l,a]]

Print[Activate[wavep[2,0,0,0.1]]]
Print[Activate[wavek[2,0,0,0.1]]]
Print[Activate[omegaMech[2,0,0,0.1]]]


(* constituent functions *)
aux\[Phi][m_,n_,l_,a_,r_,\[Theta]_,\[Phi]_]:=SphericalBesselJ[m,wavep[m,n,l,a]*r]*SphericalHarmonicY[m,l,\[Theta],\[Phi]];
aux\[Phi]T[m_,n_,l_,a_,r_,\[Theta]_,\[Phi]_]:=SphericalBesselY[m,wavep[m,n,l,a]*r]*SphericalHarmonicY[m,l,\[Theta],\[Phi]];
aux\[Psi][m_,n_,l_,a_,r_,\[Theta]_,\[Phi]_]:=SphericalBesselJ[m,wavek[m,n,l,a]*r]*SphericalHarmonicY[m,l,\[Theta],\[Phi]];
aux\[Psi]T[m_,n_,l_,a_,r_,\[Theta]_,\[Phi]_]:=SphericalBesselY[m,wavek[m,n,l,a]*r]*SphericalHarmonicY[m,l,\[Theta],\[Phi]];

(* mechanical mode displacement vector *)
garbageVariableUsedToDefineFunctions=c0*Grad[aux\[Phi][m,n,l,a,r,\[Theta],\[Phi]],{r,\[Theta],\[Phi]},"Spherical"] + d0*Grad[aux\[Phi]T[m,n,l,a,r,\[Theta],\[Phi]],{r,\[Theta],\[Phi]},"Spherical"]+c1*Cross[{r,0,0},Grad[aux\[Psi][m,n,l,a,r,\[Theta],\[Phi]],{r,\[Theta],\[Phi]},"Spherical"]]+d1*Cross[{r,0,0},Grad[aux\[Psi]T[m,n,l,a,r,\[Theta],\[Phi]],{r,\[Theta],\[Phi]},"Spherical"]]+c2*Curl[Cross[{r,0,0},Grad[aux\[Psi][m,n,l,a,r,\[Theta],\[Phi]],{r,\[Theta],\[Phi]},"Spherical"]],{r,\[Theta],\[Phi]},"Spherical"]+d2*Curl[Cross[{r,0,0},Grad[aux\[Psi]T[m,n,l,a,r,\[Theta],\[Phi]],{r,\[Theta],\[Phi]},"Spherical"]],{r,\[Theta],\[Phi]},"Spherical"];
displacementVecQnotRotated[m_,n_,l_,a_,c0_,d0_,c1_,d1_,c2_,d2_,r_,\[Theta]_,\[Phi]_]:=Evaluate[garbageVariableUsedToDefineFunctions]
Clear[garbageVariableUsedToDefineFunctions]

(* define a rotated displacement vector field *)
garbageVariableUsedToDefineFunctions=rotateSphericalVector[displacementVecQnotRotated[m,n,l,a,c0,d0,c1,d1,c2,d2,#1,#2,#3]&@@ToSphericalCoordinates[Inverse[RollPitchYawMatrix[{ang1,ang2,ang3},{3,2,3}]].FromSphericalCoordinates[{r,\[Theta],\[Phi]}]],{r,\[Theta],\[Phi]},{ang1,ang2,ang3}];
displacementVecQ[m_,n_,l_,a_,c0_,d0_,c1_,d1_,c2_,d2_,r_,\[Theta]_,\[Phi]_,ang1_,ang2_,ang3_]:=Evaluate[garbageVariableUsedToDefineFunctions]
Clear[garbageVariableUsedToDefineFunctions]




Print["defined mechanical modes..."]

checkNormFunction[mindices_]:=Module[{mvec=mindices},
coefs={"C0","D0",0.0,0.0,"C2","D2"}/.StringJoin["m",ToString[mvec[[1]]],"n",ToString[mvec[[2]]],"l",ToString[mvec[[3]]]]/.MODES[[1]];
NIntegrate[r^2*Sin[\[Theta]]*Norm[Re[Activate[displacementVecQnotRotated[mvec[[1]],mvec[[2]],mvec[[3]],innerRadius,coefs[[1]],coefs[[2]],coefs[[3]],coefs[[4]],coefs[[5]],coefs[[6]],r,\[Theta],\[Phi]]]]]^2,{r,innerRadius,outerRadius},{\[Theta],0.0,Pi},{\[Phi],-Pi,Pi},Method->"AdaptiveMonteCarlo",MaxPoints->10^5,PrecisionGoal->2]/((4/3)*Pi*(outerRadius^3-innerRadius^3))]


(* IMPORT ALL THE MECHANICAL MODES YOU WANT TO CALCULATE *)
modelist = Import["list-of-mechanical-modes.txt", "CSV"];
Print[modelist[[1]]]

wholeDataArray = ParallelTable[{modelist[[i,1]],modelist[[i,2]],modelist[[i,3]],checkNormFunction[modelist[[i]]]}, {i,1,30}];

Print[wholeDataArray]

Export[StringJoin[directory,"normalization-checks.csv"], wholeDataArray, "CSV"];
