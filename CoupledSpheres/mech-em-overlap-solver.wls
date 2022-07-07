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







(* COORDINATE ROTATION *)
(* A is a vector in spherical coordinates A = {Ar, A\[Theta], A\[Phi]} defined at point p
p is a point in spherical coordinates p = {r, theta, phi} *)
spherevec2cartvec[A_,p_]:={A[[1]]*Sin[p[[2]]]*Cos[p[[3]]]+A[[2]]*Cos[p[[2]]]*Cos[p[[3]]]-A[[3]]*Sin[p[[3]]],A[[1]]*Sin[p[[2]]]*Sin[p[[3]]]+A[[2]]*Cos[p[[2]]]*Sin[p[[3]]]+A[[3]]*Cos[p[[3]]],A[[1]]*Cos[p[[2]]]-A[[2]]*Sin[p[[2]]]}

(* X is a vector in cartesian coordinates X = {X, Y, Z} defined at point p
p is a point in spherical coordinates p = {r, theta, phi} *)
cartvec2spherevec[X_,p_]:={X[[1]]*Sin[p[[2]]]*Cos[p[[3]]]+X[[2]]*Sin[p[[2]]]*Sin[p[[3]]]+X[[3]]*Cos[p[[2]]],X[[1]]*Cos[p[[2]]]*Cos[p[[3]]]+X[[2]]*Cos[p[[2]]]*Sin[p[[3]]]-X[[3]]*Sin[p[[2]]],-X[[1]]*Sin[p[[3]]]+X[[2]]*Cos[p[[3]]]}

(* write a cartesian rotation function *)
rotateCartesianVector[vector_,eulerangles_]:=EulerMatrix[eulerangles].vector

(* write a function to rotate a spherical vector *)
rotateSphericalVector[vector_,currentpoint_,eulerangles_]:=cartvec2spherevec[rotateCartesianVector[spherevec2cartvec[vector,currentpoint],eulerangles],currentpoint]








(* GET THE MECHANICAL MODE DATA FROM DIRECTORY *)
MODES =Quiet[ToExpression[Quiet[Import[StringJoin[directory, filename]]]]];
waveNumbersK={{0},{0},Table[Quiet["wavenumberk"/.(StringJoin["m",ToString[2],"n",ToString[n],"l",ToString[0]]/.MODES[[1]])],{n, 0, 10}]};
waveNumbersP={{0},{0},Table[Quiet["wavenumberp"/.(StringJoin["m",ToString[2],"n",ToString[n],"l",ToString[0]]/.MODES[[1]])],{n, 0, 10}]};



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
garbageVariableUsedToDefineFunctions=rotateSphericalVector[displacementVecQnotRotated[m,n,l,a,c0,d0,c1,d1,c2,d2,#1,#2,#3]&@@ToSphericalCoordinates[Inverse[EulerMatrix[{ang1,ang2,ang3}]].FromSphericalCoordinates[{r,\[Theta],\[Phi]}]],{r,\[Theta],\[Phi]},{ang1,ang2,ang3}];
displacementVecQ[m_,n_,l_,a_,c0_,d0_,c1_,d1_,c2_,d2_,r_,\[Theta]_,\[Phi]_,ang1_,ang2_,ang3_]:=Evaluate[garbageVariableUsedToDefineFunctions]
Clear[garbageVariableUsedToDefineFunctions]




Print["defined mechanical modes..."]





(* WRITE DOWN EM MODES OF THE CAVITY *)
(* TE frequencies *)
wavenumberTE[m_,n_,p_,a_]:=Inactivate[SphericalBesselJZero[n,p]]/a
omegaTE[m_,n_,p_,a_]:=SPEEDOFLIGHT*wavenumberTE[m,n,p,a]
(* TM frequencies *)
wavenumberTM[m_,n_,p_,a_]:=Inactivate[SphericalBesselJDerivZero[n,p]]/a
omegaTM[m_,n_,p_,a_]:=SPEEDOFLIGHT*wavenumberTM[m,n,p,a]
(* all frequencies *)
wavenumberEM[m_,n_,p_,a_,tmorte_]:=If[tmorte==1,wavenumberTM[m,n,p,a],wavenumberTE[m,n,p,a]]
omegaEM[m_,n_,p_,a_,tmorte_]:=If[tmorte==1,omegaTM[m,n,p,a],omegaTE[m,n,p,a]]

(* Transverse electric vector potential *)
vectorPotentialTE[amplitude_,parity_, m_,n_,p_,a_,k_,r_,\[Theta]_,\[Phi]_]:={amplitude*r*SphericalBesselJ[n,k*r]*LegendreP[n,m,Cos[\[Theta]]]*If[parity==1,Cos[m*\[Phi]],Sin[m*\[Phi]]],0,0}
(* Transverse magnetic vector potential *)
vectorPotentialTM[amplitude_,parity_,m_,n_,p_,a_,k_,r_,\[Theta]_,\[Phi]_]:={amplitude*r*SphericalBesselJ[n,k*r]*LegendreP[n,m,Cos[\[Theta]]]*If[parity==1,Cos[m*\[Phi]],Sin[m*\[Phi]]],0,0}

(* calculate electric and magnetic fields for TE modes *)
garbageVariableUsedToDefineFunctions=Simplify[Evaluate[-Curl[vectorPotentialTE[amplitude,parity, m,n,p,a,k,r,\[Theta],\[Phi]],{r,\[Theta],\[Phi]},"Spherical"]]];
electricFieldTEnotrotated[amplitude_,parity_, m_,n_,p_,a_,k_,r_,\[Theta]_,\[Phi]_]:=Evaluate[garbageVariableUsedToDefineFunctions]
Clear[garbageVariableUsedToDefineFunctions]
garbageVariableUsedToDefineFunctions=Simplify[Module[{r0,\[Theta]0,\[Phi]0},Evaluate[(1/k)*Curl[Curl[vectorPotentialTE[amplitude,parity, m,n,p,a,k,r0,\[Theta]0,\[Phi]0],{r0,\[Theta]0,\[Phi]0},"Spherical"],{r0,\[Theta]0,\[Phi]0},"Spherical"]]/.{r0->r,\[Theta]0->\[Theta],\[Phi]0->\[Phi]}]];
magneticFieldTEnotrotated[amplitude_,parity_, m_,n_,p_,a_,k_,r_,\[Theta]_,\[Phi]_]:=Evaluate[garbageVariableUsedToDefineFunctions]
Clear[garbageVariableUsedToDefineFunctions]

(* calculate electric and magnetic fields for TE modes *)
garbageVariableUsedToDefineFunctions=Simplify[Module[{r0,\[Theta]0,\[Phi]0},Evaluate[(1/k)*Curl[Curl[vectorPotentialTM[amplitude,parity, m,n,p,a,k,r,\[Theta],\[Phi]],{r,\[Theta],\[Phi]},"Spherical"],{r,\[Theta],\[Phi]},"Spherical"]]/.{r0->r,\[Theta]0->\[Theta],\[Phi]0->\[Phi]}]];
electricFieldTMnotrotated[amplitude_,parity_, m_,n_,p_,a_,k_,r_,\[Theta]_,\[Phi]_]:=Evaluate[garbageVariableUsedToDefineFunctions]
Clear[garbageVariableUsedToDefineFunctions]
garbageVariableUsedToDefineFunctions=Simplify[Module[{r0,\[Theta]0,\[Phi]0},Evaluate[Curl[vectorPotentialTM[amplitude,parity, m,n,p,a,k,r,\[Theta],\[Phi]],{r,\[Theta],\[Phi]},"Spherical"]]/.{r0->r,\[Theta]0->\[Theta],\[Phi]0->\[Phi]}]];
magneticFieldTMnotrotated[amplitude_,parity_, m_,n_,p_,a_,k_,r_,\[Theta]_,\[Phi]_]:=Evaluate[garbageVariableUsedToDefineFunctions]
Clear[garbageVariableUsedToDefineFunctions]

(* calculate the rotated versions of these fields using the same method as for the mechanical mode rotations *)
garbageVariableUsedToDefineFunctions=Evaluate[rotateSphericalVector[electricFieldTEnotrotated[amplitude,parity, m,n,p,a,k,#1,#2,#3]&@@ToSphericalCoordinates[Inverse[EulerMatrix[{ang1,ang2,ang3}]].FromSphericalCoordinates[{r,\[Theta],\[Phi]}]],{r,\[Theta],\[Phi]},{ang1,ang2,ang3}]];
electricFieldTE[amplitude_,parity_, m_,n_,p_,a_,k_,r_,\[Theta]_,\[Phi]_,ang1_,ang2_,ang3_]:=Evaluate[garbageVariableUsedToDefineFunctions]
Clear[garbageVariableUsedToDefineFunctions]

garbageVariableUsedToDefineFunctions=Evaluate[rotateSphericalVector[magneticFieldTEnotrotated[amplitude,parity, m,n,p,a,k,#1,#2,#3]&@@ToSphericalCoordinates[Inverse[EulerMatrix[{ang1,ang2,ang3}]].FromSphericalCoordinates[{r,\[Theta],\[Phi]}]],{r,\[Theta],\[Phi]},{ang1,ang2,ang3}]];
magneticFieldTE[amplitude_,parity_, m_,n_,p_,a_,k_,r_,\[Theta]_,\[Phi]_,ang1_,ang2_,ang3_]:=Evaluate[garbageVariableUsedToDefineFunctions]
Clear[garbageVariableUsedToDefineFunctions]

garbageVariableUsedToDefineFunctions=Evaluate[rotateSphericalVector[electricFieldTMnotrotated[amplitude,parity, m,n,p,a,k,#1,#2,#3]&@@ToSphericalCoordinates[Inverse[EulerMatrix[{ang1,ang2,ang3}]].FromSphericalCoordinates[{r,\[Theta],\[Phi]}]],{r,\[Theta],\[Phi]},{ang1,ang2,ang3}]];
electricFieldTM[amplitude_,parity_, m_,n_,p_,a_,k_,r_,\[Theta]_,\[Phi]_,ang1_,ang2_,ang3_]:=Evaluate[garbageVariableUsedToDefineFunctions]
Clear[garbageVariableUsedToDefineFunctions]

garbageVariableUsedToDefineFunctions=Evaluate[rotateSphericalVector[magneticFieldTMnotrotated[amplitude,parity, m,n,p,a,k,#1,#2,#3]&@@ToSphericalCoordinates[Inverse[EulerMatrix[{ang1,ang2,ang3}]].FromSphericalCoordinates[{r,\[Theta],\[Phi]}]],{r,\[Theta],\[Phi]},{ang1,ang2,ang3}]];
magneticFieldTM[amplitude_,parity_, m_,n_,p_,a_,k_,r_,\[Theta]_,\[Phi]_,ang1_,ang2_,ang3_]:=Evaluate[garbageVariableUsedToDefineFunctions]
Clear[garbageVariableUsedToDefineFunctions]

(* define functions for total electric and magnetic fields *)
electricField[amplitude_,parity_,TMorTE_, m_,n_,p_,a_,k_,r_,\[Theta]_,\[Phi]_,ang1_,ang2_,ang3_]:=If[TMorTE==1,electricFieldTM[amplitude,parity, m,n,p,a,k,r,\[Theta],\[Phi],ang1,ang2,ang3],electricFieldTE[amplitude,parity, m,n,p,a,k,r,\[Theta],\[Phi],ang1,ang2,ang3]]
magneticField[amplitude_,parity_,TMorTE_, m_,n_,p_,a_,k_,r_,\[Theta]_,\[Phi]_,ang1_,ang2_,ang3_]:=If[TMorTE==1,magneticFieldTM[amplitude,parity, m,n,p,a,k,r,\[Theta],\[Phi],ang1,ang2,ang3],magneticFieldTE[amplitude,parity, m,n,p,a,k,r,\[Theta],\[Phi],ang1,ang2,ang3]]







Print["defined EM modes..."]







(* GET EM MODE NORMALIZATION FACTOR *)
getNormalizationFactor[parity_,mem_,nem_,pem_,a_,tmorte_,k_]:=(*Sqrt[((4/3)*Pi*a^2)/NIntegrate[Norm[If[tmorte==1,electricFieldTM[1.0,parity,mem,nem,pem,a,k,r,theta,phi],electricFieldTE[1.0,parity,mem,nem,pem,a,k,r,theta,phi]]]^2*r^2*Sin[theta],{r,0.0,a},{phi,0.0,2.0*Pi},{theta,0.0,Pi},Method\[Rule]"AdaptiveMonteCarlo",MaxPoints->10^5]]*)1.0





Print["defined normalization factor..."]





(* WRITE DOWN POTENTIAL ENERGY OF EM MODES *)
calculateEMFieldEnergy[TMorTE_,amplitude_,parity_,m_,n_,p_,a_,k_]:=(1/2)*NIntegrate[r^2*Sin[\[Theta]]*Norm[electricField[amplitude,parity,TMorTE,m,n,p,a,k,r,\[Theta],\[Phi],0.0,0.0,0.0]]^2*getNormalizationFactor[parity,m,n,p,a,TMorTE,k]^2,{r,0.0,a},{\[Theta],0.0,Pi},{\[Phi],0.0,2.0*Pi},Method->"AdaptiveMonteCarlo",MaxPoints->10^4,PrecisionGoal->2]






Print["defined energy..."]





(* WRITE DOWN THE OVERLAP INTEGRAND *)
(* mVec = {inclination, radial, azimuth} numbers for mechanical modes
iVec = {azimuth, inclination, radial} numbers for EM mode i 1
jVec= ----------------''------------------------
a  = innder radius of sphere
modeCoefs = {c0, d0, c1, d1, c2, d2} for mechanical modes
EMamplitude = amplitude of electric field for both fields, should be set to 1
TMorTEField = 1 for TM field and 0 for TE field
parityOfField = 1 for even parity 0 for odd (I get the stupidity here)
ki, kj, wi, wj = wave numbers and frequencies for the waves
theta, phi = overall integration variables
offsets = three sets of three Euler angles to define coordinate rotations for each of the EM fields and the mechanical mode *)
overlapIntegrand[mVec_,iVec_,jVec_,a_,modeCoefs_,EMamplitude_,TMorTEFeildi_,TMorTEFeildj_,parityOfFieldi_,parityOfFieldj_,ki_,kj_,wi_,wj_,\[Theta]_,\[Phi]_,offsets_]:=(a^2*Sin[\[Theta]] * (displacementVecQ[mVec[[1]],mVec[[2]],mVec[[3]],a,modeCoefs[[1]],modeCoefs[[2]],modeCoefs[[3]],modeCoefs[[4]],modeCoefs[[5]],modeCoefs[[6]],a,\[Theta],\[Phi],offsets[[1]],offsets[[2]],offsets[[3]]].{1,0,0}) * ((wi/wj) * magneticField[EMamplitude,parityOfFieldj,TMorTEFeildj, jVec[[1]], jVec[[2]], jVec[[3]],a,kj,a,\[Theta],\[Phi],offsets[[4]],offsets[[5]],offsets[[6]]].magneticField[EMamplitude,parityOfFieldi,TMorTEFeildi,iVec[[1]], iVec[[2]], iVec[[3]],a,ki,a,\[Theta],\[Phi],offsets[[4]],offsets[[5]],offsets[[6]]]- electricField[EMamplitude,parityOfFieldj,TMorTEFeildj,jVec[[1]], jVec[[2]], jVec[[3]],a,kj,a,\[Theta],\[Phi],offsets[[4]],offsets[[5]],offsets[[6]]].electricField[EMamplitude,parityOfFieldi,TMorTEFeildi,iVec[[1]], iVec[[2]], iVec[[3]],a,ki,a,\[Theta],\[Phi],offsets[[4]],offsets[[5]],offsets[[6]]]))

overlapfunction[modes_,a_,TMorTEFeildi_,TMorTEFeildj_,parityOfFieldi_,parityOfFieldj_,ki_,kj_,wi_,wj_,offsets_]:=((4*Pi*a^3/3)^(1/3)/(2*calculateEMFieldEnergy[TMorTEFeildi,1.0,parityOfFieldi,modes[[4]],modes[[5]],modes[[6]],a,ki])) *NIntegrate[Re[Activate[overlapIntegrand[{modes[[1]],modes[[2]],modes[[3]]},{modes[[4]],modes[[5]],modes[[6]]},{modes[[7]],modes[[8]],modes[[9]]},a,{"C0","D0",0.0,0.0,"C2","D2"}/.Quiet[(StringJoin["m",ToString[modes[[1]]],"n",ToString[modes[[2]]],"l",ToString[modes[[3]]]]/.MODES[[1]])],1.0,TMorTEFeildi,TMorTEFeildj,parityOfFieldi,parityOfFieldj,ki,kj,wi,wj,\[Theta],\[Phi],offsets]]], {\[Theta],0.0,Pi},{\[Phi],-Pi,Pi}, Method->{"AdaptiveMonteCarlo"},MaxPoints->10^4,PrecisionGoal->2]*getNormalizationFactor[parityOfFieldi,modes[[4]],modes[[5]],modes[[6]],a,TMorTEFeildi,ki]*getNormalizationFactor[parityOfFieldj,modes[[7]],modes[[8]],modes[[9]],a,TMorTEFeildj,kj]

overlap1Sphere2Modes[modes_,a_,TMorTEFeildi_,TMorTEFeildj_,parityOfFieldi_,parityOfFieldj_,offsets_]:=overlapfunction[modes,a,TMorTEFeildi,TMorTEFeildj,parityOfFieldi,parityOfFieldj,Activate[wavenumberEM[modes[[4]],modes[[5]],modes[[6]],a,TMorTEFeildi]],Activate[wavenumberEM[modes[[7]],modes[[8]],modes[[9]],a,TMorTEFeildj]],Activate[wavenumberEM[modes[[4]],modes[[5]],modes[[6]],a,TMorTEFeildi]],Activate[wavenumberEM[modes[[7]],modes[[8]],modes[[9]],a,TMorTEFeildj]],offsets]







Print["defined overlap function..."]









(* GET FREQUENCY SPLITTINGS *)







Print["defined frequency splittings"]










(* WRITE DOWN OVERLAP FUNCTION FOR COUPLED SPHERES *)
overlapCoupledSpheres[modes_,a_,TMorTEFeild_,parityOfField_,d_,offsets_]:=overlapfunction[{modes[[1]],modes[[2]],modes[[3]],modes[[4]],modes[[5]],modes[[6]],modes[[4]],modes[[5]],modes[[6]]},a,TMorTEFeild,TMorTEFeild,parityOfField,parityOfField,(wavenumberEM[modes[[4]],modes[[5]],modes[[6]],a,TMorTEFeild]-0.001*wavenumberEM[modes[[4]],modes[[5]],modes[[6]],a,TMorTEFeild])/SPEEDOFLIGHT,
(wavenumberEM[modes[[4]],modes[[5]],modes[[6]],a,TMorTEFeild]+0.001*wavenumberEM[modes[[4]],modes[[5]],modes[[6]],a,TMorTEFeild])/SPEEDOFLIGHT,(wavenumberEM[modes[[4]],modes[[5]],modes[[6]],a,TMorTEFeild]-0.001*wavenumberEM[modes[[4]],modes[[5]],modes[[6]],a,TMorTEFeild]),
(wavenumberEM[modes[[4]],modes[[5]],modes[[6]],a,TMorTEFeild]+0.001*wavenumberEM[modes[[4]],modes[[5]],modes[[6]],a,TMorTEFeild]),{offsets[[1]],offsets[[2]],offsets[[3]],offsets[[4]],offsets[[5]],offsets[[6]],offsets[[4]],offsets[[5]],offsets[[6]]}]+overlapfunction[{modes[[1]],modes[[2]],modes[[3]],modes[[4]],modes[[5]],modes[[6]],modes[[4]],modes[[5]],modes[[6]]},a,TMorTEFeild,TMorTEFeild,parityOfField,parityOfField,(wavenumberEM[modes[[4]],modes[[5]],modes[[6]],a,TMorTEFeild]-0.001*wavenumberEM[modes[[4]],modes[[5]],modes[[6]],a,TMorTEFeild])/SPEEDOFLIGHT,
(wavenumberEM[modes[[4]],modes[[5]],modes[[6]],a,TMorTEFeild]+0.001*wavenumberEM[modes[[4]],modes[[5]],modes[[6]],a,TMorTEFeild])/SPEEDOFLIGHT,(wavenumberEM[modes[[4]],modes[[5]],modes[[6]],a,TMorTEFeild]-0.001*wavenumberEM[modes[[4]],modes[[5]],modes[[6]],a,TMorTEFeild]),
(wavenumberEM[modes[[4]],modes[[5]],modes[[6]],a,TMorTEFeild]+0.001*wavenumberEM[modes[[4]],modes[[5]],modes[[6]],a,TMorTEFeild]),{offsets[[1]],offsets[[2]],offsets[[3]],offsets[[7]],offsets[[8]],offsets[[9]],offsets[[7]],offsets[[8]],offsets[[9]]}]








Print["defined coupled overlap function"]



(* DO SOME BITCHIN' CALCULATIONS *)
overlapCoupledSpheres[{2,0,0,0,1,2},innerRadius,0,1,0.1,{0,0,0,0,0,0,0,0,0}]//AbsoluteTiming



(* make a list of all possible mode combinations *)
wholeDataArray = ParallelTable[{theta1, theta2, overlapCoupledSpheres[{2,0,0,0,1,2},innerRadius,0,1,0.1,{0,theta1,theta2,0,-Pi/4.0,0,0,2.0*Pi - 3.0*Pi/4.0,0}]}, {theta1,0.0,2.0*Pi,0.3}, {theta2,0.0,2.0*Pi,0.3}];

Print[wholeDataArray]

Export[StringJoin[directory,"overlaps/overlaps-mech-EM-fixedEM-vary-mech-z-and-xy.txt"], wholeDataArray, "CSV"];
