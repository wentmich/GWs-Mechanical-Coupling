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
innerRadius=1.0;
outerRadius=1.01;
maxWaveNumber=10.0^3;







(* COORDINATE ROTATION *)
newPoint[oldpoint_,offsets_]:= ToSphericalCoordinates[EulerMatrix[offsets].FromSphericalCoordinates[oldpoint]]







(* GET THE MECHANICAL MODE DATA FROM DIRECTORY *)
MODES =Quiet[ToExpression[Quiet[Import[StringJoin[directory, filename]]]]];
waveNumbersK={{0},{0},Table[Quiet["wavenumberk"/.(StringJoin["m",ToString[2],"n",ToString[n],"l",ToString[0]]/.MODES[[1]])],{n, 0, 10}]};
waveNumbersP={{0},{0},Table[Quiet["wavenumberp"/.(StringJoin["m",ToString[2],"n",ToString[n],"l",ToString[0]]/.MODES[[1]])],{n, 0, 10}]};








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







(* WRITE DOWN THE MECHANICAL MODES OF THE CAVITY *)
(* define the necessary differential operators *)
rCrossNabla[myFnc_,r_,\[Theta]_,\[Phi]_,m_,n_,l_,a_]:={0, -(1/Sin[\[Theta]])*D[myFnc[m,n,l,a,r,\[Theta],\[Phi]],\[Phi]], D[myFnc[m,n,l,a,r,\[Theta],\[Phi]],\[Theta]]}
nablaCrossRCrossNabla[myFnc_,r_,\[Theta]_,\[Phi]_,m_,n_,l_,a_]:={(1/r)* D[myFnc[m,n,l,a,r,\[Theta],\[Phi]],\[Theta],\[Theta]] + (1/(r*Sin[\[Theta]]^2))*D[myFnc[m,n,l,a,r,\[Theta],\[Phi]],\[Phi],\[Phi]], -D[myFnc[m,n,l,a,r,\[Theta],\[Phi]],\[Theta],r]+(1/r)*D[myFnc[m,n,l,a,r,\[Theta],\[Phi]],\[Theta]], -(1/Sin[\[Theta]])*D[myFnc[m,n,l,a,r,\[Theta],\[Phi]],\[Phi],r]+(1/(r*Sin[\[Theta]]))*D[myFnc[m,n,l,a,r,\[Theta],\[Phi]],\[Phi]]}

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
garbageVariableUsedToDefineFunctions=c0*Grad[aux\[Phi][m,n,l,a,r,\[Theta],\[Phi]],{r,\[Theta],\[Phi]},"Spherical"] + d0*Grad[aux\[Phi]T[m,n,l,a,r,\[Theta],\[Phi]],{r,\[Theta],\[Phi]},"Spherical"]+c1*rCrossNabla[aux\[Psi],r,\[Theta],\[Phi],m,n,l,a]+d1*rCrossNabla[aux\[Psi]T,r,\[Theta],\[Phi],m,n,l,a]+I*c2*nablaCrossRCrossNabla[aux\[Psi],r,\[Theta],\[Phi],m,n,l,a]+I*d2*nablaCrossRCrossNabla[aux\[Psi]T,r,\[Theta],\[Phi],m,n,l,a];
displacementVecQ[m_,n_,l_,a_,c0_,d0_,c1_,d1_,c2_,d2_,r_,\[Theta]_,\[Phi]_]:=Evaluate[garbageVariableUsedToDefineFunctions]
Clear[garbageVariableUsedToDefineFunctions]







(* WRITE DOWN EM MODES OF THE CAVITY *)
(* TE frequencies *)
wavenumberTE[m_,n_,p_,a_]:=Inactivate[SphericalBesselJZero[n,p]]/a
omegaTE[m_,n_,p_,a_]:=SPEEDOFLIGHT*wavenumberTE[m,n,p,a]
(* TM frequencies *)
wavenumberTM[m_,n_,p_,a_]:=Inactivate[SphericalBesselJDerivZero[n,p]]/a
omegaTM[m_,n_,p_,a_]:=SPEEDOFLIGHT*wavenumberTM[m,n,p,a]

(* Transverse electric vector potential *)
vectorPotentialTE[amplitude_,parity_, m_,n_,p_,a_,k_,r_,\[Theta]_,\[Phi]_]:={amplitude*r*SphericalBesselJ[n,k*r]*LegendreP[n,m,Cos[\[Theta]]]*If[parity==1,Cos[m*\[Phi]],Sin[m*\[Phi]]],0,0}
(* Transverse magnetic vector potential *)
vectorPotentialTM[amplitude_,parity_,m_,n_,p_,a_,k_,r_,\[Theta]_,\[Phi]_]:={amplitude*r*SphericalBesselJ[n,k*r]*LegendreP[n,m,Cos[\[Theta]]]*If[parity==1,Cos[m*\[Phi]],Sin[m*\[Phi]]],0,0}
(* calculate electric and magnetic fields for TE modes *)
garbageVariableUsedToDefineFunctions=Evaluate[-Curl[vectorPotentialTE[amplitude,parity, m,n,p,a,k,r,\[Theta],\[Phi]],{r,\[Theta],\[Phi]},"Spherical"]];
electricFieldTE[amplitude_,parity_, m_,n_,p_,a_,k_,r_,\[Theta]_,\[Phi]_]:=Evaluate[garbageVariableUsedToDefineFunctions]
Clear[garbageVariableUsedToDefineFunctions]
garbageVariableUsedToDefineFunctions=Module[{r0,\[Theta]0,\[Phi]0},Evaluate[(1/k)*Curl[Curl[vectorPotentialTE[amplitude,parity, m,n,p,a,k,r0,\[Theta]0,\[Phi]0],{r0,\[Theta]0,\[Phi]0},"Spherical"],{r0,\[Theta]0,\[Phi]0},"Spherical"]]/.{r0->r,\[Theta]0->\[Theta],\[Phi]0->\[Phi]}];
magneticFieldTE[amplitude_,parity_, m_,n_,p_,a_,k_,r_,\[Theta]_,\[Phi]_]:=Evaluate[garbageVariableUsedToDefineFunctions]
Clear[garbageVariableUsedToDefineFunctions]
(* calculate electric and magnetic fields for TE modes *)
garbageVariableUsedToDefineFunctions=Module[{r0,\[Theta]0,\[Phi]0},Evaluate[(1/k)*Curl[Curl[vectorPotentialTM[amplitude,parity, m,n,p,a,k,r,\[Theta],\[Phi]],{r,\[Theta],\[Phi]},"Spherical"],{r,\[Theta],\[Phi]},"Spherical"]]/.{r0->r,\[Theta]0->\[Theta],\[Phi]0->\[Phi]}];
electricFieldTM[amplitude_,parity_, m_,n_,p_,a_,k_,r_,\[Theta]_,\[Phi]_]:=Evaluate[garbageVariableUsedToDefineFunctions]
Clear[garbageVariableUsedToDefineFunctions]
garbageVariableUsedToDefineFunctions=Module[{r0,\[Theta]0,\[Phi]0},Evaluate[Curl[vectorPotentialTM[amplitude,parity, m,n,p,a,k,r,\[Theta],\[Phi]],{r,\[Theta],\[Phi]},"Spherical"]]/.{r0->r,\[Theta]0->\[Theta],\[Phi]0->\[Phi]}];
magneticFieldTM[amplitude_,parity_, m_,n_,p_,a_,k_,r_,\[Theta]_,\[Phi]_]:=Evaluate[garbageVariableUsedToDefineFunctions]
Clear[garbageVariableUsedToDefineFunctions]








(* GET EM MODE NORMALIZATION FACTOR *)
getNormalizationFactor[parity_,mem_,nem_,pem_,a_,tmorte_,k_]:=(*Sqrt[((4/3)*Pi*a^2)/NIntegrate[Norm[If[tmorte==1,electricFieldTM[1.0,parity,mem,nem,pem,a,k,r,theta,phi],electricFieldTE[1.0,parity,mem,nem,pem,a,k,r,theta,phi]]]^2*r^2*Sin[theta],{r,0.0,a},{phi,0.0,2.0*Pi},{theta,0.0,Pi},Method->"AdaptiveMonteCarlo",MaxPoints->10^5]]*) 1.0







(* WRITE DOWN POTENTIAL ENERGY OF EM MODES *)
calculateEMFieldEnergy[TMorTE_,amplitude_,parity_,m_,n_,p_,a_,k_]:=(1/2)*NIntegrate[r^2*Sin[\[Theta]]*Norm[If[TMorTE==0,electricFieldTE[amplitude,parity,m,n,p,a,k,r,\[Theta],\[Phi]],electricFieldTM[amplitude,parity,m,n,p,a,k,r,\[Theta],\[Phi]]]]^2*getNormalizationFactor[parity,m,n,p,a,TMorTE,k]^2,{r,0.0,a},{\[Theta],0.0,Pi},{\[Phi],0.0,2.0*Pi},Method->"AdaptiveMonteCarlo",MaxPoints->10^5,PrecisionGoal->2]







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

overlapIntegrand[mVec_,iVec_,jVec_,a_,modeCoefs_,EMamplitude_,TMorTEFeildi_,TMorTEFeildj_,parityOfFieldi_,parityOfFieldj_,ki_,kj_,wi_,wj_,\[Theta]_,\[Phi]_,offsets_]:=(a^2*Sin[\[Theta]] * (displacementVecQ[mVec[[1]],mVec[[2]],mVec[[3]],a,modeCoefs[[1]],modeCoefs[[2]],modeCoefs[[3]],modeCoefs[[4]],modeCoefs[[5]],modeCoefs[[6]],#1,#2,#3]&@@newPoint[{a,\[Theta],\[Phi]},{offsets[[1]],offsets[[2]],offsets[[3]]}].{1,0,0}) * ((wi/wj) * (If[TMorTEFeildj==1,magneticFieldTM[EMamplitude,parityOfFieldj, jVec[[1]], jVec[[2]], jVec[[3]],a,kj,#1,#2,#3]&@@newPoint[{a,\[Theta],\[Phi]},{offsets[[4]],offsets[[5]],offsets[[6]]}],magneticFieldTE[EMamplitude,parityOfFieldj, jVec[[1]], jVec[[2]], jVec[[3]],a,kj,#1,#2,#3]&@@newPoint[{a,\[Theta],\[Phi]},{offsets[[4]],offsets[[5]],offsets[[6]]}]].  Conjugate[If[TMorTEFeildi==1,magneticFieldTM[EMamplitude,parityOfFieldi, iVec[[1]], iVec[[2]], iVec[[3]],a,ki,#1,#2,#3]&@@newPoint[{a,\[Theta],\[Phi]},{offsets[[7]],offsets[[8]],offsets[[9]]}],magneticFieldTE[EMamplitude,parityOfFieldi, iVec[[1]], iVec[[2]], iVec[[3]],a,ki,#1,#2,#3]&@@newPoint[{a,\[Theta],\[Phi]},{offsets[[7]],offsets[[8]],offsets[[9]]}]]]) - (If[TMorTEFeildj==1,electricFieldTM[EMamplitude,parityOfFieldj, jVec[[1]], jVec[[2]], jVec[[3]],a,kj,#1,#2,#3]&@@newPoint[{a,\[Theta],\[Phi]},{offsets[[4]],offsets[[5]],offsets[[6]]}],electricFieldTE[EMamplitude,parityOfFieldj, jVec[[1]], jVec[[2]], jVec[[3]],a,kj,#1,#2,#3]&@@newPoint[{a,\[Theta],\[Phi]},{offsets[[4]],offsets[[5]],offsets[[6]]}]].Conjugate[If[TMorTEFeildi==1,electricFieldTM[EMamplitude,parityOfFieldi, iVec[[1]], iVec[[2]], iVec[[3]],a,ki,#1,#2,#3]&@@newPoint[{a,\[Theta],\[Phi]},{offsets[[7]],offsets[[8]],offsets[[9]]}],electricFieldTE[EMamplitude,parityOfFieldi, iVec[[1]], iVec[[2]], iVec[[3]],a,ki,#1,#2,#3]&@@newPoint[{a,\[Theta],\[Phi]},{offsets[[7]],offsets[[8]],offsets[[9]]}]]])))

overlapfunction[modes_,a_,TMorTEFeildi_,TMorTEFeildj_,parityOfFieldi_,parityOfFieldj_,ki_,kj_,wi_,wj_,offsets_]:=((4*Pi*a^3/3)^(1/3)/(2*calculateEMFieldEnergy[TMorTEFeildi,1.0,parityOfFieldi,modes[[4]],modes[[5]],modes[[6]],a,ki])) *NIntegrate[Re[Activate[overlapIntegrand[{modes[[1]],modes[[2]],modes[[3]]},{modes[[4]],modes[[5]],modes[[6]]},{modes[[7]],modes[[8]],modes[[9]]},a,{"C0","D0",0.0,0.0,"C2","D2"}/.Quiet[(StringJoin["m",ToString[modes[[1]]],"n",ToString[modes[[2]]],"l",ToString[modes[[3]]]]/.MODES[[1]])],1.0,TMorTEFeildi,TMorTEFeildj,parityOfFieldi,parityOfFieldj,ki,kj,wi,wj,\[Theta],\[Phi],offsets]]], {\[Theta],0.0,Pi},{\[Phi],0.0,2.0*Pi}, Method->{"AdaptiveMonteCarlo"},MaxPoints->10^4,PrecisionGoal->2]*getNormalizationFactor[parityOfFieldi,modes[[4]],modes[[5]],modes[[6]],a,TMorTEFeildi,ki]*getNormalizationFactor[parityOfFieldj,modes[[7]],modes[[8]],modes[[9]],a,TMorTEFeildj,kj]

overlap1Sphere2Modes[modes_,a_,TMorTEFeildi_,TMorTEFeildj_,parityOfFieldi_,parityOfFieldj_,offsets_]:=overlapfunction[modes,a,TMorTEFeildi,TMorTEFeildj,parityOfFieldi,parityOfFieldj,If[TMorTEFeildi==1,Activate[wavenumberTM[modes[[4]],modes[[5]],modes[[6]],a]],Activate[wavenumberTE[modes[[4]],modes[[5]],modes[[6]],a]]],If[TMorTEFeildi==1,Activate[wavenumberTM[modes[[7]],modes[[8]],modes[[9]],a]],Activate[wavenumberTE[modes[[7]],modes[[8]],modes[[9]],a]]],If[TMorTEFeildi==1,Activate[omegaTM[modes[[4]],modes[[5]],modes[[6]],a]],Activate[omegaTE[modes[[4]],modes[[5]],modes[[6]],a]]],If[TMorTEFeildi==1,Activate[omegaTM[modes[[7]],modes[[8]],modes[[9]],a]],Activate[omegaTE[modes[[7]],modes[[8]],modes[[9]],a]]],offsets]










(* GET FREQUENCY SPLITTINGS *)
couplingPoints={{innerRadius,Pi/2.0,0.0},{innerRadius,Pi/2.0,Pi}};

Fpotential[parity_, m_,n_,p_,a_,k_,r_,\[Theta]_,\[Phi]_]:=vectorPotentialTE[1.0,parity, m,n,p,a,k,r,\[Theta],\[Phi]]/getNormalizationFactor[parity,m,n,p,a,0,k]
Apotential[parity_, m_,n_,p_,a_,k_,r_,\[Theta]_,\[Phi]_]:=Curl[Fpotential[parity, m,n,p,a,k,r,\[Theta],\[Phi]],{r,\[Theta],\[Phi]},"Spherical"]/k

garbageVariableUsedToDefineFunctions=2*Norm[Fpotential[parity,m,n,p,a,k,#1,#2,#3]&@@newPoint[{couplingPoints[[1]][[1]],couplingPoints[[1]][[2]],couplingPoints[[1]][[3]]},{offsets[[4]],offsets[[5]],offsets[[6]]}]]^2-(Apotential[parity,m,n,p,a,k,#1,#2,#3][[1]]&@@newPoint[{couplingPoints[[1]][[1]],couplingPoints[[1]][[2]],couplingPoints[[1]][[3]]},{offsets[[4]],offsets[[5]],offsets[[6]]}])^2;
c\[Alpha]\[Alpha][parity_, m_,n_,p_,a_,k_,offsets_]:=Evaluate[garbageVariableUsedToDefineFunctions]
Clear[garbageVariableUsedToDefineFunctions];

garbageVariableUsedToDefineFunctions=2*Norm[Fpotential[parity,m,n,p,a,k,#1,#2,#3]&@@newPoint[{couplingPoints[[2]][[1]],couplingPoints[[2]][[2]],couplingPoints[[2]][[3]]},{offsets[[7]],offsets[[8]],offsets[[9]]}]]^2-(Apotential[parity,m,n,p,a,k,#1,#2,#3][[1]]&@@newPoint[{couplingPoints[[2]][[1]],couplingPoints[[2]][[2]],couplingPoints[[2]][[3]]},{offsets[[7]],offsets[[8]],offsets[[9]]}])^2;
c\[Beta]\[Beta][parity_, m_,n_,p_,a_,k_,offsets_]:=Evaluate[garbageVariableUsedToDefineFunctions]
Clear[garbageVariableUsedToDefineFunctions];

garbageVariableUsedToDefineFunctions=2*Fpotential[parity,m,n,p,a,k,#1,#2,#3]&@@newPoint[{couplingPoints[[1]][[1]],couplingPoints[[1]][[2]],couplingPoints[[1]][[3]]},{offsets[[4]],offsets[[5]],offsets[[6]]}].Fpotential[parity,m,n,p,a,k,#1,#2,#3]&@@newPoint[{couplingPoints[[2]][[1]],couplingPoints[[2]][[2]],couplingPoints[[2]][[3]]},{offsets[[7]],offsets[[8]],offsets[[9]]}]-(Apotential[parity,m,n,p,a,k,#1,#2,#3][[1]]&@@newPoint[{couplingPoints[[1]][[1]],couplingPoints[[1]][[2]],couplingPoints[[1]][[3]]},{offsets[[4]],offsets[[5]],offsets[[6]]}])*(Apotential[parity,m,n,p,a,k,#1,#2,#3][[1]]&@@newPoint[{couplingPoints[[2]][[1]],couplingPoints[[2]][[2]],couplingPoints[[2]][[3]]},{offsets[[7]],offsets[[8]],offsets[[9]]}]);
c\[Alpha]\[Beta][parity_, m_,n_,p_,a_,k_,offsets_]:=Evaluate[garbageVariableUsedToDefineFunctions]
Clear[garbageVariableUsedToDefineFunctions];

(* solve for the frequency shift from quadratic *)
garbageVariableUsedToDefineFunctions=(( 2 - (4/3)*(d^3/V)*(c\[Alpha]\[Alpha][parity, m,n,p,a,k,offsets]+c\[Beta]\[Beta][parity, m,n,p,a,k,offsets]))/2 + (1/2)*(( 2 - (4/3)*(d^3/V)*(c\[Alpha]\[Alpha][parity, m,n,p,a,k,offsets]+c\[Beta]\[Beta][parity, m,n,p,a,k,offsets]))^2 - 4*(1 - (4/3)*(d^3/V)*(c\[Alpha]\[Alpha][parity, m,n,p,a,k,offsets]+c\[Beta]\[Beta][parity, m,n,p,a,k,offsets]) + (16/9)*(d^6/V^2)*(c\[Alpha]\[Alpha][parity, m,n,p,a,k,offsets]*c\[Beta]\[Beta][parity, m,n,p,a,k,offsets]-c\[Alpha]\[Beta][parity, m,n,p,a,k,offsets]^2)))^(1/2));
frequencyCoefPlus[parity_, m_,n_,p_,a_,k_,d_,V_,offsets_]:=Evaluate[garbageVariableUsedToDefineFunctions]
Clear[garbageVariableUsedToDefineFunctions];

garbageVariableUsedToDefineFunctions=(( 2 - (4/3)*(d^3/V)*(c\[Alpha]\[Alpha][parity, m,n,p,a,k,offsets]+c\[Beta]\[Beta][parity, m,n,p,a,k,offsets]))/2 - (1/2)*(( 2 - (4/3)*(d^3/V)*(c\[Alpha]\[Alpha][parity, m,n,p,a,k,offsets]+c\[Beta]\[Beta][parity, m,n,p,a,k,offsets]))^2 - 4*(1 - (4/3)*(d^3/V)*(c\[Alpha]\[Alpha][parity, m,n,p,a,k,offsets]+c\[Beta]\[Beta][parity, m,n,p,a,k,offsets]) + (16/9)*(d^6/V^2)*(c\[Alpha]\[Alpha][parity, m,n,p,a,k,offsets]*c\[Beta]\[Beta][parity, m,n,p,a,k,offsets]-c\[Alpha]\[Beta][parity, m,n,p,a,k,offsets]^2)))^(1/2));
frequencyCoefMinus[parity_, m_,n_,p_,a_,k_,d_,V_,offsets_]:=Evaluate[garbageVariableUsedToDefineFunctions]
Clear[garbageVariableUsedToDefineFunctions];

getShiftedFrequencies[parity_,m_,n_,p_,a_,d_,tmorte_,offsets_]:=Activate[If[tmorte==1,omegaTM[m,n,p,a],omegaTE[m,n,p,a]]*{frequencyCoefPlus[parity, m,n,p,a,If[tmorte==1,wavenumberTM[m,n,p,a],wavenumberTE[m,n,p,a]],d,(4.0/3.0)*Pi*a^3,offsets],frequencyCoefMinus[parity, m,n,p,a,If[tmorte==1,wavenumberTM[m,n,p,a],wavenumberTE[m,n,p,a]],d,(4.0/3.0)*Pi*a^3,offsets] }]








(* WRITE DOWN OVERLAP FUNCTION FOR COUPLED SPHERES *)
overlapCoupledSpheres[modes_,a_,TMorTEFeildi_,TMorTEFeildj_,parityOfFieldi_,parityOfFieldj_,d_,offsets_]:=overlapfunction[{modes[[1]],modes[[2]],modes[[3]],modes[[4]],modes[[5]],modes[[6]],modes[[4]],modes[[5]],modes[[6]]},a,TMorTEFeildi,TMorTEFeildj,parityOfFieldi,parityOfFieldj,Evaluate[getShiftedFrequencies[parityOfFieldi,modes[[4]],modes[[5]],modes[[6]],a,d,TMorTEFeildi,offsets][[1]]]/SPEEDOFLIGHT,
Evaluate[getShiftedFrequencies[parityOfFieldi,modes[[4]],modes[[5]],modes[[6]],a,d,TMorTEFeildi,offsets][[2]]]/SPEEDOFLIGHT,Evaluate[getShiftedFrequencies[parityOfFieldi,modes[[4]],modes[[5]],modes[[6]],a,d,TMorTEFeildi,offsets][[1]]],
Evaluate[getShiftedFrequencies[parityOfFieldi,modes[[4]],modes[[5]],modes[[6]],a,d,TMorTEFeildi,offsets][[2]]],{offsets[[1]],offsets[[2]],offsets[[3]],offsets[[4]],offsets[[5]],offsets[[6]],offsets[[4]],offsets[[5]],offsets[[6]]}]+overlapfunction[{modes[[1]],modes[[2]],modes[[3]],modes[[4]],modes[[5]],modes[[6]],modes[[4]],modes[[5]],modes[[6]]},a,TMorTEFeildi,TMorTEFeildj,parityOfFieldi,parityOfFieldj,Evaluate[getShiftedFrequencies[parityOfFieldi,modes[[4]],modes[[5]],modes[[6]],a,d,TMorTEFeildi,offsets][[1]]]/SPEEDOFLIGHT,
Evaluate[getShiftedFrequencies[parityOfFieldi,modes[[4]],modes[[5]],modes[[6]],a,d,TMorTEFeildi,offsets][[2]]]/SPEEDOFLIGHT,Evaluate[getShiftedFrequencies[parityOfFieldi,modes[[4]],modes[[5]],modes[[6]],a,d,TMorTEFeildi,offsets][[1]]],
Evaluate[getShiftedFrequencies[parityOfFieldi,modes[[4]],modes[[5]],modes[[6]],a,d,TMorTEFeildi,offsets][[2]]],{offsets[[1]],offsets[[2]],offsets[[3]],offsets[[7]],offsets[[8]],offsets[[9]],offsets[[7]],offsets[[8]],offsets[[9]]}]







(* DO SOME BITCHIN' CALCULATIONS *)
Print[overlapCoupledSpheres[{2,0,0,0,1,2},innerRadius,0,0,1,1,0.1,{0,0,0,0,0,0,0,0,0}]//AbsoluteTiming]



(* make a list of all possible mode combinations *)
(*wholeDataArray = ParallelTable[{2,ToExpression[argv[[7]]],lm,mi,ni,pi,overlapCoupledSpheres[{2,ToExpression[argv[[7]]],lm,mi,ni,pi}, innerRadius, 0, 0, 1, 1]}, {lm,0,2,1}, {ni,0,2,1}, {mi,0,2,1}, {pi,1,3,1}];

Print[wholeDataArray]

Export[StringJoin[directory,"overlaps/overlaps-coupled-spheres-",ToString[argv[[7]]],".txt"], wholeDataArray];*)
