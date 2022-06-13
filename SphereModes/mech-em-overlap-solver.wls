(*Limits Mathematica to requested resources*)
Unprotect[$ProcessorCount];$ProcessorCount = 40;

(*Prints the machine name that each kernel is running on*)
Print[ParallelEvaluate[$MachineName]];


argv = Rest @ $CommandLine;
argc = Length @ argv;


(* define all physical constants (not cavity parameters) *)
SPEEDOFLIGHT=3.0*10.0^8.0;
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
directory=argv[[5]];
filename=argv[[6]];
maxWaveNumber=10.0^3;


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
garbageVariableUsedToDefineFunctions=c0*Grad[aux\[Phi][m,n,l,a,r,\[Theta],\[Phi]],{r,\[Theta],\[Phi]},"Spherical"] + d0*Grad[aux\[Phi]T[m,n,l,a,r,\[Theta],\[Phi]],{r,\[Theta],\[Phi]},"Spherical"]+c1*Cross[{r,\[Theta],\[Phi]},Grad[aux\[Psi][m,n,l,a,r,\[Theta],\[Phi]],{r,\[Theta],\[Phi]},"Spherical"]]+d1*Cross[{r,\[Theta],\[Phi]},Grad[aux\[Psi]T[m,n,l,a,r,\[Theta],\[Phi]],{r,\[Theta],\[Phi]},"Spherical"]]+I*c2*Curl[Cross[{r,\[Theta],\[Phi]},Grad[aux\[Psi][m,n,l,a,r,\[Theta],\[Phi]],{r,\[Theta],\[Phi]},"Spherical"]],{r,\[Theta],\[Phi]},"Spherical"]+I*d2*Curl[Cross[{r,\[Theta],\[Phi]},Grad[aux\[Psi]T[m,n,l,a,r,\[Theta],\[Phi]],{r,\[Theta],\[Phi]},"Spherical"]],{r,\[Theta],\[Phi]},"Spherical"];
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
vectorPotentialTE[amplitude_,parity_, m_,n_,p_,a_,k_,r_,\[Theta]_,\[Phi]_]:={amplitude*r*SphericalBesselJ[n,Inactivate[SphericalBesselJZero[n,p]]*r/a]*LegendreP[n,m,Cos[\[Theta]]]*If[parity==1,Cos[m*\[Phi]],Sin[m*\[Phi]]],0,0}
(* Transverse magnetic vector potential *)
vectorPotentialTM[amplitude_,parity_,m_,n_,p_,a_,k_,r_,\[Theta]_,\[Phi]_]:={amplitude*r*SphericalBesselJ[n,Inactivate[SphericalBesselJDerivZero[n,p]]*r/a]*LegendreP[n,m,Cos[\[Theta]]]*If[parity==1,Cos[m*\[Phi]],Sin[m*\[Phi]]],0,0}
(* calculate electric and magnetic fields for TE modes *)
garbageVariableUsedToDefineFunctions=Evaluate[-Curl[vectorPotentialTE[amplitude,parity, m,n,p,a,k,r,\[Theta],\[Phi]],{r,\[Theta],\[Phi]},"Spherical"]];
electricFieldTE[amplitude_,parity_, m_,n_,p_,a_,k_,r_,\[Theta]_,\[Phi]_]:=Evaluate[garbageVariableUsedToDefineFunctions]
Clear[garbageVariableUsedToDefineFunctions]
garbageVariableUsedToDefineFunctions=Module[{r0,\[Theta]0,\[Phi]0},Evaluate[(1/(wavenumberTE[m,n,p,a]))*Curl[Curl[vectorPotentialTE[amplitude,parity, m,n,p,a,k,r0,\[Theta]0,\[Phi]0],{r0,\[Theta]0,\[Phi]0},"Spherical"],{r0,\[Theta]0,\[Phi]0},"Spherical"]]/.{r0->r,\[Theta]0->\[Theta],\[Phi]0->\[Phi]}];
magneticFieldTE[amplitude_,parity_, m_,n_,p_,a_,k_,r_,\[Theta]_,\[Phi]_]:=Evaluate[garbageVariableUsedToDefineFunctions]
Clear[garbageVariableUsedToDefineFunctions]
(* calculate electric and magnetic fields for TE modes *)
garbageVariableUsedToDefineFunctions=Module[{r0,\[Theta]0,\[Phi]0},Evaluate[(1/(wavenumberTM[m,n,p,a]))*Curl[Curl[vectorPotentialTM[amplitude,parity, m,n,p,a,k,r,\[Theta],\[Phi]],{r,\[Theta],\[Phi]},"Spherical"],{r,\[Theta],\[Phi]},"Spherical"]]/.{r0->r,\[Theta]0->\[Theta],\[Phi]0->\[Phi]}];
electricFieldTM[amplitude_,parity_, m_,n_,p_,a_,k_,r_,\[Theta]_,\[Phi]_]:=Evaluate[garbageVariableUsedToDefineFunctions]
Clear[garbageVariableUsedToDefineFunctions]
garbageVariableUsedToDefineFunctions=Module[{r0,\[Theta]0,\[Phi]0},Evaluate[Curl[vectorPotentialTM[amplitude,parity, m,n,p,a,k,r,\[Theta],\[Phi]],{r,\[Theta],\[Phi]},"Spherical"]]/.{r0->r,\[Theta]0->\[Theta],\[Phi]0->\[Phi]}];
magneticFieldTM[amplitude_,parity_, m_,n_,p_,a_,k_,r_,\[Theta]_,\[Phi]_]:=Evaluate[garbageVariableUsedToDefineFunctions]
Clear[garbageVariableUsedToDefineFunctions]



(* WRITE DOWN POTENTIAL ENERGY OF EM MODES *)
calculateEMFieldEnergy[TMorTE_,amplitude_,parity_,m_,n_,p_,a_]:=(1/2)*NIntegrate[r^2*Sin[\[Theta]]*Norm[If[TMorTE==0,electricFieldTE[amplitude,parity,m,n,p,a,If[TMorTE==0,wavenumberTE[m,n,p,a],wavenumberTM[m,n,p,a]],r,\[Theta],\[Phi]],electricFieldTM[amplitude,parity,m,n,p,a,If[TMorTE==0,wavenumberTE[m,n,p,a],wavenumberTM[m,n,p,a]],r,\[Theta],\[Phi]]]]^2,{r,0.0,a},{\[Theta],0.0,Pi},{\[Phi],0.0,2.0*Pi},AccuracyGoal->10]



(* WRITE DOWN THE OVERLAP INTEGRAND *)
overlapIntegrand[mVec_,iVec_,jVec_,a_,modeCoefs_,EMamplitude_,TMorTEFeildi_,TMorTEFeildj_,parityOfFieldi_,parityOfFieldj_,\[Theta]_,\[Phi]_]:=((4*Pi*a^3/3)^(1/3)/(2*calculateEMFieldEnergy[TMorTEFeildi,EMamplitude,parityOfFieldi,iVec[[1]],iVec[[2]],iVec[[3]],a])) *(a^2*Sin[\[Theta]] * (displacementVecQ[mVec[[1]],mVec[[2]],mVec[[3]],a,modeCoefs[[1]],modeCoefs[[2]],modeCoefs[[3]],modeCoefs[[4]],modeCoefs[[5]],modeCoefs[[6]],a,\[Theta],\[Phi]].{1,0,0}) * ((If[TMorTEFeildi==1,omegaTM[iVec[[1]],iVec[[2]],iVec[[3]],a],omegaTE[iVec[[1]],iVec[[2]],iVec[[3]],a]]/If[TMorTEFeildj==1,omegaTM[jVec[[1]],jVec[[2]],jVec[[3]],a],omegaTE[jVec[[1]],jVec[[2]],jVec[[3]],a]]) * (If[TMorTEFeildj==1,magneticFieldTM[EMamplitude,parityOfFieldj, jVec[[1]], jVec[[2]], jVec[[3]],a,wavenumberTM[jVec[[1]], jVec[[2]], jVec[[3]],a],a,\[Theta],\[Phi]],magneticFieldTE[EMamplitude,parityOfFieldj, jVec[[1]], jVec[[2]], jVec[[3]],a,wavenumberTE[jVec[[1]], jVec[[2]], jVec[[3]],a],a,\[Theta],\[Phi]]].  Conjugate[If[TMorTEFeildi==1,magneticFieldTM[EMamplitude,parityOfFieldi, iVec[[1]], iVec[[2]], iVec[[3]],a,wavenumberTM[iVec[[1]], iVec[[2]], iVec[[3]],a],a,\[Theta],\[Phi]],magneticFieldTE[EMamplitude,parityOfFieldi, iVec[[1]], iVec[[2]], iVec[[3]],a,wavenumberTE[iVec[[1]], iVec[[2]], iVec[[3]],a],a,\[Theta],\[Phi]]]]) - (If[TMorTEFeildj==1,electricFieldTM[EMamplitude,parityOfFieldj, jVec[[1]], jVec[[2]], jVec[[3]],a,wavenumberTM[jVec[[1]], jVec[[2]], jVec[[3]],a],a,\[Theta],\[Phi]],electricFieldTE[EMamplitude,parityOfFieldj, jVec[[1]], jVec[[2]], jVec[[3]],a,wavenumberTE[jVec[[1]], jVec[[2]], jVec[[3]],a],a,\[Theta],\[Phi]]].Conjugate[If[TMorTEFeildi==1,electricFieldTM[EMamplitude,parityOfFieldi, iVec[[1]], iVec[[2]], iVec[[3]],a,wavenumberTM[iVec[[1]], iVec[[2]], iVec[[3]],a],a,\[Theta],\[Phi]],electricFieldTE[EMamplitude,parityOfFieldi, iVec[[1]], iVec[[2]], iVec[[3]],a,wavenumberTE[iVec[[1]], iVec[[2]], iVec[[3]],a],a,\[Theta],\[Phi]]]])))



(* WRITE DOWN OVERLAP FUNCTION *)
overlapfunction[modes_]:=NIntegrate[Re[Activate[overlapIntegrand[{modes[[1]],modes[[2]],modes[[3]]},{modes[[4]],modes[[5]],modes[[6]]},{modes[[7]],modes[[8]],modes[[9]]},1.0,{"C0","D0",0.0,0.0,"C2","D2"}/.Quiet[(StringJoin["m",ToString[modes[[1]]],"n",ToString[modes[[2]]],"l",ToString[modes[[3]]]]/.MODES[[1]])],1.0,0,0,1,1,\[Theta],\[Phi]]]], {\[Theta],0.0,Pi},{\[Phi],0.0,2.0*Pi}, Method->{"GlobalAdaptive"},PrecisionGoal->2]


Print[overlapfunction[{2, 0, 0, 0, 1, 2, 0, 1, 2}]//AbsoluteTiming]




(* make a list of all possible mode combinations *)
wholeDataArray = {};
ParallelDo[Append[wholeDataArray,{2,arvg[[4]],lm,mi,ni,pi,mj,nj,pj,overlapfunction[{2,nm,lm,mi,ni,pi,mj,nj,pj}]}], {lm,0,2,1}, {ni,0,3,1}, {mi,0,ni,1}, {pi,1,3,1}, {nj,0,3,1}, {mj,0,nj,1}, {pj,1,3,1}];

Export[StringJoin[directory,"overlaps-",ToString[argv[[4]]],".csv"], wholeDataArray, "CSV"];





(*mmech = 2;
nmech = 0;
lmech = 0;
c0mech=Quiet["C0"/.(StringJoin["m",ToString[mmech],"n",ToString[nmech],"l",ToString[lmech]]/.MODES[[1]])];
c1mech=0.0;
c2mech=Quiet["C2"/.(StringJoin["m",ToString[mmech],"n",ToString[nmech],"l",ToString[lmech]]/.MODES[[1]])];(*0.0;*)
d0mech=Quiet["D0"/.(StringJoin["m",ToString[mmech],"n",ToString[nmech],"l",ToString[lmech]]/.MODES[[1]])];
d1mech=0.0;
d2mech=Quiet["D2"/.(StringJoin["m",ToString[mmech],"n",ToString[nmech],"l",ToString[lmech]]/.MODES[[1]])];(*0.0;*)

Print[NIntegrate[Re[Activate[overlapIntegrand[{mmech,nmech,lmech},{0,1,2},{0,1,2},innerRadius,{c0mech,d0mech,c1mech,d1mech,c2mech,d2mech},1.0,0,0,1,1,\[Theta],\[Phi]]]], {\[Theta],0.001,Pi},{\[Phi],0.0,2.0*Pi}, Method->{"GlobalAdaptive"},PrecisionGoal->2]]*)
