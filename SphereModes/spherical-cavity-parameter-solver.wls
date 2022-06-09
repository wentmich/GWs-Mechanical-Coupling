(*Limits Mathematica to requested resources*)
Unprotect[$ProcessorCount];$ProcessorCount = 1;

(*Prints the machine name that each kernel is running on*)
Print[ParallelEvaluate[$MachineName]];





(* IMPORT CONSTANTS FROM COMMAND LINE *)
(* argv should be "inner radius, thickness" *)
argv = Rest @ $CommandLine;
argc = Length @ argv;
Print[ToString[argc]]




(* DEFINE CONSTANTS *)
SPEEDOFLIGHT = 3.0*10.0^8.0;
shearModulus = 38.0*10^9;
youngsModulus = 105.0*10^9;
matieralDensity = 8.57*10^(-3)/((10^(-2))^3);

lambdaMaterial = -shearModulus*(youngsModulus - 2*shearModulus)/(youngsModulus - 3*shearModulus);
vt = (matieralDensity/shearModulus)^(-1/2);
vl = (matieralDensity/(lambdaMaterial + 2*shearModulus))^(-1/2);

innerRadius = ToExpression[argv[[1]]];
outerRadius = ToExpression[argv[[1]]] + ToExpression[argv[[2]]];
maxWaveNumber = 10.0^3.0;

Print[innerRadius]
Print[outerRadius]



(* DEFINE THE MODE DATA STRUCTURE *)
makeModeVariable[m_, n_, l_, \[Omega]_, k_, p_, C0_, C2_, D0_, D2_] := StringJoin["m", ToString[m], "n", ToString[n], "l", ToString[l]] -> {"wavenumberk" -> k, "wavenumberp" -> p, "frequency" -> \[Omega], "C0" -> C0, "C2" -> C2, "D0" -> D0, "D2" -> D2}


Print["defined make mode variable"]


(* DEFINE FUNCTIONS TO CALCULATE RESONANT FREQUENCIES *)
smallMmatrix[j1y0_, r_, m_, n_, a_, pmn_, kmn_] := {{(m/(pmn*r) - lambdaMaterial/(2*shearModulus))*If[j1y0 == 1, SphericalBesselJ[m, pmn*r], SphericalBesselY[m, pmn*r]] - If[j1y0 == 1, SphericalBesselJ[m + 1, pmn*r], SphericalBesselY[m + 1, pmn*r]], -((m*(m + 1)*((m - 1)*If[j1y0 == 1, SphericalBesselJ[m, kmn*r], SphericalBesselY[m, kmn*r]] - kmn*r*If[j1y0 == 1, SphericalBesselJ[m + 1, kmn*r], SphericalBesselY[m + 1, kmn*r]]))/(pmn^2*r^2))}, {((m - 1)*If[j1y0 == 1, SphericalBesselJ[m, pmn*r], SphericalBesselY[m, pmn*r]] - pmn*r*If[j1y0 == 1, SphericalBesselJ[m + 1, pmn*r], SphericalBesselY[m + 1, pmn*r]])/(pmn^2*r^2), ((kmn^2 * r^2*If[j1y0 == 1, SphericalBesselJ[m + 1, kmn*r], SphericalBesselY[m + 1, kmn*r]] - (kmn*m*r + m^2 + m - 2)*If[j1y0 == 1, SphericalBesselJ[m, kmn*r], SphericalBesselY[m, kmn*r]])/(2*pmn^2*r^2))}};

bigMmatrix[m_,n_,a_,outerShell_,innerShell_,pmn_,kmn_]:=Evaluate[ArrayFlatten[{{smallMmatrix[1,outerShell,m,n,a,pmn,kmn],smallMmatrix[0,outerShell,m,n,a,pmn,kmn]},{smallMmatrix[1,innerShell,m,n,a,pmn,kmn],smallMmatrix[0,innerShell,m,n,a,pmn,kmn]}}]]

myDeterminant=Simplify[Det[bigMmatrix[m,1,a,b,a,pmn,kmn]]];
matrixMdeterminent[m_,a_,b_,pmn_,kmn_]:=Evaluate[myDeterminant]
Clear[myDeterminant];

Print["right before solving"]
waveNumbersK=Table[kmn /. Quiet[Solve[{Activate[matrixMdeterminent[m,innerRadius,outerRadius,kmn*Sqrt[shearModulus/(lambdaMaterial+2*shearModulus)],kmn]]==0,0<=kmn<=maxWaveNumber},kmn]],{m,0,2}]//N;
waveNumbersP=Sqrt[shearModulus/(lambdaMaterial+2*shearModulus)]*waveNumbersK;
resonantFrequencies=vt*Sqrt[waveNumbersK];


Print["defined functions to calculate frequencies"]


(* CHECKS OF RESONANT WAVE NUMBERS *)
matrixDeterminantValues=Table[{m,kmn,Log10[Abs[matrixMdeterminent[m,innerRadius,outerRadius,kmn*Sqrt[shearModulus/(lambdaMaterial+2*shearModulus)],kmn]]]},{kmn,0.1,10.0,0.01},{m,0.0,2.0,0.1}]//AbsoluteTiming;
matrixDeterminantValues2=ArrayFlatten[matrixDeterminantValues[[2]],1];
lowKplot = ListDensityPlot[matrixDeterminantValues2, PlotTheme -> "Scientific", PlotLabel -> "det M_sph",FrameLabel -> {"m", "det(|M|)"}, ColorFunction -> "BlueGreenYellow", PlotLegends -> Automatic, InterpolationOrder -> 4];

Print["made 1st plot"]

matrixDeterminantValues=Table[{m,kmn,Log10[Abs[matrixMdeterminent[m,innerRadius,outerRadius,kmn*Sqrt[shearModulus/(lambdaMaterial+2*shearModulus)],kmn]]]},{kmn,1.0,1000.0,1.0},{m,0.0,2.0,0.1}]//AbsoluteTiming
matrixDeterminantValues2=ArrayFlatten[matrixDeterminantValues[[2]],1];
highKplot = ListDensityPlot[matrixDeterminantValues2, PlotTheme -> "Scientific", PlotLabel -> "det M_sph",FrameLabel -> {"m", "det(|M|)"}, ColorFunction -> "BlueGreenYellow", PlotLegends -> Automatic, InterpolationOrder -> 4];

Print["made 2nd plot"]

logDeterminantValuesM1=Table[{kmn,Log10[Abs[matrixMdeterminent[1,innerRadius,outerRadius,kmn*Sqrt[shearModulus/(lambdaMaterial+2*shearModulus)],kmn]]]},{kmn,0.01,3.0,0.01}];
logDeterminantValuesM2=Table[{kmn,Log10[Abs[matrixMdeterminent[2,innerRadius,outerRadius,kmn*Sqrt[shearModulus/(lambdaMaterial+2*shearModulus)],kmn]]]},{kmn,0.01,3.0,0.01}];
zerosPlot =  ListPlot[{logDeterminantValuesM1, logDeterminantValuesM2}, PlotTheme -> "Scientific", FrameLabel -> {"m", "Log(det(|M|))"}, Joined -> True];

Print["made 3rd plot"]

checkGrid = GraphicsGrid[{lowKplot, highKplot, zerosPlot}];
Print["made grid"]
Export[StringJoin[argv[[3]], "plots-of-mechanical-mode-frequencies.png"], checkGrid, "PNG"];


Print["checked resonant frequencies"]


(* WRITE DOWN MECHANICAL MODES *)
(* mechanical frequencies *)
vt = (matieralDensity/shearModulus)^(-1/2);
vl = (matieralDensity/(lambdaMaterial+2*shearModulus))^(-1/2);
(* wave numbers *)
wavek[m_,n_,l_,a_]:=Inactivate[waveNumbersK[[m+1,n+1]]];
wavep[m_,n_,l_,a_]:=Inactivate[waveNumbersP[[m+1,n+1]]];
omegaMech[m_,n_,l_,a_]:=vt*Sqrt[wavek[m,n,l,a]];
(* constituent functions *)
aux\[Phi][m_,n_,l_,a_,r_,\[Theta]_,\[Phi]_]:=SphericalBesselJ[m,wavep[m,n,l,a]*r]*SphericalHarmonicY[m,l,\[Theta],\[Phi]];
aux\[Phi]T[m_,n_,l_,a_,r_,\[Theta]_,\[Phi]_]:=SphericalBesselY[m,wavep[m,n,l,a]*r]*SphericalHarmonicY[m,l,\[Theta],\[Phi]];
aux\[Psi][m_,n_,l_,a_,r_,\[Theta]_,\[Phi]_]:=SphericalBesselJ[m,wavek[m,n,l,a]*r]*SphericalHarmonicY[m,l,\[Theta],\[Phi]];
aux\[Psi]T[m_,n_,l_,a_,r_,\[Theta]_,\[Phi]_]:=SphericalBesselY[m,wavek[m,n,l,a]*r]*SphericalHarmonicY[m,l,\[Theta],\[Phi]];
(* mechanical mode displacement vector *)
garbageVariableUsedToDefineFunctions=c0*(Grad[aux\[Phi][m,n,l,a,r,\[Theta],\[Phi]],{r,\[Theta],\[Phi]},"Spherical"] + (d0/c0)*Grad[aux\[Phi]T[m,n,l,a,r,\[Theta],\[Phi]],{r,\[Theta],\[Phi]},"Spherical"]+I*(c2/c0)*Curl[Cross[{r,\[Theta],\[Phi]},Grad[aux\[Psi][m,n,l,a,r,\[Theta],\[Phi]],{r,\[Theta],\[Phi]},"Spherical"]],{r,\[Theta],\[Phi]},"Spherical"]+I*(d2/c0)*Curl[Cross[{r,\[Theta],\[Phi]},Grad[aux\[Psi]T[m,n,l,a,r,\[Theta],\[Phi]],{r,\[Theta],\[Phi]},"Spherical"]],{r,\[Theta],\[Phi]},"Spherical"]);
displacementVecQ[m_,n_,l_,a_,c0_,d0_,c2_,d2_,r_,\[Theta]_,\[Phi]_]:=Evaluate[garbageVariableUsedToDefineFunctions];
Clear[garbageVariableUsedToDefineFunctions];
garbageVariableUsedToDefineFunctions=(Grad[aux\[Phi][m,n,l,a,r,\[Theta],\[Phi]],{r,\[Theta],\[Phi]},"Spherical"] + (d0ratio)*Grad[aux\[Phi]T[m,n,l,a,r,\[Theta],\[Phi]],{r,\[Theta],\[Phi]},"Spherical"]+I*(c2ratio)*Curl[Cross[{r,\[Theta],\[Phi]},Grad[aux\[Psi][m,n,l,a,r,\[Theta],\[Phi]],{r,\[Theta],\[Phi]},"Spherical"]],{r,\[Theta],\[Phi]},"Spherical"]+I*(d2ratio)*Curl[Cross[{r,\[Theta],\[Phi]},Grad[aux\[Psi]T[m,n,l,a,r,\[Theta],\[Phi]],{r,\[Theta],\[Phi]},"Spherical"]],{r,\[Theta],\[Phi]},"Spherical"]);
displacementVecQnorm[m_,n_,l_,a_,d0ratio_,c2ratio_,d2ratio_,r_,\[Theta]_,\[Phi]_]:=Evaluate[garbageVariableUsedToDefineFunctions];
Clear[garbageVariableUsedToDefineFunctions];





(* SOLVE FOR FREQUENCIES AND COEFFICIENTS AND EXPORT DATA *)
findCoefRatios[m_,n_]:=ArrayFlatten[LinearSolve[bigMmatrix[m,n,innerRadius,outerRadius,innerRadius,(waveNumbersK[[m+1]][[n+1]])*Sqrt[shearModulus/(lambdaMaterial+2*shearModulus)],waveNumbersK[[m+1]][[n+1]]][[1;;4,2;;4]],-bigMmatrix[m,n,innerRadius,outerRadius,innerRadius,(waveNumbersK[[m+1]][[n+1]])*Sqrt[shearModulus/(lambdaMaterial+2*shearModulus)],waveNumbersK[[m+1]][[n+1]]][[1;;4,1;;1]]],1]
findNormalization[m_,n_,l_]:=(NIntegrate[r^2*Sin[\[Theta]]*Norm[displacementVecQnorm[m,n,l,innerRadius,findCoefRatios[[2]],findCoefRatios[[1]],findCoefRatios[[3]],r,\[Theta],\[Phi]]]^2,{r,0.0,innerRadius},{\[Theta],0.0,Pi},{\[Phi],0.0,2.0*Pi}])^(-1/2)
modeInfoDataStructure={};
For[m=2,m<3,m++,For[n=0,n<6,n++,For[l=0,l<=n,l++,ratios=findCoefRatios[m,n];c0=Quiet[(Quiet[NIntegrate[Activate[r^2*Sin[\[Theta]]*Quiet[Norm[Quiet[displacementVecQnorm[m,n,l,innerRadius,ratios[[2]],ratios[[1]],ratios[[3]],r,\[Theta],\[Phi]]]]]^2],{r,innerRadius,outerRadius},{\[Theta],0.01,Pi},{\[Phi],0.01,2.0*Pi},Method->"AdaptiveMonteCarlo",MaxPoints->10^5]]/((4/3)*Pi*outerRadius^3 - (4/3)*Pi*innerRadius^3))^(-1/2)];modeInfoDataStructure=Append[modeInfoDataStructure,makeModeVariable[m,n,l,resonantFrequencies[[m+1,n+1]],waveNumbersK[[m+1,n+1]],waveNumbersP[[m+1,n+1]],c0,ratios[[2]]*c0,ratios[[1]]*c0,ratios[[3]]*c0]]]]]

Quiet[Export[StringJoin[argv[[3]],"NbCavity-r=1.0-1.01.txt"], InputForm[modeInfoDataStructure]]];
