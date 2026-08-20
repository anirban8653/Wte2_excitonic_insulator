(* ::Package:: *)

SetDirectory["~/Desktop/exciton/bulk"];
(*Fourier Transformation from real to k-space*)
i = 0;
Nx=12;
Ny=12;
iter=10;
kernum = 48;

u0 = -0.64;
Nband=8;
flagx=1;
flagy=1;
Te = 1.65; (*onsite telurium potential for semimetal*)
For[iy = 1, iy <= Ny, iy++,
For[ix = 1, ix <= Nx, ix++, i++;
    xpos[i] = ix;
    ypos[i] = iy;
    ]];
      NN=i;
gso[rxs_,rys_]:={
{0.74 DiracDelta[rxs] DiracDelta[rys]-(0.41+0.008 I) DiracDelta[-1.+rxs,rys]-(0.41-0.008 I) DiracDelta[1.+rxs,rys],0.14 DiracDelta[1.-rxs,0.39+rys]-(0.+0.012 I) DiracDelta[rxs,0.61-rys]-(0.+0.012 I) DiracDelta[rxs,0.39+rys]-0.14 DiracDelta[1.+rxs,0.39+rys],0.51 DiracDelta[0.5-rxs,0.36-rys]+0.51 DiracDelta[0.5+rxs,0.36-rys],0.39 DiracDelta[0.5-rxs,0.25+rys]+0.29 DiracDelta[1.5-rxs,0.25+rys]-0.39 DiracDelta[0.5+rxs,0.25+rys]-0.29 DiracDelta[1.5+rxs,0.25+rys],-0.031 DiracDelta[1.-rxs,rys]+0.031 DiracDelta[1.+rxs,rys],-0.05 DiracDelta[rxs,0.61-rys]-0.051 DiracDelta[rxs,0.39+rys],0,-0.011 DiracDelta[0.5-rxs,0.25+rys]-0.011 DiracDelta[0.5+rxs,0.25+rys]},
{-0.14 DiracDelta[1.-rxs,0.39-rys]+(0.+0.012 I) DiracDelta[rxs,0.39-rys]+(0.+0.012 I) DiracDelta[rxs,0.61+rys]+0.14 DiracDelta[1.+rxs,0.39-rys],-Te DiracDelta[rxs] DiracDelta[rys]+(1.13-0.01 I) DiracDelta[-1.+rxs,rys]+0.13 DiracDelta[rxs,1.-rys]+0.13 DiracDelta[rxs,1.+rys]+(1.13+0.01 I) DiracDelta[1.+rxs,rys],-0.39 DiracDelta[0.5-rxs,0.25+rys]-0.29 DiracDelta[1.5-rxs,0.25+rys]+0.39 DiracDelta[0.5+rxs,0.25+rys]+0.29 DiracDelta[1.5+rxs,0.25+rys],0.4 DiracDelta[0.5-rxs,0.14-rys]+0.4 DiracDelta[0.5+rxs,0.14-rys],0.051 DiracDelta[rxs,0.39-rys]+0.05 DiracDelta[rxs,0.61+rys],-0.04 DiracDelta[1.-rxs,rys]+0.04 DiracDelta[1.+rxs,rys],-0.011 DiracDelta[0.5-rxs,0.25+rys]-0.011 DiracDelta[0.5+rxs,0.25+rys],0},
{0.51 DiracDelta[0.5-rxs,0.36+rys]+0.51 DiracDelta[0.5+rxs,0.36+rys],0.39 DiracDelta[0.5-rxs,0.25-rys]+0.29 DiracDelta[1.5-rxs,0.25-rys]-0.39 DiracDelta[0.5+rxs,0.25-rys]-0.29 DiracDelta[1.5+rxs,0.25-rys],0.74 DiracDelta[rxs] DiracDelta[rys]-(0.41-0.008 I) DiracDelta[-1.+rxs,rys]-(0.41+0.008 I) DiracDelta[1.+rxs,rys],0.14 DiracDelta[1.-rxs,0.39-rys]+(0.+0.012 I) DiracDelta[rxs,0.39-rys]+(0.+0.012 I) DiracDelta[rxs,0.61+rys]-0.14 DiracDelta[1.+rxs,0.39-rys],0,0.011 DiracDelta[0.5-rxs,0.25-rys]+0.011 DiracDelta[0.5+rxs,0.25-rys],0.031 DiracDelta[1.-rxs,rys]-0.031 DiracDelta[1.+rxs,rys],0.051 DiracDelta[rxs,0.39-rys]+0.05 DiracDelta[rxs,0.61+rys]},
{-0.39 DiracDelta[0.5-rxs,0.25-rys]-0.29 DiracDelta[1.5-rxs,0.25-rys]+0.39 DiracDelta[0.5+rxs,0.25-rys]+0.29 DiracDelta[1.5+rxs,0.25-rys],0.4 DiracDelta[0.5-rxs,0.14+rys]+0.4 DiracDelta[0.5+rxs,0.14+rys],-0.14 DiracDelta[1.-rxs,0.39+rys]-(0.+0.012 I) DiracDelta[rxs,0.61-rys]-(0.+0.012 I) DiracDelta[rxs,0.39+rys]+0.14 DiracDelta[1.+rxs,0.39+rys],-Te DiracDelta[rxs] DiracDelta[rys]+(1.13+0.01 I) DiracDelta[-1.+rxs,rys]+0.13 DiracDelta[rxs,1.-rys]+0.13 DiracDelta[rxs,1.+rys]+(1.13-0.01 I) DiracDelta[1.+rxs,rys],0.011 DiracDelta[0.5-rxs,0.25-rys]+0.011 DiracDelta[0.5+rxs,0.25-rys],0,-0.05 DiracDelta[rxs,0.61-rys]-0.051 DiracDelta[rxs,0.39+rys],0.04 DiracDelta[1.-rxs,rys]-0.04 DiracDelta[1.+rxs,rys]},
{0.031 DiracDelta[1.-rxs,rys]-0.031 DiracDelta[1.+rxs,rys],0.05 DiracDelta[rxs,0.61-rys]+0.051 DiracDelta[rxs,0.39+rys],0,0.011 DiracDelta[0.5-rxs,0.25+rys]+0.011 DiracDelta[0.5+rxs,0.25+rys],0.74 DiracDelta[rxs] DiracDelta[rys]-(0.41-0.008 I) DiracDelta[-1.+rxs,rys]-(0.41+0.008 I) DiracDelta[1.+rxs,rys],0.14 DiracDelta[1.-rxs,0.39+rys]+(0.+0.012 I) DiracDelta[rxs,0.61-rys]+(0.+0.012 I) DiracDelta[rxs,0.39+rys]-0.14 DiracDelta[1.+rxs,0.39+rys],0.51 DiracDelta[0.5-rxs,0.36-rys]+0.51 DiracDelta[0.5+rxs,0.36-rys],0.39 DiracDelta[0.5-rxs,0.25+rys]+0.29 DiracDelta[1.5-rxs,0.25+rys]-0.39 DiracDelta[0.5+rxs,0.25+rys]-0.29 DiracDelta[1.5+rxs,0.25+rys]},
{-0.051 DiracDelta[rxs,0.39-rys]-0.05 DiracDelta[rxs,0.61+rys],0.04 DiracDelta[1.-rxs,rys]-0.04 DiracDelta[1.+rxs,rys],0.011 DiracDelta[0.5-rxs,0.25+rys]+0.011 DiracDelta[0.5+rxs,0.25+rys],0,-0.14 DiracDelta[1.-rxs,0.39-rys]-(0.+0.012 I) DiracDelta[rxs,0.39-rys]-(0.+0.012 I) DiracDelta[rxs,0.61+rys]+0.14 DiracDelta[1.+rxs,0.39-rys],-Te DiracDelta[rxs] DiracDelta[rys]+(1.13+0.01 I) DiracDelta[-1.+rxs,rys]+0.13 DiracDelta[rxs,1.-rys]+0.13 DiracDelta[rxs,1.+rys]+(1.13-0.01 I) DiracDelta[1.+rxs,rys],-0.39 DiracDelta[0.5-rxs,0.25+rys]-0.29 DiracDelta[1.5-rxs,0.25+rys]+0.39 DiracDelta[0.5+rxs,0.25+rys]+0.29 DiracDelta[1.5+rxs,0.25+rys],0.4 DiracDelta[0.5-rxs,0.14-rys]+0.4 DiracDelta[0.5+rxs,0.14-rys]},
{0,-0.011 DiracDelta[0.5-rxs,0.25-rys]-0.011 DiracDelta[0.5+rxs,0.25-rys],-0.031 DiracDelta[1.-rxs,rys]+0.031 DiracDelta[1.+rxs,rys],-0.051 DiracDelta[rxs,0.39-rys]-0.05 DiracDelta[rxs,0.61+rys],0.51 DiracDelta[0.5-rxs,0.36+rys]+0.51 DiracDelta[0.5+rxs,0.36+rys],0.39 DiracDelta[0.5-rxs,0.25-rys]+0.29 DiracDelta[1.5-rxs,0.25-rys]-0.39 DiracDelta[0.5+rxs,0.25-rys]-0.29 DiracDelta[1.5+rxs,0.25-rys],0.74 DiracDelta[rxs] DiracDelta[rys]-(0.41+0.008 I) DiracDelta[-1.+rxs,rys]-(0.41-0.008 I) DiracDelta[1.+rxs,rys],0.14 DiracDelta[1.-rxs,0.39-rys]-(0.+0.012 I) DiracDelta[rxs,0.39-rys]-(0.+0.012 I) DiracDelta[rxs,0.61+rys]-0.14 DiracDelta[1.+rxs,0.39-rys]},
{-0.011 DiracDelta[0.5-rxs,0.25-rys]-0.011 DiracDelta[0.5+rxs,0.25-rys],0,0.05 DiracDelta[rxs,0.61-rys]+0.051 DiracDelta[rxs,0.39+rys],-0.04 DiracDelta[1.-rxs,rys]+0.04 DiracDelta[1.+rxs,rys],-0.39 DiracDelta[0.5-rxs,0.25-rys]-0.29 DiracDelta[1.5-rxs,0.25-rys]+0.39 DiracDelta[0.5+rxs,0.25-rys]+0.29 DiracDelta[1.5+rxs,0.25-rys],0.4 DiracDelta[0.5-rxs,0.14+rys]+0.4 DiracDelta[0.5+rxs,0.14+rys],-0.14 DiracDelta[1.-rxs,0.39+rys]+(0.+0.012 I) DiracDelta[rxs,0.61-rys]+(0.+0.012 I) DiracDelta[rxs,0.39+rys]+0.14 DiracDelta[1.+rxs,0.39+rys],-Te DiracDelta[rxs] DiracDelta[rys]+(1.13-0.01 I) DiracDelta[-1.+rxs,rys]+0.13 DiracDelta[rxs,1.-rys]+0.13 DiracDelta[rxs,1.+rys]+(1.13+0.01 I) DiracDelta[1.+rxs,rys]}};
(*chempot=SparseArray[{x_,x_}->-mu,{NN*Nband,NN*Nband}];*)
r[1]={0.,0.};
r[2]={0.,-0.39};
r[3]={0.5,-0.64};
r[4]={0.5,-0.25};
r[5]=r[1];r[6]=r[2];r[7]=r[3];r[8]=r[4];
For[m=1,m<=8,m++,
For[n=1,n<=8,n++,
For[cx=-2,cx<=2,cx++,
For[cy=-2,cy<=2,cy++,
rxs=cx+r[m][[1]]-r[n][[1]];
rys=cy-r[m][[2]]+r[n][[2]];
If[(cx==0&&cy==0&&m==n),t[-cx][-cy][m][n]=Chop[gso[rxs,rys][[m,n]]/DiracDelta[0]^2],
t[-cx][-cy][m][n]=Chop[gso[rxs,rys][[m,n]]/DiracDelta[0,0]]];]]]];
(*t[0][0][1][2]*)
(*(*******chemical potential matrix formation using sparse array*******)
icd=0;
For[iy=1,iy<=Ny,iy++,
For[ix=1,ix<=Nx,ix++,icd++;
xposcd[icd]=ix;
yposcd[icd]=iy;
]]*)
Clear[Chempot,mu]
Chempot[mu_]:= SparseArray[{x_, x_} -> mu, {NN*Nband, NN*Nband}](*Dcdw is the CDW order parameter*)
cdwMatrix[dco_]:= SparseArray[Band[{1,1}]->dco,{NN*Nband,NN*Nband}];
(*MatrixForm[Chempot[0,0.0]];*)
(*cdwMatrix[dco]//MatrixForm*)
For[lx=1,lx<=NN*Nband,lx++,
For[ly=1,ly<=NN*Nband,ly++,
If[lx<= NN*Nband/2,
m=Mod[lx-1/10,Nband/2]+1/10;i=Quotient[lx-0.1,Nband/2]+1;,
m=Mod[lx-1/10,Nband/2]+1/10+4;i=Quotient[lx-0.1,Nband/2]+1-NN;];
If[ly<= NN*Nband/2,
n=Mod[ly-1/10,Nband/2]+1/10;j=Quotient[ly-0.1,Nband/2]+1;,
n=Mod[ly-1/10,Nband/2]+1/10+4;j=Quotient[ly-0.1,Nband/2]+1-NN;];

(*For[lx=1,lx<=NN*Nband,lx++,
For[ly=1,ly<=NN*Nband,ly++,
m=Mod[lx-1/10,Nband]+1/10;i=Quotient[lx-0.1,Nband]+1;
n=Mod[ly-1/10,Nband]+1/10;j=Quotient[ly-0.1,Nband]+1;*)

tijmn[lx][ly]=If[NumericQ[t[xpos[i]-(xpos[j]-Nx)][ypos[i]-(ypos[j]-Ny)][m][n]],t[xpos[i]-(xpos[j]-Nx)][ypos[i]-(ypos[j]-Ny)][m][n]*flagx*flagy,
If[NumericQ[t[xpos[i]-(xpos[j]-Nx)][ypos[i]-(ypos[j])][m][n]],t[xpos[i]-(xpos[j]-Nx)][ypos[i]-(ypos[j])][m][n]*flagx,
If[NumericQ[t[xpos[i]-(xpos[j]-Nx)][ypos[i]-(ypos[j]+Ny)][m][n]],t[xpos[i]-(xpos[j]-Nx)][ypos[i]-(ypos[j]+Ny)][m][n]*flagx*flagy,
If[NumericQ[t[xpos[i]-(xpos[j])][ypos[i]-(ypos[j]-Ny)][m][n]],t[xpos[i]-(xpos[j])][ypos[i]-(ypos[j]-Ny)][m][n]*flagy,
If[NumericQ[t[xpos[i]-(xpos[j])][ypos[i]-(ypos[j])][m][n]],t[xpos[i]-(xpos[j])][ypos[i]-(ypos[j])][m][n],
If[NumericQ[t[xpos[i]-(xpos[j])][ypos[i]-(ypos[j]+Ny)][m][n]],t[xpos[i]-(xpos[j])][ypos[i]-(ypos[j]+Ny)][m][n]*flagy,
If[NumericQ[t[xpos[i]-(xpos[j]+Nx)][ypos[i]-(ypos[j]-Ny)][m][n]],t[xpos[i]-(xpos[j]+Nx)][ypos[i]-(ypos[j]-Ny)][m][n]*flagx*flagy,
If[NumericQ[t[xpos[i]-(xpos[j]+Nx)][ypos[i]-(ypos[j])][m][n]],t[xpos[i]-(xpos[j]+Nx)][ypos[i]-(ypos[j])][m][n]*flagx,
If[NumericQ[t[xpos[i]-(xpos[j]+Nx)][ypos[i]-(ypos[j]+Ny)][m][n]],t[xpos[i]-(xpos[j]+Nx)][ypos[i]-(ypos[j]+Ny)][m][n]]*flagx*flagy
]]]]]]]];
]];
Hop=SparseArray[{x_,y_}/;NumericQ[tijmn[x][y]]->tijmn[x][y],{NN*Nband,NN*Nband}];
kT=0.001;
f[x_]:=If[Re[x]>0.02,0.0,
If[Re[x]<-0.02,1.0,
1.0/(Exp[Re[x]/kT]+1.0)]];
(*initialization*)

For[ mu = -0.29, mu<= -0.2, mu = mu +0.01,
nexp = 1.0;
dcoup =  Table[Table[0.1*Cos[2.0*Pi/6.0 * xpos[k]],{i,Nband/2},{j,Nband/2}],{k,1,NN}];
dcodn =  Table[Table[0.1*Cos[2.0*Pi/6.0 * xpos[k]],{i,Nband/2},{j,Nband/2}],{k,1,NN}];

dco = Join[dcoup,dcodn];(*vc = -0.185;*)
(*Table[dco[[i,j,j]]=0,{j,4},{i,2*NN}];*)
(*nexp = 1.0;*)



(*****Formation of Hamiltonian*****)
mdata = {};
 CloseKernels[];
 LaunchKernels[kernum];
For[a =1, a <= iter, a++, 
 H = Hop - Chempot[mu]+cdwMatrix[dco];
 ES = Eigensystem[H];(*calculating eigen vector and eigen values*)
 fermi = Map[f, ES[[1]]];(*calculating the fermi energy*)
 nutcalc = Table[(Abs[(ES[[2]][[;; , i]])]^2.) . fermi, {i, 1, NN*Nband/2}];
 ndtcalc = Table[((Abs[(ES[[2]][[;; , i + NN*Nband/2]])]^2.)) .  fermi, {i, 1, NN*Nband/2}];
 ni = nutcalc +ndtcalc;
 mi = nutcalc - ndtcalc;
 nt = Total[ni]/(NN*Nband/2);
 mt = Total[Abs[mi]]/(NN*Nband/2);

(*CloseKernels[];
LaunchKernels[kernum];
 *)

 (*order parameter for spin up*)
deltaABup = ParallelTable[Cos[ (2.0*Pi/6.0 *xpos[k1]) ]*2*Total[Table[Table[u0*(*KroneckerDelta[i,j]**)((Cos[ (2.0*Pi/6.0 *xpos[k2+1]) ]+1)*(Conjugate[ES[[2]][[;;,i]]]*ES[[2]][[;;,j]])) . fermi,
 {i,k2*Nband/2+1,k2*Nband/2+Nband/2},{j,k2*Nband/2+1,k2*Nband/2+Nband/2}],{k2,0,NN-1}]]/NN,{k1,1,NN}];

(*order parameter for spin down*)
deltaABdn =  ParallelTable[Cos[ (2.0*Pi/6.0 *xpos[k1]) ]*2*Total[Table[Table[u0*(*KroneckerDelta[i,j]**)((Cos[ (2.0*Pi/6.0 *xpos[k2+1-NN]) ]+1)*(Conjugate[ES[[2]][[;;,i]]]*ES[[2]][[;;,j]])) . fermi,
 {i,k2*Nband/2+1,k2*Nband/2+Nband/2},{j,k2*Nband/2+1,k2*Nband/2+Nband/2}],{k2,NN,2*NN-1}]]/NN,{k1,1,NN}];
 
 (*total order parameter that goes back into the hamiltonian*)
 dco = Join[deltaABup,deltaABdn];
 delta = Total[Table[Abs[deltaABup[[i]]],{i,Ny/2+1,Ny/2+Nx}]/Nx,4]/4; 
  
  (*try to keep the value of delta with in the range of 0.9 - 1.0 , 
  that generates a gap around 80-90 meV *)
 Print[a, " mu=", mu, " mt=", Chop[mt], " nt=", nt," deltaAB=",delta];
 (*mu = mu+0.25*(nexp-nt);*)
AppendTo[mdata,{"mu",mu,"mt",Chop[mt],"nt",nt,"deltaAB",delta}]];
Export["mfield_values.dat",mdata];
Export["cdwmat.wdx",cdwMatrix[dco]];
Export[StringJoin["WTe2_Hamiltonian_Mtarix_CDW_sm_2pi6_",ToString[mu],".wdx"],H];
]

Exit[];
Chop[cdwMatrix[dco]]//MatrixForm
