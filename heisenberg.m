#!/usr/local/bin/MathematicaScript -script

(* generate CGC *)
Print["Heisenberg"]

(*  AppendTo[$Path, ToFileName[$HomeDirectory,$]];
  *)
SetDirectory[NotebookDirectory[]];
Get["clebsh`"];

v={{1,0}};
vbas={{p},{m}};
n=4;
Clebsch[n,v,vbas,False];
imat=SparseArray[{{1,1}->0,{2,2}->0}];
H=(Nest[KroneckerProduct[#,imat]&,imat,n-1]);

(*  
    calculate Hamil
    H = (S-S+ + S+S-)/2 + S_z*S_z
  *)

Do[

    H=H + 1/2( Dot[up[[i]],do[[i+1]]] + Dot[do[[i]],up[[i+1]]])
            + 1/2 ( Dot[jz[[i]],jz[[i+1]]] + Dot[jz[[i]],jz[[i+1]]] )

,{i,n-1}];

pl=MatrixPlot[N@H];
Export["mat.pdf",pl];
{evals,evec}=Eigensystem[N[H]];
{alam,alpha}=Transpose@SortBy[Transpose[{evals,evec}],First];

Print@alam;
(*  
    calculate S_z
    <ψ|S_z|ψ>
  *)

Do[
  tmp=Plus @@ (Map[Dot[#,alpha[[i]]]&,jz]);
  tmp2=Dot[alpha[[i]],tmp];
  Print@tmp2;
,{i,Power[2,n]}];

Print["-------"];
(*  
    calculate S^2
    <ψ|S^2|ψ> = <ψ|S-S+ + S_z(S_z + 1)|ψ>
  *)

Do[
  tmpsz=  Plus @@ (Map[Dot[#,alpha[[i]]]&,jz]);
  tmpszsz=Plus @@ (Map[Dot[#,tmpsz]&,jz]);
  tmpsup= Plus @@ (Map[Dot[#,alpha[[i]]]&,up]);
  tmpsupsdo=Plus @@ (Map[Dot[#,tmpsup]&,do]);
  tmptot=tmpszsz+tmpsz+tmpsupsdo;
  tmpfin=Dot[alpha[[i]],tmptot];
  Print@tmpfin;
,{i,Power[2,n]}];













(*  
Run["rm inp"];
Export["inp",Normal[tmp1-tmp2],"Table"];
Run["./print.py"];
  Print@Normal[tmp1-tmp2];
tmp=Normalize@{1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0};
Print[Plus@@(Map[Dot[#,tmp]&,jz])];
Print[Dot[tmp,Plus@@(Map[Dot[#,tmp]&,jz])]];
eigen=Eigenvalues[N[H]]

Print@eigen[[1]]
Print@Eigenvectors[N[H]][[2]]
Print@eigen[[2]]
Print@(eigen[[2]]-eigen[[1]])

Print@H
H=1/2( Dot[up[[1]],do[[2]]] + Dot[do[[1]],up[[2]]]  + Dot[up[[2]],do[[3]]] + Dot[do[[2]],up[[3]]] )
sz = {{1/2,0},{0,-1/2}};
sm = {{0,0},{1,0}};
sp = {{0,1},{0,0}};
imat = {{1,0},{0,1}};

sz2=KroneckerProduct[sz,sz];
smsp=KroneckerProduct[sm,sp];
spsm=KroneckerProduct[sp,sm];
imat=KroneckerProduct[imat,imat];
sz2BB=ConstantArray[(Nest[KroneckerProduct[#,imat]&,imat,2-1]),{2}];

Orot=imat;

a=ConstantArray[imat,{2}];

Do[
    smsp=Dot[Dot[Orot,smsp],Transpose[Orot]];
    spsm=Dot[Dot[Orot,spsm],Transpose[Orot]];
    sz2=Dot[Dot[Orot,sz2],Transpose[Orot]];

(*      Calculate Generators J_                     *)
    Do[
        Clear[a];
        a=ConstantArray[imat,{i}];
        a[[j]]=sz2;
        sz2BB[[j]]=Apply[KroneckerProduct,a];
    ,{j,2}];
(*      Fin - Store generators in list d            *)

,{i,2,2}]
  *)
