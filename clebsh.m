  BeginPackage[ "Clebsch`"]

  Clebsch::usage = 
    "Clebsch[ n, v, vbas] computes clebsch gordan coeff .(HALF FILLING)
    n: number of particles
    v: M representation (Ref: ?)
    vbas: some vector to represent the basis
    "
  do::usage =
  "
  stores the J_ operator
  "
  up::usage =
  "
  stores the J+ operator
  "
  jz::usage =
  "
  stores the Jz operator
  "
  sirrp::usage =
  "
  2S+1
  "
  cirrp::usage =
  "
  multiplicities of 2S+1
  "
  Begin[ "Private`"]

(* ******************************************************************
        The Main Program
   ******************************************************************
*)

  Clebsch[ n_,v_,vbas_,print_] := (
(*       Times @@ Map[Power @@ # &,f]
   *)
(*       FixedPointList [(# + 2/#)/2 &,x]
 *)
(*      FoldList [KroneckerProduct, First[x], Rest[x]];
  *)
    Module[
(*      Local variables                             *)
        {i,a,b,c,e,f,vir,bas,tmp,tmp2} ,

(*      Initialize variables                        *)
(*  
        c={{1,0},{0,1}};
        b={{0,0},{1,0}};
        e={{0,1},{0,0}};
        f={{1/2,0},{0,-1/2}};
  *)
        c=SparseArray[{{1,1}->1,{2,2}->1}];
        b=SparseArray[{{1,2}->0,{2,1}->1}];
        e=SparseArray[{{1,2}->1,{2,1}->0}];
        f=SparseArray[{{1,1}->1/2,{2,2}->-1/2}];
        do=ConstantArray[(Nest[KroneckerProduct[#,c]&,c,n-1]),{n}];
        up=ConstantArray[(Nest[KroneckerProduct[#,c]&,c,n-1]),{n}];
        jz=ConstantArray[(Nest[KroneckerProduct[#,c]&,c,n-1]),{n}];

(*      generate the basis                          *)
        bas=Nest[KroneckerProduct[(#+1),vbas]&,vbas,n-1];

(*      Find the Irreps for V(X)^N                  *)
        vir=Irreps[v,1,n];

(*      Irreps found stored in Vir                  *)

(*      Calculate Generators J_                     *)
        Do[
        Clear[a];
        a=ConstantArray[c,{n}];
        a[[i]]=b;
        do[[i]]=Apply[KroneckerProduct,a];
        ,{i,n}];
(*      Fin - Store generators in list d            *)

(*      Calculate Generators J+                     *)
        Do[
        Clear[a];
        a=ConstantArray[c,{n}];
        a[[i]]=e;
        up[[i]]=Apply[KroneckerProduct,a];
        ,{i,n}];
(*      Fin - Store generators in list up            *)

(*      Calculate Generators Jz                     *)
        Do[
        Clear[a];
        a=ConstantArray[c,{n}];
        a[[i]]=f;
        jz[[i]]=Apply[KroneckerProduct,a];
        ,{i,n}];
(*      Fin - Store generators in list up            *)

(*      Seperate the irreps generated, read Ms       *)
        sirrp=ConstantArray[-1,{Length[vir]}];
        cirrp=ConstantArray[1,{Length[vir]}];
        j=1;
        k=0;
        Do[
            If [sirrp[[j]]!=-1,
                k=0;
                Do[
                If[sirrp[[l]]!=vir[[i,1]],
                    k+=1,
                    cirrp[[l]]+=1;
                ];
                ,{l,j}];
                If[k==j,
                    j=j+1;
                    sirrp[[j]]=vir[[i,1]];
                ],
                sirrp[[j]]=vir[[1,1]];
            ];
        ,{i,Length[vir]}];

(*      seperated, purge the -1's                   *)
        sirrp2=sirrp;
        sirrp=DeleteCases[sirrp2, Alternatives @@ {-1}];

(*      sirrp contains the 2S+1 of different irreps *)
        sirrp2=sirrp;
        sirrp=Map[(#+1)&,sirrp2];
        cirrp[[Length[sirrp]]]=1;
        sirrp2=Take[cirrp,Length[sirrp]];
        cirrp=sirrp2;


(*      Do some output to the user                  *)

        Print ["
            Direct sum of Irreducible representatios :
            "];
        Print [vir];
        If[print,

        Print ["
            *** PROGRAM Clebsch ***\n 
            Version: 0.0.1\n
            HALF FILLING\n\n
            Number of spins: ",n];
        
        Print ["
            Size of the Hilbert space: ",Length[bas]];

        Print ["
            The basis of Hilbert space: "];
(*        Print [bas];*)

        Print ["
            Direct sum of Irreducible representatios :
            "];
        Print [vir];

        Print ["
            Dimension of Ms_min sub-space: ",Length[vir]];

        Print ["
            Irreps appearing in the Direct sum representation (2S+1): "];
        Print[sirrp];

        Print ["
            Multiplicities of above irreps in order as they appear : "];
        Print[cirrp];

        Print ["
            Generators (Only 1) jz
            "];
(*      Print [jz[[1]]];*)

        Print ["
            Generators (Only 1) j_
            "];
(*      Print [do[[1]]];*)

        Print ["
            Generators (Only 1) j+
            "];
(*      Print [up[[1]]];*)

        Print ["
            *** Fin ***"]

        ];

    ];
  );

(* ******************************************************************
        End of main program
   ******************************************************************
*)



(*      This function finds the Irreps
        it works recursively coupling
        one spin with the previous block
        this might not work for larger J            *)

(*      WARN WARN WARN WARN WARN WARN WARN          *)
(*      The current Irreps is designed only
        for S=1/2, it might need slight modifs
        for the general case                        *)
(*      WARN WARN WARN WARN WARN WARN WARN          *)


  Irreps[u_,j_,m_]:= (
      If[ j==m,
(*      Return the first (1,0)
        stopping condition J==M                     *)

          Return[u],

(*      Call Irreps recursively                     *)

          utmp=Irreps[u,j+1,m];

(*      Just add the newly found irrps to the list  *)
(*      should i use flatten !?                     *)

          utmp2=Apply[Join,Map[Join[{# - {1, 0}}, {# + {1, 0}}] &, utmp]];

(*      Here i remove the negative irreps 
        which are not physical                      *)

          tmp=utmp2;
          Do[
              If [i>Length[tmp],
                  Print["vijay",Length[tmp],i];
                  Break[];
              ];
              If [ tmp[[i,1]] <0,
                  tmp2=Drop[tmp,{i}];
                  tmp=tmp2;
              ];
          ,{i,Length[utmp2],1,-1}];

(*      Return the irreps at this level             *)

          utmp2=tmp;
          Return[utmp2];
        ];
  );

  End[]

  EndPackage[]
