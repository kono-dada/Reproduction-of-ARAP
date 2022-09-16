(* ::Package:: *)

(*generate the mesh*)
\[ScriptCapitalR] = MeshRegion[{{0, 0}, {2, -1}, {1, 0}, {2, 1}}, Polygon[{1, 2, 3, 4}]];
region = BoundaryMeshRegion[{{0,0},{0,1},{1,1},{1,0}},Line[{1,2,3,4,1}]];
mesh = TriangulateMesh[\[ScriptCapitalR], MeshQualityGoal->"Maximal"];
nv = Length[MeshCells[mesh,0]];  (*number of vertice*)
nt = Length[MeshCells[mesh,2]];  (*number of triangles*)

(*Compute the coefficients in each triangle. 
Notice that the order of vertice in each triangle is counterclockwise, 
so if we store the coefficients in an right order, it will work well.*)
computeXAndY = {(#2-#1).(#3-#1)/Norm[#2-#1]^2,Det[{#2-#1,#3-#1}]/Norm[#2-#1]^2}&;
computeCoef[v1_,v2_,v3_]:={computeXAndY[v2,v3,v1],computeXAndY[v3,v1,v2],computeXAndY[v1,v2,v3]};
coefficientsInEachTriangle = (computeCoef@@Part[#,1])&/@MeshPrimitives[mesh,2];

(*In following steps, we aim to obtain matrix G of (5) in the paper*)
fixedVertice = {2,4,74};
nfixed=Length[fixedVertice];  (*number of fixed vertice*)
nfree=nv-nfixed;

error[v0_,v1_,v2_,x_,y_]:=(v0+x (v1-v0)+y({{0,-1},{1,0}}.(v1-v0))-v2).(v0+x (v1-v0)+y({{0,-1},{1,0}}.(v1-v0))-v2);  (*(4) in the paper*)
verticeIndexes=Array[#&,nv];
(*use 2*nv variables to obtain the polynomial of error, and compute its coefficients matrix*)
variablesOfFreeVertice=Flatten[{Subscript[v,#,"x"],Subscript[v,#,"y"]}&/@Select[verticeIndexes,!ContainsAll[fixedVertice,{#}]&]];
variablesOfFixedVertice=Flatten[{Subscript[v,#,"x"],Subscript[v,#,"y"]}&/@Select[verticeIndexes,ContainsAll[fixedVertice,{#}]&]];
errorOf3OrderedVertice[i_,j_,k_,x_,y_]:=error[
	{Subscript[v,i,"x"],Subscript[v,i,"y"]},
	{Subscript[v,j,"x"],Subscript[v,j,"y"]},
	{Subscript[v,k,"x"],Subscript[v,k,"y"]},
	x,y
];
errorOfATriangle[index_]:=Module[
	{coefi=coefficientsInEachTriangle[[index]],indexesOfVertice=MeshCells[mesh,2][[index]][[1]]},
	errorOf3OrderedVertice@@Join[RotateLeft[indexesOfVertice,1],coefi[[1]]]+
	errorOf3OrderedVertice@@Join[RotateRight[indexesOfVertice,1],coefi[[2]]]+
	errorOf3OrderedVertice@@Join[indexesOfVertice,coefi[[3]]]
];
thePolynomial=Total[errorOfATriangle/@Array[#&,nt]];
gInFormular5=CoefficientArrays[thePolynomial,Join[variablesOfFreeVertice,variablesOfFixedVertice]][[3]];

(*obtain G00 and G01*)
G=gInFormular5+Transpose[gInFormular5];
G00=Take[G,{1,2nfree},{1,2nfree}];
G01=Take[G,{1,2nfree},{2nfree+1,2nv}];
amazingMetrix=Inverse[G00].G01;

(*map new vertice to the mesh*)

DynamicModule[
	{newFixedVerticeCoordinate=MeshPrimitives[mesh,0][[#]][[1]]&/@fixedVertice},
	LocatorPane[
		Dynamic[newFixedVerticeCoordinate],
		Dynamic[
			u=-amazingMetrix.Flatten[newFixedVerticeCoordinate];
			newVerticeCoordinate=ArrayReshape[u,{nfree,2}];
			For[i=1,i<=nfixed,i++,newVerticeCoordinate=Insert[newVerticeCoordinate,newFixedVerticeCoordinate[[i]],fixedVertice[[i]]];];
			MeshRegion[
				newVerticeCoordinate,
				Triangle/@(#[[1]]&/@MeshCells[mesh,2]),
				MeshCellHighlight->(({0,#}->Red)&)/@fixedVertice,
				MeshCellLabel->(({0,#}->"Index")&)/@fixedVertice,
				MeshCellShapeFunction->(({0,#}->(Disk[#,0.03]&))&)/@fixedVertice,
				PlotRange->{{-0.5,2.5},{-1.5,1.5}}
			]
		]
	]
]






(* ::InheritFromParent:: *)
(**)
