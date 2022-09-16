(* ::Package:: *)

(*generate the mesh*)
\[ScriptCapitalR] = MeshRegion[{{0, 0}, {2, -1}, {1, 0}, {2, 1}}, Polygon[{1, 2, 3, 4}]];
region = BoundaryMeshRegion[{{-0.29,1.14},{-0.35,1.46},{0,1.7},{0.36,1.43},{0.22,1.1},{1,0.7},{0.9,0.42},{0.46,0.57},{0.44,0.25},{0.8,0},{0.64,-0.33},{0.2,0},{-0.2,0},{-0.43,-0.34},{-0.8,0},{-0.45,0.21},{-0.41,0.6},{-0.86,0.5},{-1,0.8}},Line[Append[Array[#&,19],1]]];
mesh = TriangulateMesh[region]
nv = Length[MeshCells[mesh,0]];  (*number of vertice*)
nt = Length[MeshCells[mesh,2]];  (*number of triangles*)
triangleCells=MeshCells[mesh,2];
verticeCells=MeshCells[mesh,0];


(* ::InheritFromParent:: *)
(**)


(* ::InheritFromParent:: *)
(**)


(*Compute the coefficients in each triangle. 
Notice that the order of vertice in each triangle is counterclockwise, 
so if we store the coefficients in an right order, it will work well.*)
computeXAndY = {(#2-#1).(#3-#1)/Norm[#2-#1]^2,Det[{#2-#1,#3-#1}]/Norm[#2-#1]^2}&;
computeCoef[v1_,v2_,v3_]:={computeXAndY[v2,v3,v1],computeXAndY[v3,v1,v2],computeXAndY[v1,v2,v3]};
coefficientsInEachTriangle = (computeCoef@@Part[#,1])&/@MeshPrimitives[mesh,2];

fixedVertice = {23,57,61,93,95,107,128,140,159,171};
nfixed=Length[fixedVertice];  (*number of fixed vertice*)
nfree=nv-nfixed;


(*In following steps, we aim to obtain matrix G of (5) in the paper*)
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


(*functions in step 2*)
(*the coordinates had been worked out*)
fittedVerticeInATriangle[x_,y_,v0newx_,v0newy_,v1newx_,v1newy_,v2newx_,v2newy_]:=Module[{
	v0x=-(-v0newx-v2newx+v1newx x+v2newx x-v0newx x^2-v1newx x^2-v1newy y+v2newy y-v0newx y^2-v1newx y^2)/(2 (1-x+x^2+y^2)),
	v0y=-(-v0newy-v2newy+v1newy x+v2newy x-v0newy x^2-v1newy x^2+v1newx y-v2newx y-v0newy y^2-v1newy y^2)/(2 (1-x+x^2+y^2)),
	v1x=-(-2 v1newx+v0newx x+2 v1newx x-v2newx x-v0newx x^2-v1newx x^2+v0newy y-v2newy y-v0newx y^2-v1newx y^2)/(2 (1-x+x^2+y^2)),
	v1y=-(-2 v1newy+v0newy x+2 v1newy x-v2newy x-v0newy x^2-v1newy x^2-v0newx y+v2newx y-v0newy y^2-v1newy y^2)/(2 (1-x+x^2+y^2))
},
	{{v0x,v0y},{v1x,v1y},{v0x,v0y}+x ({v1x,v1y}-{v0x,v0y})+y {{0,-1},{1,0}}.({v1x,v1y}-{v0x,v0y})}
]
normOfFirstVectors=Array[Norm[MeshPrimitives[mesh,2][[#]][[1]][[2]]-MeshPrimitives[mesh,2][[#]][[1]][[1]]]&,nt];
scaledTriangle[fittedTriangle_,i_]:=ScalingTransform[{normOfFirstVectors[[i]],normOfFirstVectors[[i]]}/Norm[fittedTriangle[[2]]-fittedTriangle[[1]]],Total[fittedTriangle]/3][fittedTriangle];



variablesOfFreeVertice2=Flatten[{Subscript[v',#,"x"],Subscript[v',#,"y"]}&/@Select[verticeIndexes,!ContainsAll[fixedVertice,{#}]&]];
variablesOfFixedVertice2=Flatten[{Subscript[v',#,"x"],Subscript[v',#,"y"]}&/@Select[verticeIndexes,ContainsAll[fixedVertice,{#}]&]];
variablesOfVerticefitted=Flatten[Array[{Subscript[vf,#,0,"x"],Subscript[vf,#,0,"y"],Subscript[vf,#,1,"x"],Subscript[vf,#,1,"y"],Subscript[vf,#,2,"x"],Subscript[vf,#,2,"y"]}&,nt]];
square = (#.#)&;
E2[i_]:=Module[
{indexOfv0=triangleCells[[i]][[1]][[1]],indexOfv1=triangleCells[[i]][[1]][[2]],indexOfv2=triangleCells[[i]][[1]][[3]]},
	vp0={Subscript[v',indexOfv0,"x"],Subscript[v',indexOfv0,"y"]};
	vp1={Subscript[v',indexOfv1,"x"],Subscript[v',indexOfv1,"y"]};
	vp2={Subscript[v',indexOfv2,"x"],Subscript[v',indexOfv2,"y"]};
	vf0={Subscript[vf,i,0,"x"],Subscript[vf,i,0,"y"]};
	vf1={Subscript[vf,i,1,"x"],Subscript[vf,i,1,"y"]};
	vf2={Subscript[vf,i,2,"x"],Subscript[vf,i,2,"y"]};
square[vp1-vp0-vf1+vf0]+square[vp2-vp1-vf2+vf1]+square[vp0-vp2-vf0+vf2]
]

coef2=CoefficientArrays[Total[Array[E2,nt]],Join[variablesOfFreeVertice2,variablesOfFixedVertice2]];
H=coef2[[3]];
matrixHp=Take[H+Transpose[H],{1,2nfree},{1,2nfree}];  (*important*)
inverseHp=Inverse[matrixHp];
matrixD=G01=Take[H+Transpose[H],{1,2nfree},{2nfree+1,2nv}];  (*important*)
f=coef2[[2]];
f0=Take[f,2nfree];  (*important*)


DynamicModule[
	{newFixedVerticeCoordinate=MeshPrimitives[mesh,0][[#]][[1]]&/@fixedVertice},  (*map new vertice to the mesh*)
	LocatorPane[
		Dynamic[newFixedVerticeCoordinate],
		Dynamic[
			u=-amazingMetrix.Flatten[newFixedVerticeCoordinate];
			newVerticeCoordinate=ArrayReshape[u,{nfree,2}];
			For[i=1,i<=nfixed,i++,newVerticeCoordinate=Insert[newVerticeCoordinate,newFixedVerticeCoordinate[[i]],fixedVertice[[i]]];];
			fittedTriangleCoordinate=Array[(fittedVerticeInATriangle@@Flatten[Join[coefficientsInEachTriangle[[#]][[3]],(newVerticeCoordinate[[#]]&)/@triangleCells[[#]][[1]]]])&,nt];
			scaledTriangles=Array[(scaledTriangle[fittedTriangleCoordinate[[#]],#]&),nt];
			uu=-inverseHp.(matrixD.Flatten[newFixedVerticeCoordinate]+Normal[f0]/.AssociationThread[variablesOfVerticefitted,Flatten[scaledTriangles]]);
			newVerticeCoordinate2=ArrayReshape[uu,{nfree,2}];
			For[i=1,i<=nfixed,i++,newVerticeCoordinate2=Insert[newVerticeCoordinate2,newFixedVerticeCoordinate[[i]],fixedVertice[[i]]];];
			Show[
			MeshRegion[
				newVerticeCoordinate2,
				Triangle/@(#[[1]]&/@triangleCells),
				MeshCellHighlight->(({0,#}->Red)&)/@fixedVertice,
				MeshCellLabel->(({0,#}->"Index")&)/@fixedVertice,
				MeshCellShapeFunction->(({0,#}->(Disk[#,0.03]&))&)/@fixedVertice,
				PlotRange->{{-2,2},{-1,2.5}},
				ImageSize->Full
			]
			]
		]
	]
]
(*the following are optional graphic*)

(*Graphics[{Opacity[0.2],Gray,Triangle/@scaledTriangles}]*)
(*MeshRegion[
				newVerticeCoordinate,
				Triangle/@(#[[1]]&/@triangleCells),
				MeshCellHighlight->(({0,#}->Red)&)/@fixedVertice,
				MeshCellLabel->(({0,#}->"Index")&)/@fixedVertice,
				MeshCellShapeFunction->(({0,#}->(Disk[#,0.03]&))&)/@fixedVertice,
				PlotRange->{{-2.5,2.5},{-1.5,2.5}},
				ImageSize->Full,
				MeshCellStyle\[Rule]{{2,All}->Opacity[0.5,Pink]}
			],*)

