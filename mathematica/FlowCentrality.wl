(* ::Package:: *)

(* :Title: FlowCentrality *)
(* :Context: FlowCentrality` *)
(* :Author: Jacob A. Elliott *)
(* :Summary: Spectral-tilt occupancy, Action--Betweenness, and LDP rate functions on quotient flows. *)
(* :Package Version: 0.1 *)
(* :Mathematica Version: 13.0+ *)
(* :Keywords: flow centrality, spectral tilt, large deviations, quotient flows, action--betweenness *)

BeginPackage["FlowCentrality`"];

QuotientKernel::usage =
  "QuotientKernel[P_, blocks_] -> {Pdown, disp}, where P is row-stochastic SparseArray and blocks is a list of index lists.";

TiltedSpectral::usage =
  "TiltedSpectral[P_, f_, t_] -> {logRho, rt, lt}, Perron log spectral radius and right/left eigenvectors of diag(Exp[t f]) . P.";

OccCentrality::usage =
  "OccCentrality[P_, f_, t_] -> occupancy derivative d/dt log rho at tilt t.";

ActionBetweenness::usage =
  "ActionBetweenness[P_, f_, t_] -> Association of edge -> sensitivity scores lt[u]*rt[v] at tilt t.";

RateFunction::usage =
  "RateFunction[P_, f_, theta_, tGuess_:0.0] -> I(theta) via Legendre dual sup_t {t theta - log rho(t)}.";

Begin["`Private`"];

rowNormalize[m_] := Module[{rs = Total[m, {2}]},
  SparseArray[If[#1 == 0, 0, #2/#1] & @@@ Thread[{rs, m}], Dimensions[m]]
];

QuotientKernel[P_SparseArray, blocks_List] := Module[
  {n = Length[blocks], m, idx, Pdown, disp, rows, blk, tot, pRow, targ, numer, denom, i, j},
  Pdown = ConstantArray[0., {n, n}];
  disp = 0.;
  Do[
    blk = blocks[[i]];
    rows = P[[blk, All]];
    tot = Total[rows, {1}];
    (* For each target block j, average transitions from blk to blocks[[j]] *)
    Do[
      idx = blocks[[j]];
      numer = Total[rows[[All, idx]]];
      denom = Total[rows, 2]; (* per-row sums = 1 *)
      pRow = If[Length[blk] > 0, numer/Length[blk], 0.];
      Pdown[[i, j]] = Total[pRow],
      {j, 1, n}
    ];
    (* Dispersion: max L1 row discrepancy to block average *)
    disp = Max[disp, Max[Total[Abs[rows - ConstantArray[Mean[rows], Dimensions[rows]]], {2}]]],
    {i, 1, n}
  ];
  {rowNormalize@SparseArray[Pdown], disp}
];

TiltedSpectral[P_SparseArray, f_List, t_?NumericQ] := Module[
  {d = Exp[t f], Pt, esR, esL, rt, lt, rho},
  Pt = SparseArray[DiagonalMatrix[d]].P;
  (* right Perron eigenvector *)
  {rho, rt} = Eigensystem[Normal[Pt], 1, Method -> {"Arnoldi"}][[All, 1]];
  (* left Perron eigenvector = Perron of transpose *)
  {rho, lt} = Eigensystem[Normal[Transpose[Pt]], 1, Method -> {"Arnoldi"}][[All, 1]];
  lt = lt/ (lt.rt); (* normalize lt^T rt = 1 *)
  {Log[rho], rt, lt}
];

OccCentrality[P_SparseArray, f_List, t_?NumericQ] := Module[
  {logR, rt, lt},
  {logR, rt, lt} = TiltedSpectral[P, f, t];
  lt.(f*rt)
];

ActionBetweenness[P_SparseArray, f_List, t_?NumericQ] := Module[
  {logR, rt, lt, edges, scores, dims = Dimensions[P], r, c},
  {logR, rt, lt} = TiltedSpectral[P, f, t];
  edges = ArrayRules[P];
  scores = Association[];
  Do[
    {r, c} = edges[[k, 1]];
    If[edges[[k, 2]] > 0,
      scores[{r, c}] = lt[[r]]*rt[[c]];
    ],
    {k, 1, Length[edges]}
  ];
  scores
];

RateFunction[P_SparseArray, f_List, theta_?NumericQ, tGuess_:0.0] := Module[
  {obj, tOpt},
  obj[t_] := t*theta - First[TiltedSpectral[P, f, t]];
  tOpt = Quiet@FindMaximum[{obj[t], -10 <= t <= 10}, {t, tGuess}, 
           Method -> {"Newton", "PrincipalAxis"}][[2, 1, 2]];
  obj[tOpt]
];

End[];
EndPackage[];