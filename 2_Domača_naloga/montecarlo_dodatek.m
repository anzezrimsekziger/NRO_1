(* ::Package:: *)

MonteCarloPiSimulation[n_Integer] := Module[{pointsInsideCircle, totalPoints},
  pointsInsideCircle = 0;
  totalPoints = n;
  Do[
   {x, y} = RandomReal[{-1, 1}, 2];
   If[x^2 + y^2 <= 1, pointsInsideCircle++],
   {totalPoints}
   ];
  4 * pointsInsideCircle / totalPoints
]

