ClearAll[closure];
closure[g_] := FixedPoint[
 Function[s, Union[Join[s, Flatten[Table[Simplify[x.y], {x, s}, {y, g}], 1]]]],
 {IdentityMatrix[2]}];
x = {{0, 1}, {1, 0}};
z = {{1, 0}, {0, -1}};
p = {{1, 0}, {0, I}};
h = (1 + I) {{1, 1}, {1, -1}}/2;
w = closure[{I IdentityMatrix[2], x, z}];
r = closure[{p, h}];
refs = Select[DeleteCases[w, IdentityMatrix[2]], Det[# - IdentityMatrix[2]] == 0 &];
cents = Table[Select[r, Simplify[#.s - s.#] == ConstantArray[0, {2, 2}] &], {s, refs}];
<|"WOrder" -> Length[w], "ROrder" -> Length[r],
 "P2Z" -> (p.p == z), "H2i" -> (Simplify[h.h] == I IdentityMatrix[2]),
 "HP3" -> (Simplify[MatrixPower[h.p, 3]] == -IdentityMatrix[2]),
 "ReflectionCount" -> Length[refs],
 "CentralizerOrders" -> (Length /@ cents),
 "CentralizersAbelian" -> And @@ (Function[c, And @@ Flatten[
 Table[Simplify[a.b - b.a] == ConstantArray[0, {2, 2}], {a, c}, {b, c}]]] /@ cents)|>

