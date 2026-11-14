"""Independent round-three check of the missing small-field 3'-parameters.

This is a finite root-datum computation, not a block-theoretic calculation.
It constructs the F4 roots, their character lattice, and all 1152 Weyl
matrices from scratch.  It does not import a previous F4 certificate.

The simple roots, in the usual orthonormal Euclidean basis, are
  a1=e2-e3, a2=e3-e4, a3=e4, a4=(e1-e2-e3-e4)/2.
Matrices act on column vectors of character coordinates in this basis of
simple roots.  A torus point y modulo the prime r represents the character
lattice homomorphism x |-> zeta_r^(x^T y).  For the rational twist with
F*=epsilon*q*w on characters, rationality is
  (epsilon*q*w-I)^T y = 0 (mod r).

The two choices w=s1*s2 and w=s3*s4 rotate the long and the short A2,
respectively.  For epsilon=-1 the actual Weyl twist is w0*w=-w.
The chosen points have orders 7 at q=4 and 19 at q=8.  Their connected
root centralizers are the perpendicular A2, and the entire point
stabilizer in W is exactly its Weyl group, with all six elements rational.

The AHL rational type labels in the output are a comparison with its
Tables 23 and 25: long A2 is type 13, short A2 is type 17; k=5 at q=4
and k=4 at q=8.  Identifying the resulting blocks and their *dual* actual
defect groups still uses AHL Tables 13 and 17 and Proposition 6.18.
"""

from itertools import combinations, product
from pathlib import Path
import json

import sympy as sp


IDENTITY = sp.eye(4)
SIMPLE = sp.Matrix.hstack(
    *[
        sp.Matrix(v)
        for v in (
            (0, 1, -1, 0),
            (0, 0, 1, -1),
            (0, 0, 0, 1),
            (sp.Rational(1, 2), sp.Rational(-1, 2),
             sp.Rational(-1, 2), sp.Rational(-1, 2)),
        )
    ]
)
SIMPLE_INVERSE = SIMPLE.inv()


def reflection(vector):
    vector = sp.Matrix(vector)
    return IDENTITY - 2 * vector * vector.T / vector.dot(vector)


def generated_group(generators):
    seen = {tuple(IDENTITY)}
    frontier = [IDENTITY]
    while frontier:
        element = frontier.pop()
        for generator in generators:
            new_element = generator * element
            key = tuple(new_element)
            if key not in seen:
                seen.add(key)
                frontier.append(new_element)
    return [sp.Matrix(4, 4, list(key)) for key in sorted(seen)]


def nullspace_mod_prime(matrix, prime):
    rows = [[int(matrix[i, j]) % prime for j in range(matrix.cols)]
            for i in range(matrix.rows)]
    pivots = []
    active_row = 0
    for column in range(matrix.cols):
        pivot_row = next((j for j in range(active_row, matrix.rows)
                          if rows[j][column]), None)
        if pivot_row is None:
            continue
        rows[active_row], rows[pivot_row] = rows[pivot_row], rows[active_row]
        scalar = pow(rows[active_row][column], -1, prime)
        rows[active_row] = [entry * scalar % prime for entry in rows[active_row]]
        for j in range(matrix.rows):
            if j != active_row:
                scalar = rows[j][column]
                rows[j] = [(left - scalar * right) % prime
                           for left, right in zip(rows[j], rows[active_row])]
        pivots.append(column)
        active_row += 1
    basis = []
    for free_column in range(matrix.cols):
        if free_column in pivots:
            continue
        vector = [0] * matrix.cols
        vector[free_column] = 1
        for row, column in enumerate(pivots):
            vector[column] = -rows[row][free_column] % prime
        basis.append(sp.Matrix(vector))
    return basis


generators = [SIMPLE_INVERSE * reflection(SIMPLE[:, i]) * SIMPLE
              for i in range(4)]
assert all(entry.is_Integer for matrix in generators for entry in matrix)
weyl = generated_group(generators)
assert len(weyl) == 1152

euclidean_roots = []
for i in range(4):
    for sign in (-1, 1):
        vector = [0] * 4
        vector[i] = sign
        euclidean_roots.append(sp.Matrix(vector))
for i, j in combinations(range(4), 2):
    for first_sign, second_sign in product((-1, 1), repeat=2):
        vector = [0] * 4
        vector[i], vector[j] = first_sign, second_sign
        euclidean_roots.append(sp.Matrix(vector))
euclidean_roots += [sp.Matrix(vector) / 2
                    for vector in product((-1, 1), repeat=4)]
roots = [SIMPLE_INVERSE * root for root in euclidean_roots]
assert len({tuple(root) for root in roots}) == 48
assert all(entry.is_Integer for root in roots for entry in root)

cases = []
for name, twist in (
    ("long_A2_twist", generators[0] * generators[1]),
    ("short_A2_twist", generators[2] * generators[3]),
):
    assert twist**3 == IDENTITY and twist != IDENTITY
    for q, epsilon, prime in ((4, 1, 7), (8, -1, 19)):
        assert sp.isprime(prime) and prime != 3
        rational_matrix = (epsilon * q * twist - IDENTITY).T
        basis = nullspace_mod_prime(rational_matrix, prime)
        assert len(basis) == 1
        point = basis[0]
        assert any(int(entry) % prime for entry in point)
        assert all(int(entry) % prime == 0 for entry in rational_matrix * point)
        torus_order = abs(int((epsilon * q * twist - IDENTITY).det()))
        assert torus_order == (q - epsilon) * (q**3 - epsilon)

        centralizer_indices = [i for i, root in enumerate(roots)
                               if int((root.T * point)[0]) % prime == 0]
        centralizer_roots = [roots[i] for i in centralizer_indices]
        assert len(centralizer_roots) == 6
        assert sp.Matrix.hstack(*centralizer_roots).rank() == 2
        lengths = {int(2 * euclidean_roots[i].dot(euclidean_roots[i]))
                   for i in centralizer_indices}
        assert len(lengths) == 1
        assert all(twist * root == root for root in centralizer_roots)

        root_weyl = generated_group(
            [SIMPLE_INVERSE * reflection(euclidean_roots[i]) * SIMPLE
             for i in centralizer_indices]
        )
        stabilizer = [element for element in weyl
                      if all(int(entry) % prime == 0
                             for entry in element.T * point - point)]
        assert len(root_weyl) == len(stabilizer) == 6
        assert {tuple(element) for element in root_weyl} == {
            tuple(element) for element in stabilizer
        }
        assert any(left * right != right * left
                   for left in root_weyl for right in root_weyl)
        assert all(element * twist == twist * element for element in stabilizer)
        geometric_type = 13 if lengths == {4} else 17
        cases.append({
            "q": q,
            "epsilon": epsilon,
            "parameter_order": prime,
            "twist": name,
            "twist_matrix_on_characters": [[int(entry) for entry in row]
                                            for row in twist.tolist()],
            "torus_order": torus_order,
            "torus_point_mod_parameter_order": [int(entry) for entry in point],
            "annihilating_roots_simple_coordinates": [
                [int(entry) for entry in root] for root in centralizer_roots
            ],
            "twice_root_squared_length": next(iter(lengths)),
            "centralizer_root_rank": 2,
            "point_stabilizer_in_W_order": len(stabilizer),
            "rational_relative_W": "S3",
            "ahl_rational_type_comparison": [geometric_type, 5 if q == 4 else 4],
        })

result = {
    "simple_roots_euclidean_columns": [[str(entry) for entry in SIMPLE[:, i]]
                                       for i in range(4)],
    "root_count": len(roots),
    "weyl_group_order": len(weyl),
    "cases": cases,
    "scope": (
        "Exact rational 3'-parameter existence, root centralizers and rational "
        "Weyl inertia at q=4,8. This does not identify Brauer blocks, prove "
        "ordinary-character exhaustivity, or prove iAM. Those use the separate "
        "arguments and sources in type_f4.tex."
    ),
}
Path(__file__).with_suffix(".json").write_text(
    json.dumps(result, indent=2) + "\n", encoding="utf-8"
)
print(json.dumps(result, indent=2))
