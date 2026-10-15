# Si density reference

Generated from the Si diamond k444 tutorial rerun on 2026-09-05. ABACUS
1648a8a344427ae1b6394912bf677c4a20e053f2 and PyATB
9fb9028c59b1dbaf9cf66965280961fc2225d9eb produced the IBZ and full-grid
wave functions. This is an independent occupied-density comparison, not a
reference constructed by the symmetry routine under test.

Two atoms at fractional (0,0,0) and (1/4,1/4,1/4), FCC primitive cell;
22 orbitals per atom with l shells 0,0,0,1,1,1,2,2. Four occupied states;
D = C_occ^T C_occ^*. The reference file contains three cases, each with
9 row-convention rotation entries, 3 fractional translations, IBZ k, BZ k,
then IBZ and direct full-grid density matrices (44x44 complex, row-major).
Cases use zero-based (IBZ index, symmetry operation) (1,43), (3,1), (7,5).

The 1e-3 relative tolerance accommodates finite-grid rotational errors
(~2.4e-5 in the atom-swap case); the original library returns ~1.33.
The two other cases give ~2.2e-7 and ~4.4e-8. These checks do not establish
GW convergence or numerical acceptance of any corrected GW run.

The test checks atom-block and dense operator restoration as well as the
row-major wave-function transformation against the same direct full-grid
density. With `C_BZ = C_IBZ * R` and `D = C^T C^*`, the latter must satisfy
`D_BZ = R^T D_IBZ R^*`. It therefore checks the atom-swap block placement
against producer data, independently of the synthetic unit-test expectation.

The test also closes each fixture operation into a cyclic subgroup, prepares
its IBZ representatives on the full 4x4x4 mesh (including time reversal), and
calls `generate_kstars` and `build_kstar_member_rotations`. The automatically
chosen route is checked against the same producer density on all three
restoration paths. This exercises route selection and inversion as well as
atom routing; it does not claim coverage of all 48 diamond operations or of
the energy-weighted Green function and RPA contractions.
