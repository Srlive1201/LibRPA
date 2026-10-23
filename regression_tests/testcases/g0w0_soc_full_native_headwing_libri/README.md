# Full-BZ SOC native optical control

This finite four-state model supplies all three BZ points, original complex
spinor wavefunctions and H/S/r operators. The native public optical API retains
the supplied SCF eigenframe. The complete complex Cartesian velocity is checked
in Ha Bohr, including all diagonal and off-diagonal entries: one spin channel,
three k points, three axes and four by four bands (144 complex entries).

The velocity reference comes from an independent NumPy/SciPy construction of
this nonorthogonal, time-reversal-covariant four-state model. Generalized SCF
eigenvectors include a complex rotation within each Kramers pair. Direct Fourier
derivatives of H(R) and S(R), together with r(R), are projected into that retained
SCF frame; the Hermitian upper triangle completes the velocity. All source
dataset files match the original independent fixture. The reference is not an
output of the native LibRPA velocity kernel. The existing table comparator
checks every index and real/imaginary component with absolute tolerance 1e-10,
and rejects missing rows and nonfinite values.

Analytic head/wing and the complete complex Sigma at the two common reference
points are compared with the independently retained full-grid reference.
No material-convergence conclusion follows.

The IBZ-only companion is introduced with magnetic restoration later in the series.
