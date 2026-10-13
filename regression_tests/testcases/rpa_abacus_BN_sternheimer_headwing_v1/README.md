# Compact Delta-ST RPA regression

This case reuses the established BN reader-v1 head/wing dataset and adds a
small synthetic Delta-Sternheimer fixture under `librpa/fixtures`.

The fixture contains two q points, two imaginary frequencies per q, and 34 by
34 Hermitian response matrices in the existing two-atom auxiliary basis. The
matrix values are deterministic tridiagonal software-test data; they are not a
physical BN response. Gamma uses analytic q-averaged head/wing while the second
q point uses the ordinary body. The checked result is the summed
`Total Sternheimer EcRPA`.

Run `generate_sternheimer_fixture.py` from this directory to regenerate the
visible fixture copy. Pass an extracted `dataset/` path to update the files that
are packed into the shared BN archive.
