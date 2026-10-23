"""Load exported scalar or SOC QSGW Hamiltonians with the existing PyATB readers."""


def load_qsgw_model(
    hr_file,
    sr_file="srs1_nao.csr",
    rr_file="rr.csr",
    *,
    lattice_constant_angstrom,
    lattice_vector,
    n_kpoints,
    rr_unit="Bohr",
    nspin=1,
):
    """Return a scalar-channel or full SOC PyATB model from a saved QSGW HR.

    ``hr_file`` is the complete complex QSGW CSR written in Ry. ``sr_file``
    and ``rr_file`` are the matching original ABACUS S(R) and r(R),
    with r's length unit specified by ``rr_unit`` (``Bohr`` or ``Angstrom``).
    ``lattice_vector`` contains three dimensionless row vectors; multiplying
    them by the positive ``lattice_constant_angstrom`` gives the cell in
    Angstrom. ``n_kpoints`` sets PyATB's maximum k-point batch size.

    Use ``nspin=1`` for a scalar model or one selected collinear spin channel,
    and ``nspin=4`` for the complete SOC spinor model. SOC H/S/r must all use
    the producer's full spinor dimension and AO-major/spinor-minor ordering;
    no scalar overlap or position matrix is replicated implicitly. H always
    uses complex token parsing; SOC S also uses the native complex reader.
    ABACUS r uses real CSR entries for both scalar and SOC data. All H/S/r
    R vectors are retained, with exact zero padding
    to their common union. Unit conversions follow the installed PyATB readers.
    No file is modified and no interpolation policy,
    occupied-subspace choice, or topology algorithm is introduced here.

    NumPy, SciPy and PyATB are needed only when this function is called.
    """
    import operator
    import os

    import numpy as np
    from scipy.sparse import csc_matrix
    from pyatb.io.abacus_read_xr import abacus_readHR, abacus_readSR, abacus_readrR
    from pyatb.tb import tb
    from pyatb.tb.multixr import multiXR

    lattice_constant = float(lattice_constant_angstrom)
    lattice = np.asarray(lattice_vector, dtype=float)
    if not np.isfinite(lattice_constant) or lattice_constant <= 0:
        raise ValueError("lattice_constant_angstrom must be finite and positive")
    if lattice.shape != (3, 3) or not np.isfinite(lattice).all():
        raise ValueError("lattice_vector must be a finite 3 by 3 array of row vectors")
    if np.linalg.matrix_rank(lattice) != 3:
        raise ValueError("lattice_vector must describe a nonsingular cell")
    if isinstance(n_kpoints, (bool, np.bool_)):
        raise ValueError("n_kpoints must be a positive integer")
    try:
        max_kpoints = operator.index(n_kpoints)
    except TypeError as error:
        raise ValueError("n_kpoints must be a positive integer") from error
    if max_kpoints <= 0:
        raise ValueError("n_kpoints must be a positive integer")
    if rr_unit not in ("Bohr", "Angstrom"):
        raise ValueError("rr_unit must be Bohr or Angstrom")
    if isinstance(nspin, (bool, np.bool_)):
        raise ValueError("nspin must be 1 for one scalar channel or 4 for SOC")
    try:
        nspin = operator.index(nspin)
    except TypeError as error:
        raise ValueError("nspin must be 1 for one scalar channel or 4 for SOC") from error
    if nspin not in (1, 4):
        raise ValueError("nspin must be 1 for one scalar channel or 4 for SOC")

    h = abacus_readHR(4, os.fspath(hr_file), "Ry")
    s = abacus_readSR(nspin, os.fspath(sr_file))
    r = abacus_readrR(os.fspath(rr_file), rr_unit)
    if len(r) != 3:
        raise ValueError("r(R) must contain three Cartesian components")
    sources = [h, s, *r]
    n = h.basis_num
    if n <= 0 or any(source.basis_num != n for source in sources):
        raise ValueError("H(R), S(R) and r(R) must have the same positive AO dimension")
    if nspin == 4 and n % 2:
        raise ValueError("SOC H/S/r must have an even complete spinor dimension")
    packed_dimension = n * (n + 1) // 2
    coordinates = []
    values = []
    for source in sources:
        vectors = np.asarray(source.R_direct_coor)
        data = source.XR
        if (
            vectors.shape != (source.R_num, 3)
            or not np.isfinite(vectors).all()
            or not np.equal(vectors, np.rint(vectors)).all()
        ):
            raise ValueError("H/S/r must use finite integer real-space cell vectors")
        vectors = vectors.astype(np.int64)
        if len(set(map(tuple, vectors))) != len(vectors):
            raise ValueError("H/S/r real-space cell vectors must be unique")
        if data.shape != (source.R_num, packed_dimension) or not np.isfinite(data.data).all():
            raise ValueError("H/S/r upper-triangle arrays must have matching finite entries")
        coordinates.append(vectors)
        values.append(data)

    union = np.array(
        sorted(set().union(*(set(map(tuple, vectors)) for vectors in coordinates))),
        dtype=np.int64,
    )
    if union.size == 0:
        raise ValueError("H/S/r require a nonempty real-space support")
    union_rows = {tuple(cell): index for index, cell in enumerate(union)}
    aligned_operators = []
    for name, vectors, data in zip(["H", "S", "r_x", "r_y", "r_z"], coordinates, values):
        row_map = np.array([union_rows[tuple(cell)] for cell in vectors], dtype=np.int64)
        entries = data.tocoo(copy=False)
        aligned = csc_matrix(
            (entries.data, (row_map[entries.row], entries.col)),
            shape=(len(union), packed_dimension), dtype=complex,
        )
        operator_matrix = multiXR(name)
        operator_matrix.set_XR(len(union), union, n, aligned)
        aligned_operators.append(operator_matrix)

    model = tb(nspin, lattice_constant, lattice, max_kpoints)
    model.set_solver_HSR(*aligned_operators[:2], isSparse=False)
    model.set_solver_rR(*aligned_operators[2:], isSparse=False)
    return model
