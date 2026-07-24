import math
import re


__all__ = [
    "eigenvalue_trace",
    "iteration_summary",
]


HA2EV = 27.211386245988

EIGENVALUE_HEADER = (
    "# iter channel spin kpoint kx ky kz band energy_eV"
)
SUMMARY_HEADER = (
    "# iter max_delta_eV residual_l2_Ha residual_max_Ha "
    "efermi_eV gap_eV electron_count requested_mode applied_mode beta "
    "fallback rcond coefficient_l1 coefficient_count converged "
    "coefficients fallback_reason"
)

CONTRACT_KEYS = frozenset((
    "qsgw_contract_version",
    "fixed_basis",
    "live_update",
    "velocity",
    "head",
    "wing",
    "symmetry",
    "hartree",
    "band",
    "h_qsgw_cut",
    "qsgw_band0_unoccupied_keep",
    "qsgw_band0_cut_mode",
    "qsgw_band0_cut_shift_ha",
    "qsgw_input_contract",
    "qsgw_input_contract_sha256",
    "qsgw_mixer",
    "qsgw_mixing_beta",
))
REQUIRED_CONTRACT_KEYS = frozenset((
    "qsgw_contract_version",
    "fixed_basis",
    "live_update",
    "velocity",
    "head",
    "wing",
    "symmetry",
    "hartree",
    "band",
    "h_qsgw_cut",
    "qsgw_input_contract",
    "qsgw_input_contract_sha256",
    "qsgw_mixer",
    "qsgw_mixing_beta",
))


def eigenvalue_trace(tolerance_ha="1e-6", coordinate_tolerance="1e-12",
                     precision="3"):
    """Compare keyed QSGW eigenvalue traces using a Hartree tolerance."""
    tolerance_ha = _positive_float(
        tolerance_ha, "tolerance_ha", allow_zero=True)
    coordinate_tolerance = _positive_float(
        coordinate_tolerance, "coordinate_tolerance", allow_zero=True)
    precision = int(precision)

    def inner(test_files, reference_files):
        try:
            pairs = _paired_file_texts(test_files, reference_files)
            maximum_energy = 0.0
            maximum_coordinate = 0.0
            maximum_location = None
            nvalues = 0
            for filename, test_text, reference_text in pairs:
                contract = _require_same_contract(
                    test_text, reference_text, filename)
                test_rows = _parse_eigenvalue_trace(
                    test_text, "test {}".format(filename))
                reference_rows = _parse_eigenvalue_trace(
                    reference_text, "reference {}".format(filename))
                _validate_eigenvalue_trajectory(
                    test_rows, contract, "test {}".format(filename))
                _validate_eigenvalue_trajectory(
                    reference_rows, contract,
                    "reference {}".format(filename))
                mismatch = _key_mismatch(test_rows, reference_rows)
                if mismatch is not None:
                    return False, (
                        "{}: eigenvalue key set mismatch; {}"
                        .format(filename, mismatch)
                    )
                for key in sorted(test_rows):
                    test_coordinate, test_energy = test_rows[key]
                    reference_coordinate, reference_energy = \
                        reference_rows[key]
                    coordinate_difference = max(
                        abs(x - y) for x, y in
                        zip(test_coordinate, reference_coordinate))
                    if coordinate_difference > coordinate_tolerance:
                        return False, (
                            "{}: k-point coordinate mismatch at {}: "
                            "{:.6E} > {:.6E}"
                            .format(filename, key, coordinate_difference,
                                    coordinate_tolerance)
                        )
                    energy_difference = (
                        abs(test_energy - reference_energy) / HA2EV
                    )
                    if energy_difference > maximum_energy:
                        maximum_energy = energy_difference
                        maximum_location = (filename, key)
                    maximum_coordinate = max(
                        maximum_coordinate, coordinate_difference)
                    nvalues += 1
            message = (
                "max abs eigenvalue diff = {:.{p}E} Ha "
                "(tol = {:.{p}E} Ha), max k-point coordinate diff = "
                "{:.{p}E} over {} values"
            ).format(maximum_energy, tolerance_ha, maximum_coordinate,
                     nvalues, p=precision)
            if maximum_location is not None:
                message += ", max at {} key {}".format(*maximum_location)
            return maximum_energy <= tolerance_ha, message
        except (TypeError, ValueError) as error:
            return False, str(error)

    return inner


def iteration_summary(energy_tolerance_ev="1e-5",
                      residual_tolerance_ha="1e-8",
                      scalar_tolerance="1e-10",
                      coefficient_tolerance="1e-12", precision="3"):
    """Compare QSGW iteration summaries with field-specific tolerances."""
    energy_tolerance_ev = _positive_float(
        energy_tolerance_ev, "energy_tolerance_ev", allow_zero=True)
    residual_tolerance_ha = _positive_float(
        residual_tolerance_ha, "residual_tolerance_ha", allow_zero=True)
    scalar_tolerance = _positive_float(
        scalar_tolerance, "scalar_tolerance", allow_zero=True)
    coefficient_tolerance = _positive_float(
        coefficient_tolerance, "coefficient_tolerance", allow_zero=True)
    precision = int(precision)

    energy_fields = ("max_delta_eV", "efermi_eV", "gap_eV")
    residual_fields = ("residual_l2_Ha", "residual_max_Ha")
    scalar_fields = ("electron_count", "beta", "rcond")
    exact_fields = (
        "requested_mode", "applied_mode", "fallback", "coefficient_count",
        "converged", "fallback_reason",
    )

    def inner(test_files, reference_files):
        try:
            pairs = _paired_file_texts(test_files, reference_files)
            maxima = {
                name: 0.0
                for name in energy_fields + residual_fields + scalar_fields
            }
            maximum_coefficient = 0.0
            niterations = 0
            for filename, test_text, reference_text in pairs:
                _require_same_contract(test_text, reference_text, filename)
                test_rows = _parse_iteration_summary(
                    test_text, "test {}".format(filename))
                reference_rows = _parse_iteration_summary(
                    reference_text, "reference {}".format(filename))
                mismatch = _key_mismatch(test_rows, reference_rows)
                if mismatch is not None:
                    return False, (
                        "{}: iteration key set mismatch; {}"
                        .format(filename, mismatch)
                    )
                for iteration in sorted(test_rows):
                    test_row = test_rows[iteration]
                    reference_row = reference_rows[iteration]
                    for field in exact_fields:
                        if test_row[field] != reference_row[field]:
                            return False, (
                                "{}: iteration {} {} mismatch: {} != {}"
                                .format(filename, iteration, field,
                                        test_row[field],
                                        reference_row[field])
                            )
                    for fields, tolerance in (
                        (energy_fields, energy_tolerance_ev),
                        (residual_fields, residual_tolerance_ha),
                        (scalar_fields, scalar_tolerance),
                    ):
                        for field in fields:
                            difference = abs(
                                test_row[field] - reference_row[field])
                            maxima[field] = max(maxima[field], difference)
                            if difference > tolerance:
                                return False, (
                                    "{}: iteration {} {} difference "
                                    "{:.6E} exceeds {:.6E}"
                                    .format(filename, iteration, field,
                                            difference, tolerance)
                                )
                    test_coefficients = test_row["coefficients"]
                    reference_coefficients = \
                        reference_row["coefficients"]
                    if len(test_coefficients) != len(
                            reference_coefficients):
                        return False, (
                            "{}: iteration {} coefficient count mismatch"
                            .format(filename, iteration)
                        )
                    for index, (test_value, reference_value) in enumerate(
                            zip(test_coefficients,
                                reference_coefficients)):
                        difference = abs(test_value - reference_value)
                        maximum_coefficient = max(
                            maximum_coefficient, difference)
                        if difference > coefficient_tolerance:
                            return False, (
                                "{}: iteration {} coefficient {} "
                                "difference {:.6E} exceeds {:.6E}"
                                .format(filename, iteration, index,
                                        difference,
                                        coefficient_tolerance)
                            )
                    coefficient_l1_difference = abs(
                        test_row["coefficient_l1"] -
                        reference_row["coefficient_l1"])
                    if coefficient_l1_difference > coefficient_tolerance:
                        return False, (
                            "{}: iteration {} coefficient_l1 difference "
                            "{:.6E} exceeds {:.6E}"
                            .format(filename, iteration,
                                    coefficient_l1_difference,
                                    coefficient_tolerance)
                        )
                    niterations += 1

            message = (
                "max gap diff = {:.{p}E} eV, "
                "max eigenvalue-change diff = {:.{p}E} eV, "
                "max residual diff = {:.{p}E} Ha, "
                "max coefficient diff = {:.{p}E} over {} iterations"
            ).format(
                maxima["gap_eV"], maxima["max_delta_eV"],
                max(maxima[name] for name in residual_fields),
                maximum_coefficient, niterations, p=precision)
            return True, message
        except (TypeError, ValueError) as error:
            return False, str(error)

    return inner


def _positive_float(value, label, allow_zero=False):
    result = float(value)
    if not math.isfinite(result) or result < 0.0 or (
            result == 0.0 and not allow_zero):
        relation = "nonnegative" if allow_zero else "positive"
        raise ValueError("{} must be finite and {}".format(label, relation))
    return result


def _paired_file_texts(test_files, reference_files):
    filenames = set(test_files) | set(reference_files)
    if not filenames:
        raise ValueError("no files found")
    pairs = []
    for filename in sorted(filenames, key=str):
        if filename not in test_files or filename not in reference_files:
            raise ValueError("missing file {}".format(filename))
        pairs.append((
            filename,
            _single_text(test_files[filename], filename),
            _single_text(reference_files[filename], filename),
        ))
    return pairs


def _single_text(raw, filename):
    if isinstance(raw, str):
        return raw
    values = list(raw)
    if len(values) != 1 or not isinstance(values[0], str):
        raise ValueError(
            "{} must contain one complete trace".format(filename))
    return values[0]


def _require_header(text, expected, label):
    comments = [
        line.strip() for line in text.splitlines()
        if line.strip().startswith("#")
    ]
    if expected not in comments:
        raise ValueError("{}: missing trace header".format(label))


def _require_same_contract(test_text, reference_text, filename):
    test_contract = _parse_contract(
        test_text, "test {}".format(filename))
    reference_contract = _parse_contract(
        reference_text, "reference {}".format(filename))
    if test_contract != reference_contract:
        differing = sorted(
            key for key in set(test_contract) | set(reference_contract)
            if test_contract.get(key) != reference_contract.get(key)
        )
        raise ValueError(
            "{}: QSGW contract differs for {}"
            .format(filename, differing))
    return test_contract


def _parse_contract(text, label):
    values = {}
    for line_number, line in enumerate(text.splitlines(), 1):
        line = line.strip()
        if not line.startswith("#"):
            continue
        fields = line[1:].strip().split(None, 1)
        if len(fields) != 2 or fields[0] not in CONTRACT_KEYS:
            continue
        key, value = fields[0], fields[1].strip()
        if key in values:
            raise ValueError(
                "{}:{}: duplicate QSGW contract key {}"
                .format(label, line_number, key))
        values[key] = value

    missing = sorted(REQUIRED_CONTRACT_KEYS - set(values))
    if missing:
        raise ValueError(
            "{}: missing QSGW contract keys {}".format(label, missing))
    try:
        version = int(values["qsgw_contract_version"])
    except ValueError as error:
        raise ValueError(
            "{}: invalid QSGW contract version".format(label)) from error
    if version != 6:
        raise ValueError(
            "{}: unsupported QSGW contract version {}"
            .format(label, version))
    values["qsgw_contract_version"] = version

    for key, required in (
        ("fixed_basis", "immutable_mf0"),
        ("live_update", "eigenvalues_wfc"),
        ("wing", "disabled_stage1"),
        ("hartree", "disabled"),
    ):
        if values[key] != required:
            raise ValueError(
                "{}: invalid {} QSGW contract".format(label, key))

    symmetry = values["symmetry"]
    if re.fullmatch(
            r"exx_(?:on|off)_gw_(?:on|off)_rpa_(?:on|off)",
            symmetry) is None:
        raise ValueError(
            "{}: invalid symmetry QSGW contract".format(label))

    head_velocity = {
        "disabled_stage1": "disabled_stage1",
        "scf_grid_analytic_live": "fixed_reference",
    }
    head = values["head"]
    if head not in head_velocity or values["velocity"] != \
            head_velocity[head]:
        raise ValueError(
            "{}: inconsistent head/velocity QSGW contract".format(label))

    band = values["band"]
    if band not in (
            "disabled_stage1", "fixed_reference_rotation_live"):
        raise ValueError(
            "{}: invalid band QSGW contract".format(label))
    cut_keys = (
        "qsgw_band0_unoccupied_keep",
        "qsgw_band0_cut_mode",
        "qsgw_band0_cut_shift_ha",
    )
    cut_contract = values["h_qsgw_cut"]
    if cut_contract == "disabled_non_band":
        if band != "disabled_stage1" or any(
                key in values for key in cut_keys):
            raise ValueError(
                "{}: inconsistent disabled H_QSGW cut contract"
                .format(label))
    elif cut_contract == "band_postprocess":
        if band != "fixed_reference_rotation_live":
            raise ValueError(
                "{}: H_QSGW cut requires band postprocessing"
                .format(label))
        missing_cut = [key for key in cut_keys if key not in values]
        if missing_cut:
            raise ValueError(
                "{}: enabled H_QSGW cut is missing contract fields {}"
                .format(label, missing_cut))
        try:
            unoccupied_keep = int(
                values["qsgw_band0_unoccupied_keep"])
            cut_mode = int(values["qsgw_band0_cut_mode"])
        except ValueError as error:
            raise ValueError(
                "{}: invalid H_QSGW cut integer contract"
                .format(label)) from error
        cut_shift = _finite_float(
            values["qsgw_band0_cut_shift_ha"])
        if unoccupied_keep < 0 or cut_mode not in (0, 1, 2):
            raise ValueError(
                "{}: invalid H_QSGW cut contract".format(label))
        values["qsgw_band0_unoccupied_keep"] = unoccupied_keep
        values["qsgw_band0_cut_mode"] = cut_mode
        values["qsgw_band0_cut_shift_ha"] = cut_shift
    else:
        raise ValueError(
            "{}: invalid H_QSGW cut contract".format(label))

    if not values["qsgw_input_contract"]:
        raise ValueError(
            "{}: empty QSGW input contract path".format(label))
    sha256 = values["qsgw_input_contract_sha256"].lower()
    if re.fullmatch(r"[0-9a-f]{64}", sha256) is None:
        raise ValueError(
            "{}: invalid QSGW input contract SHA256".format(label))
    values["qsgw_input_contract_sha256"] = sha256

    if values["qsgw_mixer"] not in ("none", "linear"):
        raise ValueError(
            "{}: invalid QSGW mixer contract".format(label))
    beta = _finite_float(values["qsgw_mixing_beta"])
    if not 0.0 < beta <= 1.0:
        raise ValueError(
            "{}: invalid QSGW mixing beta".format(label))
    values["qsgw_mixing_beta"] = beta
    return values


def _parse_eigenvalue_trace(text, label):
    _require_header(text, EIGENVALUE_HEADER, label)
    rows = {}
    for line_number, line in enumerate(text.splitlines(), 1):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        fields = line.split()
        if len(fields) != 9:
            raise ValueError(
                "{}:{}: eigenvalue row has {} columns"
                .format(label, line_number, len(fields)))
        try:
            iteration, channel, spin, kpoint = map(int, fields[:4])
            coordinate = tuple(
                _finite_float(value) for value in fields[4:7])
            band = int(fields[7])
            energy = _finite_float(fields[8])
        except ValueError as error:
            raise ValueError(
                "{}:{}: invalid eigenvalue row"
                .format(label, line_number)) from error
        if min(iteration, spin, kpoint, band) < 0 or channel not in (0, 1):
            raise ValueError(
                "{}:{}: invalid eigenvalue index"
                .format(label, line_number))
        key = (iteration, channel, spin, kpoint, band)
        if key in rows:
            raise ValueError(
                "{}:{}: duplicate eigenvalue row {}"
                .format(label, line_number, key))
        rows[key] = (coordinate, energy)
    if not rows:
        raise ValueError("{}: no eigenvalue rows found".format(label))
    return rows


def _validate_eigenvalue_trajectory(rows, contract, label):
    iterations = _continuous_iterations(
        {key[0] for key in rows}, label)
    expected_channels = {0}
    if contract["band"] == "fixed_reference_rotation_live":
        expected_channels.add(1)
    observed_channels = {key[1] for key in rows}
    if observed_channels != expected_channels:
        raise ValueError(
            "{}: eigenvalue channel set differs from the QSGW contract"
            .format(label))

    baseline = {}
    for channel in sorted(expected_channels):
        keys = {
            key[1:] for key in rows
            if key[0] == 0 and key[1] == channel
        }
        if not keys:
            raise ValueError(
                "{}: iteration zero channel {} has no eigenvalues"
                .format(label, channel))
        baseline[channel] = keys

    for iteration in iterations:
        for channel in sorted(expected_channels):
            keys = {
                key[1:] for key in rows
                if key[0] == iteration and key[1] == channel
            }
            if keys != baseline[channel]:
                raise ValueError(
                    "{}: iteration {} channel {} eigenvalue layout "
                    "differs from iteration zero"
                    .format(label, iteration, channel))


def _continuous_iterations(iterations, label):
    observed = sorted(iterations)
    if not observed or observed[0] != 0:
        raise ValueError("{}: iteration zero is missing".format(label))
    if observed[-1] < 1:
        raise ValueError(
            "{}: no completed QSGW iteration is present".format(label))
    expected = list(range(observed[-1] + 1))
    if observed != expected:
        raise ValueError(
            "{}: QSGW iterations are not continuous: {}"
            .format(label, observed))
    return observed


def _parse_iteration_summary(text, label):
    _require_header(text, SUMMARY_HEADER, label)
    rows = {}
    for line_number, line in enumerate(text.splitlines(), 1):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        fields = line.split()
        if len(fields) != 17:
            raise ValueError(
                "{}:{}: summary row has {} columns"
                .format(label, line_number, len(fields)))
        try:
            iteration = int(fields[0])
            row = {
                "max_delta_eV": _finite_float(fields[1]),
                "residual_l2_Ha": _finite_float(fields[2]),
                "residual_max_Ha": _finite_float(fields[3]),
                "efermi_eV": _finite_float(fields[4]),
                "gap_eV": _finite_float(fields[5]),
                "electron_count": _finite_float(fields[6]),
                "requested_mode": int(fields[7]),
                "applied_mode": int(fields[8]),
                "beta": _finite_float(fields[9]),
                "fallback": int(fields[10]),
                "rcond": _finite_float(fields[11]),
                "coefficient_l1": _finite_float(fields[12]),
                "coefficient_count": int(fields[13]),
                "converged": int(fields[14]),
                "coefficients": _parse_coefficients(fields[15]),
                "fallback_reason": fields[16],
            }
        except ValueError as error:
            raise ValueError(
                "{}:{}: invalid summary row"
                .format(label, line_number)) from error
        nonnegative_fields = (
            "max_delta_eV", "residual_l2_Ha", "residual_max_Ha",
            "gap_eV", "electron_count", "beta", "rcond",
            "coefficient_l1",
        )
        if iteration < 0 or any(
                row[name] < 0.0 for name in nonnegative_fields):
            raise ValueError(
                "{}:{}: invalid summary value"
                .format(label, line_number))
        if row["coefficient_count"] != len(row["coefficients"]):
            raise ValueError(
                "{}:{}: coefficient count is inconsistent"
                .format(label, line_number))
        coefficient_l1 = sum(
            abs(value) for value in row["coefficients"])
        if abs(coefficient_l1 - row["coefficient_l1"]) > \
                1.0e-12 * max(1.0, coefficient_l1):
            raise ValueError(
                "{}:{}: coefficient L1 norm is inconsistent"
                .format(label, line_number))
        if iteration in rows:
            raise ValueError(
                "{}:{}: duplicate iteration {}"
                .format(label, line_number, iteration))
        rows[iteration] = row
    if not rows:
        raise ValueError("{}: no iteration rows found".format(label))
    _continuous_iterations(rows, label)
    return rows


def _parse_coefficients(value):
    if value == "none":
        return tuple()
    coefficients = tuple(
        _finite_float(item) for item in value.split(","))
    if not coefficients:
        raise ValueError("empty coefficient list")
    return coefficients


def _finite_float(value):
    result = float(value.replace("D", "E").replace("d", "E"))
    if not math.isfinite(result):
        raise ValueError("non-finite value")
    return result


def _key_mismatch(test, reference):
    test_keys = set(test)
    reference_keys = set(reference)
    if test_keys == reference_keys:
        return None
    missing = sorted(reference_keys - test_keys, key=str)[:3]
    extra = sorted(test_keys - reference_keys, key=str)[:3]
    return "missing={}, extra={}".format(missing, extra)
