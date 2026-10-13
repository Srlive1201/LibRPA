#!/usr/bin/env python3
"""Generate deterministic reader-v1 Delta-ST matrices for this regression."""

from pathlib import Path
import struct
import sys


MARKER = -41073291
ATOM_NAUX = (17, 17)
N_AUX = sum(ATOM_NAUX)
FREQUENCIES = ((1, 0.5, 0.4), (2, 1.5, 0.6))
SCALES = {(1, 1): 0.010, (1, 2): 0.006, (2, 1): 0.008, (2, 2): 0.004}


def matrix_value(row: int, column: int, scale: float) -> complex:
    if row == column:
        return complex(-scale, 0.0)
    if abs(row - column) == 1:
        return complex(-0.1 * scale, 0.0)
    return 0.0j


def write_response(
    path: Path, iq: int, ifreq: int, omega: float, weight: float
) -> None:
    pairs = ((0, 0), (0, 1), (1, 1))
    header_bytes = 7 * 4 + 2 * 8 + len(ATOM_NAUX) * 4 + len(pairs) * (4 + 8)
    offsets = []
    next_offset = header_bytes
    for atom_i, atom_j in pairs:
        offsets.append(next_offset)
        next_offset += ATOM_NAUX[atom_i] * ATOM_NAUX[atom_j] * 16

    scale = SCALES[(iq, ifreq)]
    atom_offsets = (0, ATOM_NAUX[0])
    with path.open("wb") as output:
        output.write(struct.pack("<6i2di", MARKER, iq, ifreq, N_AUX, 1, 2,
                                 omega, weight, len(pairs)))
        output.write(struct.pack("<2i", *ATOM_NAUX))
        for pair_index, offset in enumerate(offsets):
            output.write(struct.pack("<iq", pair_index, offset))
        for atom_i, atom_j in pairs:
            for local_i in range(ATOM_NAUX[atom_i]):
                row = atom_offsets[atom_i] + local_i
                for local_j in range(ATOM_NAUX[atom_j]):
                    column = atom_offsets[atom_j] + local_j
                    value = matrix_value(row, column, scale)
                    output.write(struct.pack("<2d", value.real, value.imag))


def main() -> None:
    if len(sys.argv) > 2:
        raise SystemExit("usage: generate_sternheimer_fixture.py [output-directory]")
    fixture_dir = (
        Path(sys.argv[1]).resolve()
        if len(sys.argv) == 2
        else Path(__file__).resolve().parent / "librpa" / "fixtures"
    )
    fixture_dir.mkdir(parents=True, exist_ok=True)
    (fixture_dir / "v1_sternheimer_qpoints.dat").write_text(
        "1 0.0 0.0 0.0 0.25\n2 0.5 0.5 0.5 0.75\n", encoding="ascii"
    )
    for iq in (1, 2):
        for ifreq, omega, weight in FREQUENCIES:
            write_response(
                fixture_dir / f"v1_sternheimer_chi0_iq_{iq}_ifreq_{ifreq}.bin",
                iq,
                ifreq,
                omega,
                weight,
            )


if __name__ == "__main__":
    main()
