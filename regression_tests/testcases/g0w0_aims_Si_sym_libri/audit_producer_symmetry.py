from pathlib import Path
import json
import sys
import numpy as np

root = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parent / 'dataset'
lines = (root / 'band_out').read_text().splitlines()
nk, ns, nb, na = [int(lines[i]) for i in range(4)]
energies = np.array([[float(row.split()[3]) for row in lines[6+i*(nb+1):6+i*(nb+1)+nb]] for i in range(nk)])
k = np.array([[float(x) for x in row.split()[2:5]] for row in (root / 'bz_sampling_out').read_text().splitlines() if len(row.split()) == 10])
stru = (root / 'stru_out').read_text().splitlines()
ops = [np.array([int(x) for x in row.split()[:9]]).reshape(3,3) for row in stru[10:58]]
assert len(k) == 27 and len(ops) == 48
max_error = 0.; location = None
for isym, r in enumerate(ops, 1):
    for ik, ki in enumerate(k):
        target = ki @ np.linalg.inv(r)
        residual = k-target
        residual -= np.rint(residual)
        jk = int(np.argmin(np.linalg.norm(residual, axis=1)))
        assert np.linalg.norm(residual[jk]) < 1e-9
        diff = np.abs(energies[ik]-energies[jk])
        ib = int(np.argmax(diff))
        if diff[ib] > max_error:
            max_error = float(diff[ib]); location = {'operation':isym,'k1':ik+1,'k2':jk+1,'state':ib+1}
result = {'max_input_KS_star_spread_eV':max_error, 'location':location,
          'input_Gamma_valence_triplet_12_to_14_spread_eV':float(np.ptp(energies[0,11:14])),
          'method':'Compare sorted producer eigenvalues at symmetry-related full-grid k points using stru_out col operations.',
          'meaning':'The archived FHI-aims fixture itself is not exactly symmetry covariant. This measures its spectral discrepancy, not a proof that every GW residual is explained.'}
print(json.dumps(result, indent=2))
