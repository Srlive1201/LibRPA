#!/usr/bin/env python3
"""Convert native LibRPA model arrays to TRIQS; perform no projection or continuation."""
from pathlib import Path
import argparse
import json
import numpy as np


def load_model(path):
    path = Path(path)
    with np.load(path, allow_pickle=False) as archive:
        model = dict(archive)
    metadata = json.loads(model['metadata_json'].tobytes())
    if metadata['format_version'] != 1:
        raise ValueError('unsupported LibRPA model format')
    with np.load(path.parent / metadata['static_file'], allow_pickle=False) as archive:
        model.update(dict(archive))
    model['metadata'] = metadata
    return model


def triqs_data(model):
    """Build a band-basis dft_input dictionary, retaining every AFM site separately."""
    metadata = model['metadata']
    scale = metadata['ha_to_ev']
    ns, nk, nb = model['E'].shape
    so = int(metadata['n_spinor'] == 2)
    sites = model['output_sites']
    labels = json.loads(model['output_labels_json'].tobytes())
    shells = []
    for index, (atom, start, dim) in enumerate(sites):
        angular = {label.split(':')[0][0] for label in labels[start:start+dim]}
        if len(angular) != 1 or next(iter(angular)) not in 'spd':
            raise ValueError('mixed-angular-momentum sp sites require a custom DMFT shell definition')
        l = 'spd'.index(next(iter(angular)))
        shells.append(dict(atom=int(atom), sort=index, l=l, dim=int(dim), SO=so, irep=0))
    proj = np.zeros((nk, ns, len(sites), int(max(sites[:,2])), nb), complex)
    for index, (_, start, dim) in enumerate(sites):
        proj[:,:,index,:dim,:] = model['output_P'][:,:,start:start+dim,:].transpose(1,0,2,3)
    hopping = np.zeros((nk,ns,nb,nb), complex)
    bands = np.arange(nb)
    hopping[:,:,bands,bands] = model['E'].transpose(1,0,2)*scale
    result = dict(energy_unit=scale, dft_code='hk', n_k=nk, k_dep_projection=0,
        SP=int(ns == 2 or so), SO=so, charge_below=0.,
        density_required=float(np.einsum('k,skb->',model['k_weights'],model['occupations'])),
        symm_op=0, n_shells=len(shells), shells=shells, n_corr_shells=len(shells),
        corr_shells=shells, n_inequiv_shells=len(shells), corr_to_inequiv=list(range(len(shells))),
        inequiv_to_corr=list(range(len(shells))), use_rotations=0,
        rot_mat=[np.eye(s['dim'],dtype=complex) for s in shells], rot_mat_time_inv=[0]*len(shells),
        n_reps=[1]*len(shells), dim_reps=[[s['dim']] for s in shells],
        T=[np.eye((2*s['l']+1)*(so+1),dtype=complex) for s in shells],
        n_orbitals=np.full((nk,ns),nb,dtype=int), proj_mat=proj,
        bz_weights=model['k_weights'], kpt_weights=model['k_weights'], kpts=model['k_fractional'],
        hopping=hopping)
    # Full native tensors and their separate spin-frame conventions stay explicit.
    interaction = {k:v*scale if k.startswith(('v_','U_','continued_U_','continued_spread_')) else v
                   for k,v in model.items() if k.startswith(('v_','U_','continued_'))}
    interaction['nu'] = model['nu']*scale
    if 'continued_z' in interaction:
        interaction['continued_z'] = model['continued_z']*scale
    interaction['metadata_json'] = model['metadata_json']
    interaction['tensor_sites'] = model['tensor_sites']
    if 'A_weights' in model:
        interaction['A_weights'] = model['A_weights']
    interaction['energy_unit'] = 'eV'
    interaction['frequency_unit'] = 'eV'
    return result, interaction


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('model',type=Path)
    parser.add_argument('output',type=Path)
    args=parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    data,interaction=triqs_data(load_model(args.model))
    # Do not silently feed unstable extrapolations to an impurity solver.
    has_static = 'continued_z' in interaction and interaction['continued_z'][0] == 0
    for key,value in interaction.items():
        if has_static and key.startswith('continued_status_') and np.any(value[0]):
            raise ValueError('static continuation has flagged entries; inspect native diagnostics first')
    from h5 import HDFArchive
    with HDFArchive(str(args.output),'w') as archive:
        archive['dft_input']=data
        archive['crpa_interaction']=interaction
    print(args.output)


if __name__=='__main__':
    main()
