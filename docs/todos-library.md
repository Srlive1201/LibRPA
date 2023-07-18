
# LibRPA 库

- [ ] librpa.h APIs
  - [ ] dimension 设置: spins, kpoints, bands, basis
  - [ ] MPI communicator
  - [ ] 数组传入
    - [ ] Kohn-Sham energies/eigenvectors
    - [ ] RI coefficients
    - [ ] Coulomb matrices (truncated/untruncated)
    - [ ] Structures: number of atoms
  - [ ] 返回物理性质: RPA correlation energy, self-energy ...
- [ ] driver 和 src 之间充分分离, 主要是 input.cpp input.h
- [ ] 考虑 DFT 程序中数组数据分布和 LiRPA 中相应数组数据分布的差异
  - [ ] ABACUS
