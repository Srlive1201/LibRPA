## Regression tests

### Create a regression test case data

**Step 1**: create reference and testcase folder for the case
```bash
export CASE=case_new
mkdir -p {refs,testcases}/${CASE}/{dataset,librpa}
```

Step 2: prepare dataset
```bash
cd refs/${CASE}
cd dataset
# Create inputs and run dataset preparation program.
# For example, FHI-aims:
#     vim control.in; vim geometry.in; aims.x > aims.out
cd ..
# Compress dataset as dataset.tar.gz. Exclude all outputs, but keep all input files.
# For example, FHI-aims
tar --exclude aims.out -zchf dataset.tar.gz dataset
```

Step 3: prepare LibRPA input and run LibRPA to get reference data
```bash
mkdir librpa
cat > librpa.in << EOF
task = rpa
input_dir = ../dataset
# other parameters
# ...
EOF
mpirun -np 4 chi0_main > librpa.out
cd ..
```

Step 4: archive output to reference folder
```bash
cp -a dataset/aims.out ../../refs/${CASE}/dataset/aims.out.export
cp -a librpa/librpa.out ../../refs/${CASE}/librpa/librpa_omp1_np4.out
```
Addtional reference data (band structures, benchmark from other programs)
can be added to the reference `dataset` folder.
