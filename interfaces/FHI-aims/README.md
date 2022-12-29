## patch file to interface with FHI-aims

### Usage

At the root of FHI-aims directory, i.e. the path where `src` and `doc` live, run

```bash
$ patch -p1 < fhiaims.hash.patch
```

You can also use git command to apply the modification,
if the source codes are under git control

```bash
$ git apply fhiaims.hash.patch
```

To revert the change, run

```bash
$ patch -R -p1 < fhiaims.hash.patch
```

### Notes on patch files

The hash number indicates the git version of FHI-aims on which the patch file is created.

`fhiaims.fbbdcd540c.patch`: Base on [a commit](https://aims-git.rz-berlin.mpg.de/aims/FHIaims/-/tree/fbbdcd540c775f6feb7e37dba9f73eba67591bdb) on 2022-03-02. Tested, work on carbon diamond. Also applicable on `210716_2` release.

`fhiaims.62c7cfa45f.patch`: Base on [commit 62c7cfa45f](https://aims-git.rz-berlin.mpg.de/aims/FHIaims/-/tree/62c7cfa45f7161f9ac04f5d97a0aaa272af94f34) on 2022-06-03. Used for GW development.

`fhiaims.06f38003bc.patch`: Base on [06f38003bc](https://aims-git.rz-berlin.mpg.de/aims/FHIaims/-/tree/06f38003bccd7601f15ed8e16a218d715a16de85) on 2022-12-09. Add tricoeff derivative dump. Add Cs/dCs dump threshold to reduce the size of Cs/dCs files.

`fhiaims.009e8f35e9.patch`: Base on [009e8f35e9](https://aims-git.rz-berlin.mpg.de/aims/FHIaims/-/tree/009e8f35e99e12dc55e8c1bc222c2c1b593d9fe1) on 2022-12-26. Essentially same as the 06f38003bc patch.
