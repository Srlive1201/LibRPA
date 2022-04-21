## patch file to interface with FHI-aims

### Usage

At the root of FHI-aims directory, i.e. the path where `src` and `doc` live, run

```bash
patch -p1 < fhiaims.hash.patch
```

You can also use git command to apply the modification,
if the source codes are under git control

```bash
git apply fhiaims.hash.patch
```

To revert the change, run

```bash
patch -R -p1 < fhiaims.hash.patch
```

### Notes on patch files

The hash number indicates the git version of FHI-aims on which the patch file is created.

`fhiaims.fbbdcd540c.patch`: Base on [a commit](https://aims-git.rz-berlin.mpg.de/aims/FHIaims/-/tree/fbbdcd540c775f6feb7e37dba9f73eba67591bdb) on 2022-03-02. Tested, work on carbon diamond. Also applicable on `210716_2` release.

