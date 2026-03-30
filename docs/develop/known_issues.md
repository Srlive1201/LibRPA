# Known Issues

## Caveats

### `nan` in GW output

`nan` values in the `sigc` output indicate that the quasi-particle equation solver was unable to find a stable solution.
This issue is occasionally encountered for unoccupied states far above the Fermi level (typically beyond 20 eV).

When such cases occur only sporadically and not at consecutive k-points, they usually have little effect on the visual appearance of a plotted band structure, and may be ignored for plotting purposes.

If `nan` values occur in larger contiguous regions, the GW calculation should be examined more carefully.
This may indicate that the calculation has become unreliable, for example because the filtering threshold is too large or because the analytic continuation produces multiple nearby solutions.

## Planned Improvements

- [ ] Adapt RPA force work by Mohammad in the [backup branch](https://github.com/Srlive1201/LibRPA/tree/master-backup-240416)
