# TODO

Tracked improvements and technical debt.

## Ebi binary used instead of ebi Python package

A fairly broad Python wrapper is now available for most ebi features. The 
current requirement for ebi to be on the PATH or to override 
`probabilities.EBI_EXECUTABLE` can be obsoleted if this package will suffice.

This would also fix the CLI-shape fragility this code currently has: the
`disc occ`/`prob trac` subprocess calls in `probabilities.py` are written
against a specific Ebi CLI version's argument shape (already had to be
patched once — see CHANGELOG — after Ebi's `disc occ` gained a required
`stochastic-labelled-Petri-net` subcommand token), with no version pin or
check anywhere. A Python dependency would put Ebi under normal dependency
management instead.

