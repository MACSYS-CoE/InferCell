# Talks

Source files for talks given about this project. Each talk lives in its
own date-prefixed subfolder (e.g. `2026-04-macsys/`). Shared assets —
the beamer theme (`beamerthemeoxfordmaths.sty`), reusable `figures/`,
and the `posterior_plot/` title-logo source — sit at the root of this
directory and are referenced from each talk via relative paths.

To compile a talk:

```
cd <talk-folder>
latexmk -pdf talk.tex
```
