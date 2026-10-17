# POSYDON-MESA-INLISTS

The MESA inlists behind POSYDON's grids, versioned per MESA version and per
POSYDON release.

Note: for MESA version 11701, we use mesa_sdk: 20190503

## Layout

```
r11701/                       one folder per MESA version
  mesa/default/               MESA's own stock files for one installation:
                              defaults/, columns/ and src/ (the mesa_inlists
                              layer); mesa/cluster_patched/ later
  ZAMS_models/                ZAMS models for grids that start from ZAMS
  releases/
    v2.0/                     the inlists of POSYDON v2.0 (frozen)
      release.ini             MESA installation, frozen flag, source commit
      common/                 shared by every grid type: binary/inlist1,
                              inlist2, inlist_project, src/, column lists
      HMS-HMS/                one folder per grid type, each with a
      CO-HMS_RLO/             scenario.ini manifest and the files the grid
      CO-HeMS/                type adds on top of common/
      HeMS-HMS/
      single_HMS/
      single_HeMS/
      reruns/<name>/          physics variants of v2.0 (overlays)
    development/              the next release; new work lands here
  running_scripts/            the v2 parameter files and example ini files
templates/                    annotated ini templates (grid setup, pipeline)
v2_sources/                   archive: the inlists of every v2 production
                              commit and production_grids.csv
```

### Grid types and `scenario.ini`

Each grid-type folder has a `scenario.ini` that lists, per target, the files
the grid type uses, in order (a later file wins key by key):

- `[grid_type]`: name, description and how each star starts (`zams`,
  `formation`, `point_mass` or `none`);
- `[inlists]`: `binary`, `star1`, `star2` (both `&controls` and `&star_job`),
  `star1_formation_every_step` and `star1_formation_step<N>`;
- `[columns]`: the history, binary-history and profile column lists;
- `[extras]`: the Fortran sources, by their legacy `[mesa_extras]` key;
- `[zams]`: a ZAMS model the grid type fixes, if any.

Paths are relative to the release folder. A grid type that uses common
inlists only (HMS-HMS, single_HMS) has a folder with just its manifest; the
single-star grid types add `single_star_flag` (`x_logical_ctrl(1) = .true.`,
read by `run_star_extras.f`).

### Reruns

A rerun is a variant of a release that changes some physics. Its folder
mirrors the release and holds only what it changes; its `scenario.ini`
says how:

- `[overlays]`: an inlist with only the changed keys, applied right after
  the release file of the same path, wherever that file is used;
- `[replace]`: a whole file (a Fortran source) used instead of the
  release's.

`[rerun]` records the source branch and commit the rerun came from and the
grid types the v2 production grids ran it for.

`v2.0` is `main` at `9fca459`, identical to production commit `a060b7b`.
Its 17 reruns are the other 16 production commits of the v2 grids (see
`v2_sources/README.md`), plus `jeff_new_grid`, which no v2 grid used but
which holds HeMS-HMS as it was set up. Every one of these commits is
reproduced exactly by `v2.0` plus its rerun, for every grid type: each
target's merged keys, the column lists and the Fortran sources are the same
as those today's `posydon-setup-grid` takes from the commit itself.

### Rules

- A release folder is frozen once its `release.ini` says `frozen = true`;
  change `development/` instead.
- A new MESA installation (e.g. a patched one) is a new folder under
  `mesa/`, never an edit.
- A new MESA version is a sibling of `r11701/` with the same layout.

## Old refs

Grids set up before this layout name a commit of this repository in their
ini file (`scenario = ['posydon', '<branch>-<commit>', '<system type>']`).
Those commits still exist and still check out their own `r11701/` folders,
so they keep working with the setup code of their time.
