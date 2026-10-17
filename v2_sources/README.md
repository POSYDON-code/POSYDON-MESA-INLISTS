# v2 sources: the inlists behind the POSYDON v2 production grids

Starting material for restructuring this repository. Every POSYDON v2 grid run
on yggdrasil names a scenario commit of this repo; this folder holds the
physics inlists of each of those commits side by side, on one branch, so none
of them depends on a feature branch surviving.

Nothing here is the final layout yet. The next step is to pick a base, diff
the others against it, and turn the differences into release and rerun
folders.

## Contents

- `inlists/<short-hash>_<branch>/`: `default_common_inlists/`, `CO-H_star/`,
  `CO-He_star/` and `single_He_star/` exactly as in that commit (`r11701/`
  dropped from the path).
  `35dd7a4_jeff_new_grid/` also holds `HeMS-HMS/`; no v2 grid used it, but
  HeMS-HMS is a supported grid type.
- `zams_models/`: the ZAMS models used by all these commits. The commits
  hold one of two versions of `r11701/ZAMS_models/`, but the newer one is a
  superset: the 6 shared files are byte-identical, and it adds the 0.2, 0.45
  (2023-03-07) and 2 Z_sun (2023-04-05) models. Only the newer one is kept.
- `builds/default/`: what depends only on the MESA build: MESA's namelist
  `defaults/`, its stock column lists (`columns/`) and Fortran sources
  (`src/`). Byte-identical copy of POSYDON's `grid_params/defaults/r11701/`,
  which the v2 ini files read through `${posydon_github_root}`; see
  `builds/default/BUILD.md` for its history and how it differs from the
  patched MESA build on yggdrasil.
- `templates/`: the annotated ini templates from POSYDON's `grid_params/`:
  `grid_params.ini` (grid setup) and `pipeline_{yggdrasil,quest}.ini`
  (post-processing). Byte-identical copies; their paths still point into a
  POSYDON clone (`${posydon_github_root}/grid_params/...`) until the setup
  reads this repository directly.
- `production_grids.csv`: one row per ini file found in a v2 run folder
  (797 rows), with columns:
  - `grid_type`, `metallicity`, `run`, `ini_file`: where the ini sits
  - `cluster_folder`: `running_grids` or `data_sync_staging`, both under
    `/srv/astro/groups/posydon/`
  - `scenario_branch`, `scenario_commit`, `system_type`, `single_star_grid`,
    `grid_csv`: as written in the ini
  - `inlist_folder`, `zams_models`: the folders here that hold its inputs
  - `replay`: today's `posydon-setup-grid` re-run on that ini, compared key
    by key with the inlists the run used:
    `identical`, `hand-edited` (inlists changed after setup),
    `setup fails today`, or `not replayable` (no inlists or grid copied, or
    the ini does not belong to the run)
  - `notes`: oddities, e.g. the non-existent system type `MS-MS`
- Not included: POSYDON's submodule `grid_params/POSYDON-MESA-INLISTS`
  (this repository, pinned at `c1f75ac`).
- Not included: `r11701/running_scripts/` (64-102 MB of parameter files;
  4 versions). Read it from the commit if needed.

## Commits used

Sorted by number of ini files. "On main" means reachable from `origin/main`;
"Branches" counts remote branches containing the commit (1 = it disappears
if that branch is deleted).

| Folder | Branch | Commit | Date | Ini files | On main | Branches |
| --- | --- | --- | --- | --- | --- | --- |
| `0e736a8_fix_iniHe_pre_rerun5` | `fix_iniHe_pre_rerun5` | `0e736a853c3c2e522b76299a54aea8f5eea15d08` | 2023-12-18 | 127 | no | 33 |
| `0513438_matthias_LBV_wind` | `matthias_LBV_wind` | `051343856f50dc5518abdf28904d62df4caacaff` | 2024-07-24 | 104 | no | 31 |
| `dcc42e1_matthias_lbv_thermo_dedt_noAgeLimit` | `matthias_lbv_thermo_dedt_noAgeLimit` | `dcc42e1d6d89fc267826970643e7f4f1bf463171` | 2024-08-29 | 99 | no | 6 |
| `2737b3b_seth_lbv_thermo_dedt_hepulse` | `seth_lbv_thermo_dedt_hepulse` | `2737b3bce69066751477c802a6e515183e13d303` | 2024-09-20 | 92 | no | 1 |
| `18710e9_main` | `main` | `18710e943edd926a5653c4cdb7d6e18e5bdb35a2` | 2022-12-23 | 86 | yes | 77 |
| `a119ea9_matthias_no_age_limit` | `matthias_no_age_limit` | `a119ea9548cca9b6b7ddfe31309ced1ae40c271e` | 2024-08-21 | 80 | no | 6 |
| `936f76f_seth_lbv_thermo_dedt` | `seth_lbv_thermo_dedt` | `936f76f82d44b16139ff89e41e5b0f2232fd9b30` | 2024-08-14 | 55 | no | 11 |
| `7db33d9_matthias_fix_initial_Hestar` | `matthias_fix_initial_Hestar` | `7db33d9fe37af4b1d391e74209051ce334eb6b32` | 2023-12-15 | 42 | no | 33 |
| `a060b7b_main` | `main` | `a060b7bf1a4d94d693f77be9a1b0b3a522c1eadf` | 2023-01-11 | 21 | yes | 77 |
| `f1803af_seth_dedt_energy_eqn` | `seth_dedt_energy_eqn` | `f1803afe400ed16a4fad47b0bad6abbf0965cae1` | 2024-06-05 | 16 | no | 2 |
| `ecae5ee_seth_rerun_egrav_Pgas_dXdt_terms` | `seth_rerun_egrav_Pgas_dXdt_terms` | `ecae5eecfc3cecfb20f50b673d0f575bc5a75998` | 2024-03-21 | 16 | no | 1 |
| `0274be9_seth_dedt_hepulse` | `seth_dedt_hepulse` | `0274be9895480cae3aa2a4d14b056a2447d5921b` | 2024-06-06 | 16 | no | 1 |
| `c4f90ce_development` | `development` | `c4f90ce2d93595f66751011d2002fc3ab3d090ec` | 2023-09-21 | 15 | no | 49 |
| `afa1860_zepei_fix_implicit` | `zepei_fix_implicit` | `afa1860ddf9894aa1d82742ee2a73e8e92acd4a9` | 2023-04-06 | 9 | no | 55 |
| `91ee7ca_seth_TPAGBwind_runstar` | `seth_TPAGBwind_runstar` | `91ee7cac7111b19fc7e10e20a47f9cfec320ba6b` | 2023-06-22 | 9 | no | 54 |
| `d682283_matthias_PISN` | `matthias_PISN` | `d68228338b91fd487ef5d55c9b6ebb8cc5f0e668` | 2023-03-08 | 6 | no | 55 |
| `77640b5_seth_lbv_thermo_dedt_hepulse` | `seth_lbv_thermo_dedt_hepulse` | `77640b5f773825129d6e3bf29267e52c90f328a2` | 2024-09-19 | 4 | no | 1 |

Generated on 2026-10-16 from a read-only copy of the yggdrasil run folders,
with the replay in POSYDON's `grid_rewrite_tests/` (branch
`enhancement/rewrite_mesa_grids`).
