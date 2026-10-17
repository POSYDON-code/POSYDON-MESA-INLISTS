# MESA installation `default` (r11701)

MESA's own stock files for this installation, independent of POSYDON's
physics choices. They form the `mesa_inlists` layer of the grid setup:

- `defaults/star/controls.defaults`, `defaults/star/star_job.defaults`,
  `defaults/binary/binary_controls.defaults`,
  `defaults/binary/binary_job.defaults`: MESA's namelist defaults. The grid
  setup starts every inlist from these (layer `mesa_inlists`).
- `columns/`: MESA's own `history_columns.list`, `binary_history_columns.list`
  and `profile_columns.list`. POSYDON's grids use the lists in
  `default_common_inlists/` instead; these differ from them.
- `src/`: MESA's stock `binary_run.f`, `run_binary_extras.f` and
  `run_star_extras.f`.

## Origin

Copied byte for byte from POSYDON's `grid_params/defaults/r11701/` at
POSYDON commit `0f699da4f` (branch `enhancement/rewrite_mesa_grids`).

That folder was added in 2022 (`92b23a104`, "code migration") and changed
once since, in 2025 (`c827f2b64`), when a pre-commit hook stripped trailing
whitespace and blank lines at the end of files. Ignoring whitespace, the
files are identical, so the v2 grids, set up in 2023-2024 with
`${posydon_github_root}/grid_params/defaults/r11701/...`, used these
defaults.

## Lines left out on purpose

These files are not a verbatim copy of the `.defaults` of the MESA the v2
grids ran with on yggdrasil
(`mesa-r11701_witheoschange_andwithreverseMTchange_fiximplicitmdot`). That
MESA's files have 15 more settings, which are not here on purpose:

- `star/controls.defaults`: the 12 lines
  `xa_{central,surface,average}_{lower,upper}_limit_species(1) = ''` and
  `xa_{central,surface,average}_{lower,upper}_limit(1) = 0`;
- `star/star_job.defaults`: `remove_center_by_he4 = 0`,
  `remove_initial_surface_by_v_surf_div_cs = 0` and
  `remove_initial_surface_by_v_surf_div_v_escape = 0`.

The grid setup starts every inlist from these files, and MESA's original
lines for these settings did not work in an inlist with that version of
MESA. Leaving them out changes nothing: their values are MESA's own
defaults (off), which MESA uses whenever a setting is absent from an
inlist. Keep them out when updating these files.
