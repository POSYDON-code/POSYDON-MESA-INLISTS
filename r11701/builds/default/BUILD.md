# Build `default` (MESA r11701)

What depends only on the MESA build, not on POSYDON's physics choices:

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

## Not the cluster build

These are not exactly the defaults of the MESA build the v2 grids ran with
on yggdrasil (`mesa-r11701_witheoschange_andwithreverseMTchange_fiximplicitmdot`):
there, `controls.defaults` comments out 10 `xa_*_limit(1)` lines, and
`star_job.defaults` has `remove_center_by_he4` and two
`remove_initial_surface_by_v_surf_div_*` options of the patched MESA. That
build will be added as its own folder (`builds/cluster_patched/`), never as
an edit of this one.
