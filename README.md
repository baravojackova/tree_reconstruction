# Tree QSM comparison

Reconstruction of trees from TLS point clouds using different QSM methods, and
export of the resulting cylinder geometry to **ANSYS** for structural analysis.

The aim is to compare reconstruction methods (TreeQSM, AdQSM, AdTree, TreeGraph)
against destructively measured reference volumes, and to decide which output is
best suited as a beam model for FE analysis.

## Folder structure

```
tree-qsm-comparison/
├── scripts/
│   ├── python/          # conversion, extraction and comparison scripts
│   ├── matlab/          # run_treeqsm.m (TreeQSM reconstruction)
│   └── ansys/           # buk.mac (APDL beam model)
├── results/
│   └── volume_results.csv   # shared master results table (versioned)
├── data/                # point clouds, .ply skeletons (NOT versioned), data/<series>/<tree>/
└── output/              # geom_*.txt, ANSYS files (NOT versioned)
```

Only the scripts and `results/volume_results.csv` are tracked by git; data and
outputs are excluded via `.gitignore` because of their size.

## Scripts

| Script | Purpose |
|---|---|
| `scripts/matlab/run_treeqsm.m` | TreeQSM reconstruction; writes `volumes_<tree>_<run>.csv` and the ANSYS geometry export |
| `scripts/python/ply_to_geom.py` | Converts an AdTree skeleton `.ply` into an ANSYS `geom.txt` beam model (prune, smooth, adaptive resample, optional AdQSM radius calibration, 3D plot) |
| `scripts/python/reference_volume.py` | Extracts the destructive reference volume of one tree from the field measurement table |
| `scripts/python/qsm_volume_mean.py` | Mean/std volume from the 20 published TreeQSM realization files of one tree |
| `scripts/python/import_matlab_results.py` | Imports my own TreeQSM `volumes_*.csv` tables (several trees/runs at once) |
| `scripts/python/compare_volumes.py` | Prints all methods side by side per tree and computes Bias / MAE / RMSE / CV-RMSE against the reference |

## How the results are collected

Every script writes its result into one shared table,
`results/volume_results.csv`, with the columns:

```
tree, method, total_m3, stem_m3, branch_m3, std_m3
```

Writing uses an *upsert*: a row with the same `tree` + `method` is replaced, so
scripts can be re-run without duplicating rows. `compare_volumes.py` then reads
this table, so adding a new tree or a new method means adding rows — nothing
else changes.

## Workflow

1. Reconstruct the tree (TreeQSM in MATLAB, AdTree/AdQSM on Windows).
2. Run the extraction scripts — each appends its row to `volume_results.csv`.
3. Run `compare_volumes.py` to compare the methods against the reference.
4. Use `ply_to_geom.py` (or the MATLAB export) to produce the ANSYS geometry,
   then run `buk.mac`.

## Notes

- Volumes are in **m³** everywhere in the shared table. TreeQSM stores volumes
  in litres internally, so the MATLAB script divides by 1000.
- Tree IDs must match exactly across all scripts (e.g. `IND01_054`), otherwise
  the rows will not line up in the comparison.
- AdTree skeleton radii are not calibrated and overestimate branch thickness;
  `ply_to_geom.py` can rescale them using AdQSM taper and per-order medians.
- Reference volumes (de Tanago dataset) cover the stem and branches down to
  ~10 cm diameter, so finer branches are not included in the reference.

## Paths on a new computer (`local_paths.json`)

Project folders (`data/`, `scripts/ansys`, `scripts/matlab`, ...) are worked out
automatically from where `scripts/python/paths.py` lives, so nothing has to be
edited when the project is moved between the PC and the Mac. Only folders that
live **outside** the project (scan folders, TreeGraph results, a separate MATLAB
data folder) differ between computers. They are stored in `local_paths.json`
in the project root, which is ignored by Git (each computer has its own copy).

To create it on a new computer:

1. Copy `local_paths.example.json` to `local_paths.json`.
2. Replace every `SEM_DOPLN_CESTU/...` with a real, existing folder. Use `/` in
   paths (also on Windows). You can delete keys you do not need
   (`matlab_data_dir` is optional: without it `runsken.m` uses its own folder).
3. If a key points to a folder that does not exist, the scripts stop with a
   clear error asking you to fix the path or delete the key.

## Dependencies

Python: `numpy`, `scipy`, `matplotlib` (`pip install numpy scipy matplotlib`).
MATLAB R2025b with TreeQSM. ANSYS Mechanical APDL.
