# Cluster dispatch usage

Launch HTCondor jobs for InnOpTEM over instances and settings.

## 1) Prepare settings
1. Create a folder under `settings/`, e.g. `settings/myfolder/`.
2. Add one or more `.toml` files. Each file defines a solver configuration.

## 2) Prepare instances
1. Put test instance files in `check/instances/main/`.
   - Supports root-level `.m` files and subfolders with `*.m` case files.
   - Example path: `check/instances/main/pglib_scopf_case323.m`

## 3) Run dispatch script
```
./check/scripts/cluster/dispatch_testset.sh --settings settings/myfolder
```

Or run a single settings file:
```
./check/scripts/cluster/dispatch_testset.sh --settings settings/myfolder/single_setting.toml
```

### Optional arguments
- `--experiment_name <name>`: (default `experiment`) output subfolder under `check/results/`
- `--min_buses <n>`: only include instances with at least n buses (default `73`)
- `--max_buses <n>`: only include instances with at most n buses (default `5000`)
- `--dry-run`: show planned jobs without submitting to condor

## 4) Output structure
- Results: `check/results/<experiment_name>/<setting>/<case>.<setting>.csv`
- Condor logs: `check/results/cluster_logs/<experiment_name>_<case>.<setting>.*`
- temp .sub files: `check/results/temp/`

## 5) Behavior and policies
- Jobs are skipped if results already exist or are already queued.
- A periodic job cleanup policy is configured in the generated submit file:
  - `periodic_remove = (JobStatus == 2 && (CurrentTime - JobStartDate) > 14400)` (remove running jobs after 4 hours)

## 6) Visualization
After job completion run:
```
python3 check/scripts/visu/compare.py
```
This reads `check/results/<experiment_name>/` and produces summary tables and performance/profile plots.
