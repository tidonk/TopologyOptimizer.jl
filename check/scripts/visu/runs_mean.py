from pathlib import Path
import os
import pandas as pd

from across_all import _stage_name_from_path, read_experiment
from tables import TIME_LIMIT


SHARED = Path(os.environ.get("TOPOLOGYOPTIMIZER_EXPERIMENT_ROOT", Path(__file__).resolve().parents[4]))       # root of the experiment directories
R = SHARED

### SP solving ###
#SOURCES_FULL = [R / "2-SP_solving/check/results/2-SP_solving1",
#                R / "2-SP_solving/check/results/2-SP_solving2",
#                R / "2-SP_solving/check/results/2-SP_solving3"]

#OUTPUT_ROOT = R / "2-SP_solving/check/results/mean"

### SP selection ###
#SOURCES_FULL = [R / "3-SP_selection/check/results/2-score-diversity-full1",
#                R / "3-SP_selection/check/results/2-score-diversity-full2",
#                R / "3-SP_selection/check/results/2-score-diversity-full3",]

#SOURCES_FULL = [R / "3-SP_selection/check/results/2-score-diversity-full-hybrid1",
#                R / "3-SP_selection/check/results/2-score-diversity-full-hybrid2",
#                R / "3-SP_selection/check/results/2-score-diversity-full-hybrid3",]

#SOURCES_FULL = [R / "3-SP_selection/check/results/2-score-diversity-full-medoid1",
#                R / "3-SP_selection/check/results/2-score-diversity-full-medoid2",
#                R / "3-SP_selection/check/results/2-score-diversity-full-medoid3",]

#SOURCES_FULL = [R / "3-SP_selection/check/results/7-score-full-eval1",
#                R / "3-SP_selection/check/results/7-score-full-eval2",
#                R / "3-SP_selection/check/results/7-score-full-eval3"]

#SOURCES_FULL = [R / "3-SP_selection/check/results/8-pb-full1",
#                R / "3-SP_selection/check/results/8-pb-full2",
#                R / "3-SP_selection/check/results/8-pb-full3",]

#SOURCES_FULL = [R / "3-SP_selection/check/results/9-full-eval1",
#                R / "3-SP_selection/check/results/9-full-eval2",
#                R / "3-SP_selection/check/results/9-full-eval3"]

#OUTPUT_ROOT = R / "3-SP_selection/check/results/mean"

### Combined SP + CF ###
#SOURCES_FULL = [R / "4-cut_filtering/check/results/9-cutfiltering-full1",
#                R / "4-cut_filtering/check/results/9-cutfiltering-full2",
#                R / "4-cut_filtering/check/results/9-cutfiltering-full3"]

#SOURCES_FULL = [R / "4-cut_filtering/check/results/9-cutfiltering-full-spdiv1",
#                R / "4-cut_filtering/check/results/9-cutfiltering-full-spdiv2",
#                R / "4-cut_filtering/check/results/9-cutfiltering-full-spdiv3"]

#OUTPUT_ROOT = R / "4-cut_filtering/check/results/mean"

### Partial decomp ###
#SOURCES_FULL = [R / "5-partial_decomp/check/results/2-full-eval_num_11",
#                R / "5-partial_decomp/check/results/2-full-eval_num_12",
#                R / "5-partial_decomp/check/results/2-full-eval_num_13"]

#SOURCES_FULL = [R / "5-partial_decomp/check/results/2-full-eval_num_1_div_full1",
#                R / "5-partial_decomp/check/results/2-full-eval_num_1_div_full2",
#                R / "5-partial_decomp/check/results/2-full-eval_num_1_div_full3"]

SOURCES_FULL = [R / "5-partial_decomp/check/results/2-full-eval_num_1_no_base_div1",
                R / "5-partial_decomp/check/results/2-full-eval_num_1_no_base_div2",
                R / "5-partial_decomp/check/results/2-full-eval_num_1_no_base_div3"]

#SOURCES_FULL = [R / "5-partial_decomp/check/results/2-full-eval-scores1",
#                R / "5-partial_decomp/check/results/2-full-eval-scores2",
#                R / "5-partial_decomp/check/results/2-full-eval-scores3"]

#SOURCES_FULL = [R / "5-partial_decomp/check/results/2-full-eval-scores_div_full1",
#                R / "5-partial_decomp/check/results/2-full-eval-scores_div_full2",
#                R / "5-partial_decomp/check/results/2-full-eval-scores_div_full3"]

OUTPUT_ROOT = R / "5-partial_decomp/check/results/mean"

### Fractional Sepa ###

#SOURCES_FULL = [R / "82-fractional_sepa/check/results/2-fractional_sepa-full_div1",
#                R / "82-fractional_sepa/check/results/2-fractional_sepa-full_div2",
#                R / "82-fractional_sepa/check/results/2-fractional_sepa-full_div3"]

#OUTPUT_ROOT = R / "82-fractional_sepa/check/results/mean"

### Branching ###

#SOURCES_FULL = [R / "83-branching/check/results/2-bprio-full1",
#                R / "83-branching/check/results/2-bprio-full2",
#                R / "83-branching/check/results/2-bprio-full3"]

#SOURCES_FULL = [R / "83-branching/check/results/2-bprio-full_div1",
#                R / "83-branching/check/results/2-bprio-full_div2",
#                R / "83-branching/check/results/2-bprio-full_div3"]

#OUTPUT_ROOT = R / "83-branching/check/results/mean"

### Nested ###

#SOURCES_FULL = [R / "9-nested/check/results/3-nested-bpfi-full1",
#                R / "9-nested/check/results/3-nested-bpfi-full2",
#                R / "9-nested/check/results/3-nested-bpfi-full3"]

#OUTPUT_ROOT = R / "9-nested/check/results/mean"

# Column in the experiment DataFrames to use as the performance metric.
TARGET = "time_solve"

# copied from accross_all, we filter out all infeasible instances first
def _mean_per_instance_full(df: pd.DataFrame) -> pd.DataFrame:
    """Aggregate multiple runs per ``(instance, setting [, tal])`` preserving
    the full CSV schema.

    Numeric columns are arithmetic-mean averaged across runs.  Non-numeric
    columns take the first value (assumed identical across runs for the same
    instance/setting).  ``time_solve`` (TARGET) is penalised to ``TIME_LIMIT``
    for any non-OPTIMAL row before averaging.  The aggregated ``status`` is
    ``"OPTIMAL"`` only if every constituent row was OPTIMAL, else
    ``"TIME_LIMIT"`` (matches ``_mean_per_instance``).
    """
    if df.empty:
        return df.copy()

    df = df.copy()
    df = df[df["status"] != "INFEASIBLE"]

    # Penalise non-OPTIMAL TARGET values before averaging
    if TARGET in df.columns:
        df.loc[df["status"] != "OPTIMAL", TARGET] = TIME_LIMIT


    group_cols = ["instance"]
    if "tal" in df.columns:
        group_cols.append("tal")
    if "settings" in df.columns:
        group_cols.append("settings")
    

    numeric_cols = [
        c for c in df.select_dtypes(include="number").columns
        if c not in group_cols
    ]
    other_cols = [
        c for c in df.columns
        if c not in group_cols and c not in numeric_cols and c != "status"
    ]

    agg_spec: dict = {c: "mean" for c in numeric_cols}
    agg_spec.update({c: "first" for c in other_cols})
    if "status" in df.columns:
        agg_spec["status"] = (
            lambda x: "OPTIMAL" if (x == "OPTIMAL").all() else "TIME_LIMIT"
        )

    return df.groupby(group_cols, as_index=False).agg(agg_spec)

# this is taken fromn across all but i changed the input parameter
def _write_per_stage_final_dumps(sources_full) -> None:
    """For every stage that has 3 FINAL-repetition directories, write one
    pair of intermediate CSVs **per setting present in the source dirs** into
    ``OUTPUT_ROOT / <stage>/``:

    * ``compiled_final_<setting>.csv``       — full-schema concatenation of
      all repetitions for that setting; an extra ``repetition`` column tags
      the source directory (trailing digit, e.g. ``1``/``2``/``3``).
    * ``compiled_final_<setting>_mean.csv``  — per-instance mean across
      repetitions (TL-penalised, full schema preserved via
      :func:`_mean_per_instance_full`).

    Settings are auto-detected from the source dirs (no PICKS filter) so e.g.
    stage 2 emits CSV pairs for ``benders_FLLP``, ``benders_LEQ`` *and*
    ``benders_LEQRF`` even though PICKS only lists two of them.

    Stages with no FINAL reps are skipped silently.
    """
    # Collect: stage -> set(unique source dirs)
    stage_sources: dict[str, list[Path]] = {}
    for src_full in sources_full:
        srcs: list[Path] = (
            src_full if isinstance(src_full, list)
            else ([] if src_full is None else [src_full])
        )
        if not srcs:
            continue
        stage = _stage_name_from_path(srcs[0])
        if stage is None:
            continue
        bucket = stage_sources.setdefault(stage, [])
        for s in srcs:
            if s not in bucket:
                bucket.append(s)

    for stage, srcs in stage_sources.items():
        # Read every repetition once
        all_rows: list[pd.DataFrame] = []
        for s in srcs:
            if not s.exists():
                continue
            try:
                part = read_experiment(s)
            except Exception as e:
                print(f"  [warn] {stage}/{s.name}: {e}")
                continue
            if part.empty:
                continue
            tail = s.name.rstrip("0123456789")
            tag = s.name[len(tail):] or s.name
            part = part.copy()
            part["repetition"] = tag
            part["source_dir"] = s.name
            all_rows.append(part)

        if not all_rows:
            print(f"  [skip] {stage}: no FINAL repetitions on disk")
            continue

        combined = pd.concat(all_rows, ignore_index=True)
        if "settings" not in combined.columns:
            print(f"  [warn] {stage}: no 'settings' column in source data")
            continue

        out_dir = OUTPUT_ROOT / stage
        out_dir.mkdir(parents=True, exist_ok=True)
        print(f"writing to {out_dir}...")
        for setting, sub in combined.groupby("settings", sort=True):
            safe = str(setting).replace("/", "_")
            setting_dir = out_dir / safe
            setting_dir.mkdir(parents=True, exist_ok=True)

            compiled_path = setting_dir / f"compiled_final_{safe}.csv"
            sub.to_csv(compiled_path, index=False)

            for_mean = sub.drop(columns=["repetition", "source_dir"],
                                errors="ignore")
            mean_df = _mean_per_instance_full(for_mean)
            mean_path = setting_dir / f"compiled_final_{safe}_mean.csv"
            mean_df.to_csv(mean_path, index=False)

            print(f"  [ok]   {stage}/{safe}: {len(sub)} raw / "
                  f"{len(mean_df)} mean rows")

_write_per_stage_final_dumps(SOURCES_FULL)