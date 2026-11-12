#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# build_paper_resources.sh  —  uniform, zero-touch paper build.
#
# One mechanism, no per-stage special-casing:
#
#   1. refresh the CUT_ERROR exclusion lists (full + grid);
#   2. run every compare/heatmap script in PAPER_SCRIPTS.  Each reads the
#      read-only stage repos and writes tables + plots + heatmaps into the
#      1-baseline staging tree  check/results/paper/<stage>/{grid,full}/<exp>/
#      (raw benchmark data is only ever READ, never mutated);
#   3. mirror each produced paper leaf into  across_all/<stage>/<target>/  using
#      the PAPER_TABLES map (the paper \inputs across_all/*/comparison_short_all.tex);
#   4. run across_all.py for the cross-stage FINAL_GRID + compiled_final CSVs.
#
# To add a paper table: add its script to PAPER_SCRIPTS and, if it feeds the
# paper, one "paper_leaf | across_target" line to PAPER_TABLES.  Nothing else.
#
# Usage:  bash build_paper_resources.sh
# -----------------------------------------------------------------------------
set -euo pipefail
shopt -s nullglob extglob

VISU="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
RES="${TOPOLOGYOPTIMIZER_RESULTS_ROOT:-${VISU}/../../results}"
PAPER=$RES/paper
OUT=$RES/across_all

# -----------------------------------------------------------------------------
# Manifest
# -----------------------------------------------------------------------------
# Scripts that (re)generate paper resources, in run order.  One script may emit
# several leaves (e.g. compare_combined_grid, compare_grid); they are run once.
PAPER_SCRIPTS=(
    # headline
    FINAL/compare_full.py
    # stage 2
    2-SP_solving/compare_full.py
    # stage 3
    3-SP_selection/compare_1-across.py
    3-SP_selection/grid/compare_grid_score_individual.py
    3-SP_selection/grid/compare_grid_abs_rel.py
    3-SP_selection/grid/heatmap_grid_abs_rel.py
    3-SP_selection/grid/compare_grid_cold_diversity.py
    3-SP_selection/full/compare_full_sp_selection.py
    # stage 4  (combined grid: 5 leaves + per-leaf heatmaps)
    4-cut_filtering/compare_combined_grid.py
    # stage 5  (num + num_no_base + full + 2x2 heatmap)
    5-partial_decomp/compare_grid.py
    5-partial_decomp/compare_full.py
    5-partial_decomp/heatmaps/heatmap_5-partial_decomp_num_grid.py
    # stage 6
    6-seq_warmstarting/compare_1-tal1.py
    6-seq_warmstarting/compare_1-tal1-pd.py
    # stages 82-85
    82-fractional_sepa/compare_grid.py
    82-fractional_sepa/compare_full.py
    83-branching/compare_grid.py
    83-branching/compare_full.py
    84-in_out_strengthening/compare_grid.py
    85-benders-primal-heur/compare_grid.py
    # stage 9 (nested — extra, not in results_oli)
    9-nested/compare_grid.py
    9-nested/compare_full.py
)

# paper_leaf (under paper/)  |  across_all target (under across_all/)
# Exactly the leaves the paper \inputs get a mapping; everything else stays
# under paper/ as on-demand catalog.
PAPER_TABLES=(
    "FINAL/full|FINAL"
    "2-SP_solving/full|2-SP_solving/full"
    "3-SP_selection/grid/1-n_hits|3-SP_selection/1-n_hits"
    "3-SP_selection/grid/2-score-individual|3-SP_selection/2-score-individual"
    "4-cut_filtering/grid/6-cutfiltering-abs-spsel-rel|4-cut_filtering/6-cutfiltering-abs-spsel-rel"
    "5-partial_decomp/grid/num|5-partial_decomp/1-num"
    "5-partial_decomp/grid/num_no_base|5-partial_decomp/2-num_no-base"
    "6-seq_warmstarting/full/1-tal1|6-seq_warmstarting/6-seq_warmstart--1-tal1"
    "6-seq_warmstarting/full/1-tal1-pd|6-seq_warmstarting/6-seq_warmstart--1-tal1-pd"
    "82-fractional_sepa/grid/all|82-fractional_sepa/1-frac_sepa"
    "83-branching/grid/1-bprio_grid_div|83-branching/1-bprio_grid"
    "84-in_out_strengthening/grid/1-in_out|84-in_out_strengthening/1-in_out"
    "85-benders-primal-heur/grid/1-benders-primal-heur|85-benders-primal-heur/1-benders-primal-heur"
)
# NB: row 1 (baseline B&BC vs Compact) is deferred — a cross-repo merge of
# benders (2-SP_solving) + compact (1-baseline).  Its curated table stays in
# across_all/1-baseline/1-across and is not regenerated here.

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------
copy_artifacts() {
    local src="$1" dst="$2" f
    mkdir -p "$dst"
    for f in "$src"/comparison*.md "$src"/comparison*.tex \
             "$src"/comparison*.pdf "$src"/comparison*.png \
             "$src"/performance_profile_*.pdf "$src"/performance_profile_*.png \
             "$src"/*heatmap*.pdf "$src"/*heatmap*.png \
             "$src"/*.toml; do
        [[ -f "$f" ]] || continue
        cp -f "$f" "$dst/$(basename "$f")"
    done
    # log_metrics/ (stage 6) is a directory of derived plots — mirror wholesale.
    [[ -d "$src/log_metrics" ]] && { rm -rf "$dst/log_metrics"; cp -rf "$src/log_metrics" "$dst/log_metrics"; }
    return 0
}

# -----------------------------------------------------------------------------
# 1. CUT_ERROR exclusion lists (must precede every comparison)
# -----------------------------------------------------------------------------
echo "== CUT_ERROR exclusion lists"
python3 "$VISU/build_cut_err_exclusions.py"

# -----------------------------------------------------------------------------
# 2. Run the paper scripts (fresh paper/ so re-runs are reproducible)
# -----------------------------------------------------------------------------
rm -rf "$PAPER"
for s in "${PAPER_SCRIPTS[@]}"; do
    if [[ ! -f "$VISU/$s" ]]; then echo "  [warn] missing script: $s"; continue; fi
    echo "== run $s"
    ( cd "$VISU/$(dirname "$s")" && python3 "$(basename "$s")" )
done

# -----------------------------------------------------------------------------
# 3. Mirror the paper leaves into across_all/ (the paper's \input targets)
# -----------------------------------------------------------------------------
echo "== mirror paper/ -> across_all/"
for entry in "${PAPER_TABLES[@]}"; do
    leaf="${entry%%|*}"; target="${entry##*|}"
    if [[ ! -d "$PAPER/$leaf" ]]; then echo "  [warn] no paper leaf: $leaf"; continue; fi
    copy_artifacts "$PAPER/$leaf" "$OUT/$target"
    echo "  $leaf -> across_all/$target"
done

# -----------------------------------------------------------------------------
# 4. Cross-stage tables (FINAL_GRID) + compiled_final_<setting>{,_mean}.csv
# -----------------------------------------------------------------------------
echo "== across_all.py"
python3 "$VISU/across_all.py"

echo
echo "Done. Paper resources under $OUT/  (staging under $PAPER/)"
