#!/usr/bin/env python3
"""Export the gap-rule comparison on the seed-42 parameter-search instances.

Requires Python 3.11. Reuses completed runs only after checking their identities,
settings, source, input hashes, objective bounds, and one-hour completion.
The grid's pending configurations and original submission records are unchanged.
"""

import argparse
import json
from pathlib import Path
import random
import tomllib

from export_random20_grid import (
    LABELS, LIMIT, METHODS, SHIFT, TOL, load_record, name_for,
    read_json, sgm, sha, write_csv, write_gap_panel,
)


GAP = "violation_pb_allfeas_110"


def write_panel(path, row, rho):
    write_gap_panel(path, dict(row, ratio=row['ratio_to_unfiltered']), rho)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--grid", required=True, type=Path)
    parser.add_argument("--gap", required=True, type=Path)
    parser.add_argument("--paper", required=True, type=Path)
    args = parser.parse_args()
    grid, gap = args.grid.resolve(), args.gap.resolve()
    source = read_json(grid / "source-provenance.json")
    gap_source = read_json(gap / "source-provenance.json")
    selection = read_json(grid / "instance-selection.json")
    cases = selection["selected"]
    assert selection["seed"] == 42 and len(cases) == len(set(cases)) == 20
    assert sorted(random.Random(42).sample(sorted(selection["population"]), 20)) == cases
    for name, key in (("instance-selection.json", "selection_sha256"),
                      ("reused-results.json", "reused_results_sha256")):
        assert sha((grid / name).read_bytes()) == source[key]

    settings_path = gap / "settings" / f"{GAP}.toml"
    settings_hash = sha(settings_path.read_bytes())
    assert settings_hash == gap_source["settings_sha256"][GAP]
    settings = tomllib.loads(settings_path.read_text())
    selector = settings["BENDERS"]["CUTS"]["SELECTION"]
    assert selector["ENABLED"] and selector["METHOD"] == ":VIOLATION"
    assert selector["FEASCUT_ONLY"] and selector["MIN_OPT_CUTS"] == 0
    assert selector["VIOLATION"]["LIMIT_MODE"] == ":pb"
    assert selector["VIOLATION"]["LIMIT_PB_RHO"] == 1.1
    assert settings["SOLVER"]["SOLVING_TIME_LIMIT"] == gap_source["time_limit"] == 10800
    assert sha((gap / "source.tar").read_bytes()) == gap_source["source_archive_sha256"]

    retained = read_json(grid / "reused-results.json")
    names = ["unfiltered"] + [name_for(m, "all", .05) for m in METHODS]
    assert {(e["method"], e["case"]) for e in retained} == {
        (name, case) for name in names for case in cases}
    assert len(retained) == 120
    records = [load_record(e, source, True) for e in retained]
    controls = {e["case"]: e for e in retained if e["method"] == "unfiltered"}
    expected_gap = [e for e in read_json(gap / "expected.json") if e["case"] in cases]
    assert len(expected_gap) == 20 and {e["case"] for e in expected_gap} == set(cases)
    for entry in expected_gap:
        assert entry["method"] == GAP
        assert entry["instance_sha256"] == controls[entry["case"]]["instance_sha256"]
        assert sha((gap / entry["instance"]).read_bytes()) == entry["instance_sha256"]
        assert gap_source["source_archive_sha256"] == controls[entry["case"]]["source_sha256"]
        output = Path(entry["output"])
        retained_gap = dict(
            entry, source_sha256=gap_source["source_archive_sha256"],
            settings_sha256=settings_hash, original_solver_limit=gap_source["time_limit"],
            result_sha256=sha((output / "result.json").read_bytes()),
            metadata_sha256=sha((output / "manifest.json").read_bytes()),
        )
        records.append(load_record(retained_gap, gap_source, True))
    for row in records:
        assert row["state"] == "complete" and row["solved"], row
    control_rows = {r["instance"]: r for r in records if r["configuration"] == "unfiltered"}
    for row in records:
        control = control_rows[row["instance"]]
        tolerance = TOL * max(1, abs(control["objective"]))
        assert abs(row["objective"] - control["objective"]) <= tolerance, row
        assert row["objective"] >= control["bound"] - tolerance, row
        assert row["bound"] <= control["objective"] + tolerance, row

    labels = ["No filtering", "SCIP cut selection", *LABELS[1:], "Gap rule"]
    summary = []
    for name, label in zip(names + [GAP], labels):
        rows = [r for r in records if r["configuration"] == name]
        assert len(rows) == 20
        summary.append(dict(configuration=name, label=label, solved=20,
                            shifted_mean=sgm([r["seconds"] for r in rows])))
    for row in summary:
        row["ratio_to_unfiltered"] = row["shifted_mean"] / summary[0]["shifted_mean"]

    output = grid / "analysis-gap-comparison"
    publication = args.paper.resolve() / "resources/grid_study_20261027"
    output.mkdir(parents=True, exist_ok=True)
    publication.mkdir(parents=True, exist_ok=True)
    write_csv(output / "records.csv", records)
    write_csv(output / "summary.csv", summary)
    audit = dict(
        grid_campaign=str(grid), gap_campaign=str(gap), instances=cases,
        shift=SHIFT, comparison_cutoff=LIMIT, original_time_limit=10800,
        objective_relative_tolerance=TOL, checked_runs=len(records),
        scope="Only the 20 selected instances; no assertion about the other gap-rule runs.",
        gap_rule="All feasibility cuts and no optimality cuts if any feasibility cuts exist. "
                 "Otherwise use decreasing physical optimality violation until its sum exceeds "
                 "1.1 times max(0, incumbent minus current master objective), or exhaust the pool. "
                 "Without finite bounds retain all candidates.",
        reuse_note="All retained runs finish within one hour despite their original three-hour limit.",
        gap_source_sha256=gap_source["source_archive_sha256"],
        gap_settings_sha256=settings_hash,
        generator_sha256=sha(Path(__file__).read_bytes()),
        input_sha256={str(p): sha(p.read_bytes()) for p in (
            grid / "source-provenance.json", grid / "instance-selection.json",
            grid / "reused-results.json", gap / "source-provenance.json", gap / "expected.json")},
    )
    (output / "audit.json").write_text(json.dumps(audit, indent=2) + "\n")
    write_panel(publication / "heatmap_pb.tikz", summary[-1],
                selector["VIOLATION"]["LIMIT_PB_RHO"])
    short_ids = [case.removeprefix("pglib_scopf_case") for case in cases]
    short_ids.sort(key=lambda value: tuple(int(n) for n in value.split("_")))
    escaped = [r"\texttt{" + value.replace("_", r"\_") + "}" for value in short_ids]
    (publication / "selected_instances.tex").write_text(
        "% Generated by check/scripts/export_gap_grid.py.\n"
        + ", ".join(escaped[:-1]) + ", and " + escaped[-1] + "%\n")
    print(json.dumps(summary, indent=2))
    print(f"Checked {len(records)} completed results. Diagnostics: {output}")


if __name__ == "__main__":
    main()
