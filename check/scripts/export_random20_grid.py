#!/usr/bin/env python3
"""Export the seed-42 budget grid, including its retained 5% and control runs.

Usage:
  python3 check/scripts/export_random20_grid.py --campaign PATH --paper PATH

Only a complete set of 20 runs with checked identities and settings produces a
cell. Completed failed runs receive the time limit, as do timeouts. Unfinished
runs and mismatched inputs remain unscored. The common budget pools Violation,
Diversity, and Hybrid only. All 18 budgets for these strategies must be complete
before selecting a winner. Raw records and checks stay in the experiment directory.
"""

import argparse
from collections import Counter
import csv
import hashlib
import json
import math
from pathlib import Path
import os
import random

from grid_significance import one_sided_wilcoxon_test, holm_adjust, ratio_label, significance_stars


METHODS = ("scip_hybrid", "random", "violation", "diversity", "hybrid")
TUNING_METHODS = ("violation", "diversity", "hybrid")
LABELS = ("SCIP", "Random", "Violation", "Diversity", "Hybrid")
VALUES = {"abs": (1, 2, 5, 10, 20, 50),
          "all": (.01, .05, .10, .25, .50, .75),
          "current": (.01, .05, .10, .25, .50, .75)}
PARAMETERS = {"abs": r"\cfmaxn", "all": r"\cfmaxrel", "current": r"\cfmaxreliter"}
SHIFT, LIMIT, TOL = 10.0, 3600.0, 1e-4
SCIP_REPLACEMENT = (Path(os.environ.get("TOPOLOGYOPTIMIZER_RESULTS_ROOT", Path(__file__).resolve().parents[2] / "check/results")) / '20261027-cutfilter-scip-policy-native-s42')
GAP_CAMPAIGN = (Path(os.environ.get("TOPOLOGYOPTIMIZER_RESULTS_ROOT", Path(__file__).resolve().parents[2] / "check/results")) / '20261023-cutfilter-pb-allfeas-b8a05381')
GAP_METHOD = 'violation_pb_allfeas_110'
SIGNIFICANCE_FAMILY = 92  # 90 budget settings, SCIP default limits, primal bound rule.


def sha(data):
    return hashlib.sha256(data).hexdigest()


def read_json(path):
    return json.loads(path.read_text())


def sgm(values):
    return math.exp(math.fsum(math.log(v + SHIFT) for v in values) / len(values)) - SHIFT


def finite(value):
    return isinstance(value, (int, float)) and math.isfinite(value)


def name_for(method, mode, value):
    number = value if mode == "abs" else round(100 * value)
    return f"{method}_{mode}_{number:03d}"


def common_budget_summary(mode, value, grouped, lookup, baseline_sgm):
    names = [name_for(method, mode, value) for method in TUNING_METHODS]
    complete = all(lookup[name]["complete"] for name in names)
    mean = sgm([r["scored_seconds"] for name in names for r in grouped[name]]) if complete else None
    return dict(mode=mode, value=value, complete=complete,
                shifted_mean=mean, ratio=mean / baseline_sgm if complete else None)


def penalize_failure(row, issue):
    row.update(state="failed", solved=False, scored_seconds=LIMIT, issue=issue)
    return row


def load_record(entry, source, reused):
    row = dict(configuration=entry["method"], instance=entry["case"],
               reused=reused, output=entry["output"], state="pending",
               status="", seconds=None, objective=None, bound=None,
               solved=False, scored_seconds=None, issue="",
               result_sha256="", metadata_sha256="")
    output = Path(entry["output"])
    try:
        metadata_raw = (output / "manifest.json").read_bytes()
        meta = json.loads(metadata_raw)
    except FileNotFoundError:
        return row
    except (ValueError, OSError) as error:
        row.update(state="invalid", issue=f"Unreadable metadata: {error}")
        return row
    if not meta.get("complete") and not (meta.get("end_epoch") and meta.get("returncode") is not None):
        return row
    row["metadata_sha256"] = sha(metadata_raw)
    try:
        expected_source = entry["source_sha256"] if reused else source["source_archive_sha256"]
        expected_settings = (entry["settings_sha256"] if reused
                             else source["settings_sha256"][entry["method"]])
        expected_limit = entry["original_solver_limit"] if reused else entry.get("solver_limit", LIMIT)
        checks = {
            "identity": (meta["method"], meta["case"]) == (entry["method"], entry["case"]),
            "source": meta["source_sha256"] == expected_source,
            "settings": meta["settings_sha256"] == expected_settings,
            "input": meta["input_sha256"] == entry["instance_sha256"],
            "time limit": meta["time_limit"] == expected_limit,
            "seed": meta["seed"] == 42,
            "threads": meta["master_threads"] == 4 and meta["julia_threads"] == 1,
            "switching limit": meta["tal"] == 1,
        }
        if reused:
            checks["metadata"] = row["metadata_sha256"] == entry["metadata_sha256"]
        failed = [key for key, ok in checks.items() if not ok]
        if failed:
            raise ValueError("Mismatched " + ", ".join(failed))
        returncode = meta["returncode"]
    except (KeyError, ValueError) as error:
        row.update(state="invalid", issue=str(error))
        return row
    try:
        raw = (output / "result.json").read_bytes()
        row["result_sha256"] = sha(raw)
        if reused and row["result_sha256"] != entry["result_sha256"]:
            row.update(state="invalid", issue="Mismatched result")
            return row
        result = json.loads(raw)
        if not isinstance(result, dict):
            raise ValueError("Result is not an object")
    except (ValueError, OSError) as error:
        return penalize_failure(row, f"Completed run without a readable result: {error}")
    row.update(status=result.get("termination_status", ""), seconds=result.get("solve_time"),
               objective=result.get("objective"), bound=result.get("objective_lb"))
    if returncode != 0:
        return penalize_failure(row, f"Nonzero return code: {returncode}")
    try:
        seconds = result.get("solve_time")
        if not finite(seconds) or seconds < 0:
            raise ValueError("Missing or invalid solve time")
        status = result["termination_status"]
        if status == "OPTIMAL":
            if not finite(row["objective"]) or not finite(row["bound"]):
                raise ValueError("Missing objective or bound")
            tolerance = TOL * max(1.0, abs(row["objective"]))
            if abs(row["objective"] - row["bound"]) > tolerance:
                raise ValueError("Optimal status with inconsistent objective and bound")
            if result.get("primal_status") != "FEASIBLE_POINT":
                raise ValueError("Optimal status without a feasible primal point")
        elif status != "TIME_LIMIT":
            raise ValueError(f"Unsuccessful termination: {status}")
        if reused and (status != "OPTIMAL" or seconds > LIMIT):
            raise ValueError("Retained result is not solved within the comparison cutoff")
        row.update(state="complete", solved=status == "OPTIMAL" and seconds <= LIMIT)
        row["scored_seconds"] = seconds if row["solved"] else LIMIT
    except (KeyError, ValueError) as error:
        penalize_failure(row, str(error))
    return row


def check_against_control(row, control):
    if row["state"] != "complete":
        return
    tolerance = TOL * max(1.0, abs(control["objective"]))
    objective, bound = row["objective"], row["bound"]
    incompatible = (
        (row["status"] == "OPTIMAL" and abs(objective - control["objective"]) > tolerance)
        or (finite(objective) and objective < control["bound"] - tolerance)
        or (finite(bound) and bound > control["objective"] + tolerance))
    if incompatible:
        penalize_failure(row, "Objective or bound disagrees with the retained control")


def load_gap_records(campaign, cases, source, controls):
    provenance = read_json(campaign / 'source-provenance.json')
    assert sha((campaign / 'source.tar').read_bytes()) == provenance['source_archive_sha256']
    assert sha((campaign / 'settings' / (GAP_METHOD + '.toml')).read_bytes()) == provenance['settings_sha256'][GAP_METHOD]
    entries = [e for e in read_json(campaign / 'expected.json') if e['case'] in cases]
    assert len(entries) == 20 and {e['case'] for e in entries} == set(cases)
    assert {e['method'] for e in entries} == {GAP_METHOD}
    rows = []
    for entry in entries:
        assert entry['instance_sha256'] == source['input_sha256'][entry['case'] + '.m']
        assert sha((campaign / entry['instance']).read_bytes()) == entry['instance_sha256']
        row = load_record(dict(entry, solver_limit=10800), provenance, False)
        check_against_control(row, controls[row['instance']])
        rows.append(row)
    return rows


def annotate_significance(cells, grouped, controls):
    cases = sorted(controls)
    control_times = [controls[case]['scored_seconds'] for case in cases]
    pvalues = {}
    for cell in cells:
        cell.update(p_value=None, p_adjusted=None, significance_stars='',
                    wilcoxon_statistic=None, wilcoxon_method=None, nonzero_pairs=None, absolute_ties=None)
        name = cell['configuration']
        if name == 'unfiltered' or not cell['complete']:
            continue
        rows = grouped[name]
        assert len(rows) == len(cases) and {r['instance'] for r in rows} == set(cases)
        paired = {row['instance']: row['scored_seconds'] for row in rows}
        test = one_sided_wilcoxon_test([paired[case] for case in cases], control_times)
        cell.update(test)
        pvalues[name] = test['p_value']
    adjusted = holm_adjust(pvalues, SIGNIFICANCE_FAMILY)
    for cell in cells:
        name = cell['configuration']
        if name in pvalues:
            cell.update(p_value=pvalues[name], p_adjusted=adjusted[name],
                        significance_stars=significance_stars(adjusted[name]))


def color(ratio):
    signed = max(-1.0, min(1.0, math.log2(max(ratio, 1e-12))))
    endpoint = (33, 102, 172) if signed < 0 else (179, 88, 6)
    rgb = tuple(round(255 + abs(signed) * (component - 255)) for component in endpoint)
    foreground = "white" if sum(a * b for a, b in zip((.299, .587, .114), rgb)) < 128 else "black"
    return "".join(f"{component:02X}" for component in rgb), foreground


def write_panel(path, mode, lookup, winners, baseline):
    lines = ["% Generated by check/scripts/export_random20_grid.py.",
             r"\begin{tikzpicture}[x=1.2cm,y=.52cm,font=\small]",
             r"\path[use as bounding box] (-1.55,-.04) rectangle (7.4,5.65);"]
    lines.append(rf"\node[anchor=east] at (-.12,5.35) {{${PARAMETERS[mode]}$}};")
    for x, value in enumerate(VALUES[mode]):
        label = rf"${value:g}$"
        lines.append(rf"\node[font=\scriptsize] at ({x + .5},5.35) {{{label}}};")
    for i, (method, label) in enumerate(zip(METHODS, LABELS)):
        y = 4 - i
        lines.append(rf"\node[anchor=east] at (-.12,{y + .5}) {{{label}}};")
        for x, value in enumerate(VALUES[mode]):
            name = name_for(method, mode, value)
            cell = lookup[name]
            ratio = cell["ratio"]
            if ratio is None:
                lines.extend([
                    rf"\filldraw[fill=black!7,draw=white] ({x},{y}) rectangle ({x + 1},{y + 1});",
                    rf"\node[text=black!45] at ({x + .5},{y + .5}) {{---}};"])
                continue
            fill, foreground = color(ratio)
            lines.extend([
                rf"\definecolor{{gridtile}}{{HTML}}{{{fill}}}",
                rf"\filldraw[fill=gridtile,draw=white] ({x},{y}) rectangle ({x + 1},{y + 1});",
                rf"\node[text={foreground},font=\footnotesize] at ({x + .5},{y + .67}) {{{ratio_label(cell)}}};",
                rf"\node[text={foreground},font=\tiny] at ({x + .5},{y + .24}) {{$S={cell['solved']}$}};"])
            if name in winners.get(method, []):
                lines.append(rf"\draw[black,line width=.8pt] ({x + .025},{y + .025}) rectangle ({x + .975},{y + .975});")
    lines.extend([
        r"\filldraw[fill=white,draw=black!35] (6.2,0) rectangle (7.4,5);",
        r"\node[align=center,font=\footnotesize] at (6.8,3.4) {No\\filtering};",
        r"\node[font=\footnotesize] at (6.8,2.25) {1.00};",
        rf"\node[font=\scriptsize] at (6.8,1.25) {{$S={baseline['solved']}$}};",
        r"\end{tikzpicture}"])
    path.write_text("\n".join(lines) + "\n")


def write_legend(path, pending, invalid, winners):
    lines = ["% Generated by check/scripts/export_random20_grid.py.",
             r"\begin{tikzpicture}[x=1cm,y=1cm,font=\scriptsize]",
             r"\node[anchor=east] at (-2.15,-.15) {Relative SGM};"]
    for i, ratio in enumerate((.5, .75, 1, 1.5, 2)):
        fill, _ = color(ratio)
        x = i * .8 - 2
        label = r"$\leq 0.5$" if i == 0 else (r"$\geq 2$" if i == 4 else f"{ratio:g}")
        lines.extend([
            rf"\definecolor{{gridtile}}{{HTML}}{{{fill}}}",
            rf"\filldraw[fill=gridtile,draw=black!25] ({x},-.25) rectangle ({x + .8},-.05);",
            rf"\node at ({x + .4},-.45) {{{label}}};"])
    lines.append(r"\end{tikzpicture}")
    path.write_text("\n".join(lines) + "\n")


def write_csv(path, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def replace_scip_records(records, configs, cases, campaign):
    source = read_json(campaign / 'source-provenance.json')
    assert sha((campaign / 'source.tar').read_bytes()) == source['source_archive_sha256']
    assert read_json(campaign / 'instance-selection.json')['selected'] == cases
    entries = [e for e in read_json(campaign / 'expected.json') if e['phase'] == 'grid']
    names = {name_for('scip_hybrid', mode, value) for mode, values in VALUES.items() for value in values}
    names.add('scip_hybrid_native_grid')
    assert len(entries) == 380
    assert {(e['method'],e['case']) for e in entries} == {(name,case) for name in names for case in cases}
    for name in names:
        assert sha((campaign / 'settings' / (name+'.toml')).read_bytes()) == source['settings_sha256'][name]
    # Replace all old SCIP results, including previously reused 5% runs.
    records = [r for r in records if not r['configuration'].startswith('scip_hybrid_')]
    native_source = campaign / 'native-grid-source.json'
    if native_source.exists():
        mapping = read_json(native_source)
        assert mapping['source_archive_sha256'] == source['source_archive_sha256']
        assert mapping['selection_sha256'] == sha((campaign / 'instance-selection.json').read_bytes())
        assert mapping['sampled_instances'] == sorted(cases)
        assert mapping['solver_limit'] == 10800 and mapping['comparison_cutoff'] == LIMIT
        full = [e for e in read_json(campaign / 'expected.json')
                if e['method'] == mapping['method'] and e['case'] in cases]
        assert len(full) == 20 and {e['case'] for e in full} == set(cases)
        assert all(e['solver_limit'] == 10800 for e in full)
        assert sha((campaign / 'settings' / (mapping['method']+'.toml')).read_bytes()) == source['settings_sha256'][mapping['method']]
        entries = [e for e in entries if e['method'] != 'scip_hybrid_native_grid'] + full
    for entry in entries:
        row = load_record(entry, source, False)
        if native_source.exists() and entry['method'] == mapping['method']:
            row['configuration'] = 'scip_hybrid_native_grid'
        records.append(row)
    configs = configs + [dict(name='scip_hybrid_native_grid',method='scip_native',budget_mode='native',budget_value=None)]
    return records, configs


def write_single_panel(path, cell, row_label, column_label, parameter=None):
    lines = ['% Generated by check/scripts/export_random20_grid.py.',
        r'\begin{tikzpicture}[x=1.2cm,y=.52cm,font=\small]',
        r'\path[use as bounding box] (-1.55,-.04) rectangle (1,1.65);']
    if parameter:
        lines.append(rf'\node[anchor=east] at (-.12,1.35) {{${parameter}$}};')
    lines += [rf'\node[font=\scriptsize] at (.5,1.35) {{{column_label}}};',
              rf'\node[anchor=east] at (-.12,.5) {{{row_label}}};']
    if cell['ratio'] is None:
        lines += [r'\filldraw[fill=black!7,draw=white] (0,0) rectangle (1,1);',
                  r'\node[text=black!45] at (.5,.5) {---};']
    else:
        fill, foreground = color(cell['ratio'])
        lines += [rf'\definecolor{{gridtile}}{{HTML}}{{{fill}}}',
            r'\filldraw[fill=gridtile,draw=white] (0,0) rectangle (1,1);',
            rf"\node[text={foreground},font=\footnotesize] at (.5,.67) {{{ratio_label(cell)}}};",
            rf"\node[text={foreground},font=\tiny] at (.5,.24) {{$S={cell['solved']}$}};"]
    lines.append(r'\end{tikzpicture}')
    path.write_text('\n'.join(lines)+'\n')


def write_native_panel(path, cell):
    write_single_panel(path, cell, 'SCIP', 'Default')


def write_gap_panel(path, cell, rho=1.1):
    write_single_panel(path, cell, 'Violation', f'${rho:g}$', r'\cfpb')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campaign", type=Path, required=True)
    parser.add_argument("--paper", type=Path, required=True)
    parser.add_argument("--output", type=Path, help="Diagnostic output directory")
    parser.add_argument('--significance', action='store_true',
                        help='Opt in to significance tests and stars (excluded by default)')
    parser.add_argument("--scip-campaign", type=Path, default=SCIP_REPLACEMENT,
                        help="Replacement SCIP grid and fixed default setting")
    args = parser.parse_args()
    campaign = args.campaign.resolve()
    output = args.output or campaign / "analysis-paper-grid"
    publication = args.paper.resolve() / "resources/grid_study_20261027"
    source = read_json(campaign / "source-provenance.json")
    selection = read_json(campaign / "instance-selection.json")
    cases = selection["selected"]
    assert selection["seed"] == 42 and len(cases) == 20
    assert sorted(random.Random(42).sample(sorted(selection["population"]), 20)) == cases
    configs = read_json(campaign / "configurations.json")
    expected_names = {"unfiltered"} | {
        name_for(method, mode, value)
        for method in METHODS for mode, values in VALUES.items() for value in values}
    assert {c["name"] for c in configs} == expected_names and len(configs) == 91
    new = read_json(campaign / "expected.json")
    retained = read_json(campaign / "reused-results.json")
    assert sha((campaign / "reused-results.json").read_bytes()) == source["reused_results_sha256"]
    assert sha((campaign / "instance-selection.json").read_bytes()) == source["selection_sha256"]
    entries = new + retained
    identities = {(e["method"], e["case"]) for e in entries}
    assert len(entries) == len(identities) == 1820
    assert identities == {(name, case) for name in expected_names for case in cases}
    records = [load_record(e, source, reused)
               for group, reused in ((new, False), (retained, True)) for e in group]
    replacement = args.scip_campaign
    if replacement.exists():
        records, configs = replace_scip_records(records, configs, cases, replacement)
        expected_names.add('scip_hybrid_native_grid')
    controls = {r["instance"]: r for r in records if r["configuration"] == "unfiltered"}
    assert all(r["state"] == "complete" and r["solved"] for r in controls.values())
    for row in records:
        check_against_control(row, controls[row["instance"]])
    grouped = {name: [r for r in records if r["configuration"] == name] for name in expected_names}
    baseline_sgm = sgm([r["scored_seconds"] for r in controls.values()])
    summary = []
    for config in configs:
        rows = grouped[config["name"]]
        complete = all(r["state"] in ("complete", "failed") for r in rows)
        mean = sgm([r["scored_seconds"] for r in rows]) if complete else None
        summary.append(dict(configuration=config["name"], method=config["method"],
                            mode=config["budget_mode"], value=config["budget_value"],
                            complete=complete,
                            available=sum(r["state"] in ("complete", "failed") for r in rows),
                            failed=sum(r["state"] == "failed" for r in rows),
                            solved=sum(r["solved"] for r in rows),
                            shifted_mean=mean, ratio=mean / baseline_sgm if complete else None))
    lookup = {r["configuration"]: r for r in summary}
    gap_rows = load_gap_records(GAP_CAMPAIGN, cases, source, controls)
    gap_complete = all(row['state'] in ('complete', 'failed') for row in gap_rows)
    gap_mean = sgm([row['scored_seconds'] for row in gap_rows]) if gap_complete else None
    gap_cell = dict(configuration=GAP_METHOD, complete=gap_complete,
                    ratio=gap_mean / baseline_sgm if gap_complete else None,
                    solved=sum(row['solved'] for row in gap_rows))
    if args.significance:
        annotate_significance(summary + [gap_cell], dict(grouped, **{GAP_METHOD: gap_rows}), controls)
    budgets = []
    winners = {}
    for method in METHODS:
        options = [r for r in summary if r["method"] == method]
        if all(r["complete"] for r in options):
            best = min(r["shifted_mean"] for r in options)
            winners[method] = [r["configuration"] for r in options
                               if math.isclose(r["shifted_mean"], best, rel_tol=1e-12)]
    for mode, values in VALUES.items():
        for value in values:
            budgets.append(common_budget_summary(mode, value, grouped, lookup, baseline_sgm))
    best_budget = None
    if all(b["complete"] for b in budgets):
        best_budget = min(budgets, key=lambda b: b["shifted_mean"])
    counts = dict(Counter(r["state"] for r in records))
    audit = dict(campaign=str(campaign), instances=len(cases), configurations=len(configs),
                 sampling_seed=42, shift=SHIFT, comparison_cutoff=LIMIT,
                 objective_relative_tolerance=TOL, states=counts,
                 complete_configurations=sum(r["complete"] for r in summary),
                 retained_results=len(retained), baseline_shifted_mean=baseline_sgm,
                 best_common_budget=best_budget, strategy_winners=winners,
                 common_budget_methods=list(TUNING_METHODS),
                 common_budget_metric="Shifted geometric mean of the 60 solve times "
                                      "(20 instances times Violation, Diversity, and Hybrid), using 3600 for timeouts "
                                      "and completed failed runs.",
                 failure_policy="Completed failed runs receive 3600 seconds and do not count as solved, "
                                "including solver errors and objective or bound inconsistencies.",
                 missing_policy="Unfinished runs and runs with mismatched or unreadable metadata "
                                "remain unscored, and their configurations are omitted from plot cells. "
                                "No ranking until all competing budgets are complete.",
                 reuse_note=source["reused_time_limit_note"],
                 generator_sha256=sha(Path(__file__).read_bytes()),
                 input_sha256={name: sha((campaign / name).read_bytes()) for name in (
                     "source-provenance.json", "instance-selection.json", "configurations.json",
                     "expected.json", "reused-results.json")})
    audit['significance'] = dict(enabled=False)
    if args.significance:
        import scipy
        audit['significance'] = dict(enabled=True,
            test='One-sided Wilcoxon signed-rank test on paired scored solve times',
            comparator='No filtering on the same 20 instances', adjustment='Holm',
            planned_comparisons=SIGNIFICANCE_FAMILY,
            tested_comparisons=sum(cell['p_value'] is not None for cell in summary + [gap_cell]),
            thresholds=[.05, .01, .001],
            alternative='Configuration minus no-filtering solve time is less than zero',
            direction='Improvement over no filtering', time_transformation='None',
            zero_method='wilcox', ties='Normal approximation when zeros or tied absolute differences occur',
            missing_tests='Untested planned comparisons are retained in the Holm family size',
            interpretation='Exploratory comparisons over instances, with one run per configuration and instance',
            scipy_version=scipy.__version__,
            helper_sha256=sha(Path(__file__).with_name('grid_significance.py').read_bytes()),
            gap_campaign=str(GAP_CAMPAIGN), gap_cell=gap_cell,
            gap_input_sha256={name: sha((GAP_CAMPAIGN / name).read_bytes()) for name in (
                'source-provenance.json', 'expected.json', 'settings/' + GAP_METHOD + '.toml')})
    if replacement.exists():
        audit['scip_replacement'] = str(replacement)
        audit['scip_replacement_sha256'] = sha((replacement / 'source-provenance.json').read_bytes())
        audit['native_setting_is_tuned'] = False
        if (replacement / 'native-grid-source.json').exists():
            audit['native_grid_source'] = read_json(replacement / 'native-grid-source.json')
    output.mkdir(parents=True, exist_ok=True)
    publication.mkdir(parents=True, exist_ok=True)
    write_csv(output / "records.csv", records)
    write_csv(output / "summary.csv", summary)
    write_csv(output / "common-budgets.csv", budgets)
    write_csv(output / 'gap-records.csv', gap_rows)
    if args.significance:
        write_csv(output / 'significance.csv', [dict(configuration=cell['configuration'],
            complete=cell['complete'], p_value=cell['p_value'], p_adjusted=cell['p_adjusted'],
            stars=cell['significance_stars'], statistic=cell['wilcoxon_statistic'],
            method=cell['wilcoxon_method'], nonzero_pairs=cell['nonzero_pairs'],
            absolute_ties=cell['absolute_ties']) for cell in summary + [gap_cell]
            if cell['configuration'] != 'unfiltered'])
    (output / "audit.json").write_text(json.dumps(audit, indent=2) + "\n")
    for mode in VALUES:
        write_panel(publication / f"heatmap_{mode}.tikz", mode, lookup, winners, lookup["unfiltered"])
    write_legend(publication / "heatmap_legend.tikz", counts.get("pending", 0),
                 counts.get("invalid", 0), winners)
    if 'scip_hybrid_native_grid' in lookup:
        write_native_panel(publication / 'heatmap_scip_native.tikz', lookup['scip_hybrid_native_grid'])
    write_gap_panel(publication / 'heatmap_pb.tikz', gap_cell)
    print(json.dumps({key: audit[key] for key in (
        "states", "complete_configurations", "baseline_shifted_mean", "best_common_budget")}, indent=2))
    print(f"Plots: {publication}\nDiagnostics: {output}")


if __name__ == "__main__":
    main()
