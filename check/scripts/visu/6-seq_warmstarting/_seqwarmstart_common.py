"""Shared label/style/discovery helpers for the 6-seq_warmstarting experiments.

The ``1-tal1`` sweep varies the warmstart selection pool size (``n``) and the
weight scheme (``w``) used by `:rscots` warmstarting. Setting names follow:

    rscots_n{N}_w_{W}     N in {4, 8, 16, 32}, W in {full, rnd, sens, congprox, eta}

plus the reference settings ``no_warmstart``, ``rscots`` (defaults), and ``scopf`` (SCOPF warmstart).

Labels are rendered as ``n=N / w=W``. Colours form a darker gradient with
larger ``n``; line styles cycle with ``w``.
"""

from __future__ import annotations

import re
from pathlib import Path


_PATTERN = re.compile(r"^rscots_n(?P<n>\d+)_w_(?P<w>[a-z]+)$")

# Display order for the weight scheme.
W_ORDER = ["full", "rnd", "sens", "congprox", "eta"]


def parse_setting(name: str) -> dict | None:
    """Return ``{"n": int, "w": str}`` for a grid setting, else ``None``."""
    m = _PATTERN.match(name)
    if not m:
        return None
    return {"n": int(m.group("n")), "w": m.group("w")}


def format_label(name: str) -> str:
    info = parse_setting(name)
    if info is None:
        return rf"\texttt{{{name}}}"
    return rf"\texttt{{n={info['n']} / w={info['w']}}}"


# ---- Special reference settings ------------------------------------------
# Settings not belonging to the rscots_nN_w_W grid get their own style.
_SPECIAL_STYLES: dict[str, dict] = {
    "no_warmstart": {"color": "#000000", "linestyle": "-"},
    "rscots":       {"color": "#444444", "linestyle": "-"},
    "rscots_full":  {"color": "#444444", "linestyle": "-."},
    "scopf":        {"color": "#009E73", "linestyle": "-"},
}
_SPECIAL_LABELS: dict[str, str] = {
    "no_warmstart": r"\texttt{No warmstart}",
    "rscots":       r"\texttt{RSCOTS}",
    "rscots_full":  r"\texttt{RSCOTS (full)}",
    "scopf":        r"\texttt{SCOPF}",
}

# ---- Palettes -------------------------------------------------------------
# One colour ramp keyed by n (darker = larger pool).
_HUE_BY_N = {
    4:  "#bdbdbd",
    8:  "#6baed6",
    16: "#fb6a4a",
    32: "#252525",
}
_LINESTYLES = ["-", "--", "-.", ":", (0, (3, 1, 1, 1))]


def build_style_map(settings: list[str]) -> dict[str, dict]:
    out: dict[str, dict] = {}
    for s in settings:
        if s in _SPECIAL_STYLES:
            out[s] = _SPECIAL_STYLES[s]
            continue
        info = parse_setting(s)
        if info is None:
            out[s] = {"color": "#444444", "linestyle": "-"}
            continue
        color = _HUE_BY_N.get(info["n"], "#444444")
        idx = W_ORDER.index(info["w"]) if info["w"] in W_ORDER else 0
        out[s] = {"color": color, "linestyle": _LINESTYLES[idx % len(_LINESTYLES)]}
    return out


def build_label_map(settings: list[str]) -> dict[str, str]:
    return {s: _SPECIAL_LABELS.get(s, format_label(s)) for s in settings}


def discover_settings(results_dir: Path, fallback_dir: Path | None = None
                      ) -> list[str]:
    """List setting names from a results dir; fall back to settings dir."""
    _SKIP = {"log_metrics", "old_stats"}
    if results_dir.is_dir():
        names = sorted(d.name for d in results_dir.iterdir()
                       if d.is_dir() and not d.name.startswith(".")
                       and d.name not in _SKIP
                       and not d.name.endswith("_infeasible"))
        if names:
            return names
    if fallback_dir and fallback_dir.is_dir():
        return sorted(p.stem for p in fallback_dir.glob("*.toml"))
    return []


# ---- Shared experiment runner --------------------------------------------
# ``compare_1-tal1.py`` and ``compare_1-tal1-pd.py`` previously duplicated
# ~80 lines of boilerplate; both now defer to :func:`run_experiment`.

def ordered_grid_settings(discovered: list[str]) -> list[str]:
    """Canonical display order: baselines → n × w grid → ``rscots_full``.

    Filters down to settings actually present in ``discovered``.
    """
    disc = set(discovered)
    order = (["no_warmstart", "scopf"]
             + [f"rscots_n{n}_w_{w}" for n in [4, 8, 16, 32] for w in W_ORDER]
             + ["rscots_full"])
    return [s for s in order if s in disc]


def _add_time_with_warmstart(df):
    """Derive ``time_with_warmstart = time_solve + time_warmstart``.

    All settings except ``no_warmstart`` include warmstart overhead.
    """
    import pandas as pd
    df = df.copy()
    t_solve = pd.to_numeric(df.get("time_solve", 0), errors="coerce").fillna(0)
    t_ws = pd.to_numeric(df.get("time_warmstart", 0), errors="coerce").fillna(0)
    # include warmstart time for EVERY run (time_warmstart is 0 for
    # non-warmstarted settings, so this is a no-op there)
    df["time_with_warmstart"] = t_solve + t_ws
    return df


def run_experiment(
    experiment_dir: Path,
    settings_dir: Path,
    baseline: str = "no_warmstart",
    pd_convergence: bool = False,
    pd_top_n: int | None = None,
    pd_include_warmstart: bool = False,
    pd_per_instance: bool = False,
    min_time: float = 0,
    min_buses: int = 73,
    use_csvs_if_possible: bool = True,
    output_dir: Path | None = None,
) -> None:
    """Run log_metrics + full comparison for a 6-seq_warmstart sweep.

    Discovers settings under ``experiment_dir`` (falling back to
    ``settings_dir`` for the listing), orders them via
    :func:`ordered_grid_settings`, and dispatches to ``run_log_metrics`` and
    ``run_full_comparison`` (sharing the parsed-log cache between the two so
    Gurobi logs are parsed once).

    ``pd_include_warmstart=False`` (the default for this experiment family)
    drops the warmstart sub-MIP prefix from the PD convergence plot — the
    main solve starts at ``t=0`` for all settings — so the curves visualize
    the *response* of the main MP to the warmstart it received.  Switch to
    ``True`` to keep the wall-time-faithful timeline (warmstart period
    visible).

    ``use_csvs_if_possible=True`` (default) loads each ``.out`` log from its
    pre-parsed sibling ``logs_parsed/<stem>.csv`` when one exists, so
    re-runs of the visualization pipeline skip the (slow on NFS) log-parsing
    pass.  When a CSV is missing, the raw log is parsed and the CSV is
    written for next time.  Pass ``False`` to force a full re-parse — this
    also overwrites every existing CSV.
    """
    from log_metrics import run as run_log_metrics
    from tables import run_full_comparison

    include = ordered_grid_settings(discover_settings(experiment_dir, settings_dir))
    if not include:
        print(f"warning: no settings discovered for {experiment_dir.name}")
        return
    if not experiment_dir.is_dir() or not any(experiment_dir.iterdir()):
        print(f"warning: no results in {experiment_dir}; skipping comparison")
        return

    style_map = build_style_map(include)
    label_map = build_label_map(include)

    try:
        _, raw_cache, logs_map = run_log_metrics(
            experiment_dir, baseline=baseline,
            label_map=label_map, include_settings=include,
            use_csvs_if_possible=use_csvs_if_possible,
            output_dir=(Path(output_dir) / "log_metrics") if output_dir else None)
        from log_metrics import _scopf_threshold_per_instance
        scopf_thresholds = _scopf_threshold_per_instance(logs_map, raw_cache)
        run_full_comparison(
            experiment_dir, output_dir=output_dir,
            baseline=baseline, target="time_with_warmstart",
            pd_convergence=pd_convergence, pd_top_n=pd_top_n,
            pd_include_warmstart=pd_include_warmstart,
            pd_per_instance=pd_per_instance,
            pd_scopf_thresholds=scopf_thresholds,
            style_map=style_map, label_map=label_map,
            include_settings=include,
            min_buses=min_buses, min_time=min_time,
            log_raw_cache=raw_cache, log_files_map=logs_map,
            df_transform=_add_time_with_warmstart,
        )
    except FileNotFoundError as e:
        print(f"warning: {e}")
