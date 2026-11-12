"""Shared label/style derivation for the six 4-cut_filtering experiments.

Each experiment runs a grid over the cut-filtering (``cf``) and partial
subproblem-selection (``sp``) limits. Setting names follow:

    cf{NN|NNN}_sp{NN|NNN}    e.g. cf001_sp025, cf05_sp10
    sp{NN|NNN}               (for the two pb-based experiments where the
                              cut filter is fixed to the PB-aware variant)

3-digit numeric tokens encode the *relative* limit as a zero-padded percent
(001 → 1%, 100 → 100%), 2-digit tokens encode the *absolute* limit. Labels
are formatted accordingly. Colours form a darker gradient with larger CF,
line styles cycle with SP. The complete sweep is discovered from the
experiment results directory (or, if empty, from the settings directory).
"""

from __future__ import annotations

import re
from pathlib import Path


_TOKEN_RE = re.compile(r"^(?P<kind>cf|sp)(?P<num>\d+)$")


def parse_token(tok: str) -> tuple[str, str, str]:
    """Return (kind, mode, display) for a single 'cfNNN' / 'spNN' token.

    - mode == "rel" for 3-digit tokens (zero-padded percent, 001..100)
    - mode == "abs" for 2-digit tokens (absolute count)
    """
    m = _TOKEN_RE.match(tok)
    if not m:
        return ("?", "?", tok)
    kind, raw = m.group("kind"), m.group("num")
    if len(raw) == 3:  # relative, percent (zero-padded): 001 -> 1%, 100 -> 100%
        pct = int(raw)
        disp = f"{pct:g}\\%"
        return (kind, "rel", disp)
    # absolute
    return (kind, "abs", raw)


def parse_setting(name: str) -> dict:
    """Decompose a setting name into ``{cf?, sp?, rel?}`` info dicts.

    A bare trailing numeric token (e.g. ``005`` in ``cf001_sp05_005``) is
    captured as ``rel`` and rendered as a percentage in labels. Such tokens
    encode the ``STOP_AFTER_N_HITS_REL`` value of the partial subproblem
    solving stop rule.
    """
    info = {}
    for tok in name.split("_"):
        kind, mode, disp = parse_token(tok)
        if kind in ("cf", "sp"):
            info[kind] = {"raw": tok[2:], "mode": mode, "display": disp}
        elif kind == "?" and tok.isdigit():
            info["rel"] = {"raw": tok, "display": f"{int(tok):g}\\%"}
    return info


def format_label(name: str) -> str:
    """LaTeX label for a setting name like 'cf005_sp025' or 'cf001_sp05_005'."""
    info = parse_setting(name)
    parts = []
    if "cf" in info:
        parts.append(f"CF={info['cf']['display']}")
    if "sp" in info:
        parts.append(f"SP={info['sp']['display']}")
    if "rel" in info:
        parts.append(f"rel={info['rel']['display']}")
    body = " / ".join(parts) if parts else name
    return rf"\texttt{{{body}}}"


# ---- Palettes -------------------------------------------------------------
# One colour ramp per CF level (dark = high CF); fallback grey for SP-only.
_CF_RAMPS = {
    # cf raw token -> hex; covers the common rel & abs values seen in sweeps
    "001": "#bdbdbd", "005": "#969696", "010": "#737373",
    "025": "#525252", "050": "#252525", "075": "#000000", "100": "#000000",
    "01":  "#bdbdbd", "02":  "#969696", "05":  "#525252", "10":  "#000000",
}
_HUE_BY_CF = {
    "001": "#bdbdbd", "005": "#a1d99b", "010": "#6baed6",
    "025": "#fb6a4a", "050": "#c51b8a", "075": "#E69F00", "100": "#000000",
    "01":  "#a1d99b", "02":  "#6baed6", "05":  "#fb6a4a", "10":  "#c51b8a",
}
_LINESTYLES = ["-", "--", "-.", ":", (0, (3, 1, 1, 1)), (0, (1, 1))]


def _sp_index(sp_raw: str, all_sp: list[str]) -> int:
    try:
        return all_sp.index(sp_raw)
    except ValueError:
        return 0


def build_style_map(settings: list[str]) -> dict[str, dict]:
    """Per-setting style: colour ramps with CF, line style cycles with SP/rel."""
    all_keys: list[tuple[str, str]] = []
    for s in settings:
        info = parse_setting(s)
        sp_raw = info.get("sp", {}).get("raw", "")
        rel_raw = info.get("rel", {}).get("raw", "")
        key = (sp_raw, rel_raw)
        if key not in all_keys:
            all_keys.append(key)
    # Sort by (sp, rel) ascending numerically for stable linestyle assignment.
    all_keys.sort(key=lambda k: (int(k[0]) if k[0] else 0,
                                 int(k[1]) if k[1] else 0))

    out: dict[str, dict] = {}
    for s in settings:
        info = parse_setting(s)
        cf_raw = info.get("cf", {}).get("raw")
        sp_raw = info.get("sp", {}).get("raw", "")
        rel_raw = info.get("rel", {}).get("raw", "")
        color = _HUE_BY_CF.get(cf_raw, "#444444") if cf_raw else "#444444"
        try:
            idx = all_keys.index((sp_raw, rel_raw))
        except ValueError:
            idx = 0
        ls = _LINESTYLES[idx % len(_LINESTYLES)]
        out[s] = {"color": color, "linestyle": ls}
    return out


def build_label_map(settings: list[str]) -> dict[str, str]:
    return {s: format_label(s) for s in settings}


# ---- Discovery ------------------------------------------------------------
def discover_settings(results_dir: Path, fallback_dir: Path | None = None
                      ) -> list[str]:
    """List setting names from a results dir (subdirs); fall back to settings."""
    if results_dir.is_dir():
        names = sorted(d.name for d in results_dir.iterdir()
                       if d.is_dir() and not d.name.startswith("."))
        if names:
            return names
    if fallback_dir and fallback_dir.is_dir():
        return sorted(p.stem for p in fallback_dir.glob("*.toml"))
    return []
