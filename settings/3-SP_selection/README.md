# Experiment 3 — Subproblem Selection

Compares strategies for the partial subproblem solving loop.

## Grid Search Stages
- `1-n_hits_rel/*.toml`: RND + n hits rel (0.01, 0.05, 0.1, 0.25, 0.5, 1.0)
  -> best rel_frac
- `2-score-individual/*.toml`: test SCORE: viol 1.0, reliability 1.0, cut_share 1.0 with rel_frac (with or without decay `BENDERS.PARTIAL.SCORE_DECAY`)
  -> best two scores, with or without decay
- `3-score-hybrid/*.toml`: test SCORE: best two scores with or without decay ((0.9, 0.1), (0.75, 0.25), (0.5, 0.5), (0.25, 0.75), (0.1, 0.9)) (and we already have (1.0, 0.0) and (0.0, 1.0) from the previous stage)
  -> best score combination, with or without decay (score*)
- `3-score-hybrid-random/*.toml`: SCORE hybrid combining `reliability` with the random exploration component (`random`) at the same weight pairs ((0.9, 0.1), (0.75, 0.25), (0.5, 0.5), (0.25, 0.75), (0.1, 0.9)) (and we already have (1.0, 0.0) for reliability from Stage 2 and (0.0, 1.0) is equivalent to `SORTING_MODE = ":RND"`)
  -> best reliability/random mix (does adding exploration help over pure score?)
- `4-score-decay/*.toml`: if with decay performs better, test decay weights (0.001, 0.1, 0.25, 0.5, 0.75) (and we already have 0.9 from the previous stage)
  -> best decay weight (add to score*)
- `5-score-minprio/*.toml`: SCORE + minprio limit (0.25, 0.5, 0.75, 0.9)
  -> best minprio cutoff for minprio setting

## Full Test Variants
- `rnd_pb.toml`: RND + PB
- `rnd_nhits_rel.toml`: RND + n_hits_rel (`STOP_AFTER_N_HITS_REL` on best rel_frac from Stage 1)
- after grid searches:
  - set score to score*, i.e., best score combination with or without decay
  - `score_pb.toml`: PB (`BENDERS.PARTIAL.LIMITS.STOP_ON_PB` on)
  - `score_minprio.toml`: minprio (`BENDERS.PARTIAL.LIMITS.MINPRIO` on best cutoff from grid search)
  - `score_nhits_rel.toml`: n_hits_rel (`STOP_AFTER_N_HITS_REL` on best rel_frac from grid search)

## Carrier (carried from prior experiments)
- `BENDERS.SUBPROBLEM_SOLVER_TYPE = "LinEqSysRfSolver"` (winner of Exp. 2)
- no cut filtering (to avoid any interference with the subproblem selection strategy)
- if any score is applied, apply staleness with weight 0.01 (we just think this is a good idea)

Both are already `default.toml`; no override needed.

