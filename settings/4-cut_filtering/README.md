# Experiment 4 — Cut Filtering

## Grid Search Stages
- Template inputs for the Stage 1 grids live in `cutfiltering-gridsearch-template-base-settings/` and `cutfiltering-gridsearch-template-params/`.
- Generated Stage 1 setting matrices live in:
  `1-cutfiltering-pb-spsel-rel/`, `2-cutfiltering-pb-spsel-abs/`, `3-cutfiltering-abs-spsel-abs/`,
  `4-cutfiltering-rel-spsel-rel/`, `5-cutfiltering-rel-spsel-abs/`, `6-cutfiltering-abs-spsel-rel/`.
- For `pb` cut filtering, only the `sp` factor varies.
- Relative components always use three-digit filename tokens: `cfXXX` means `LIMIT_MODE = ":rel_num"` with `LIMIT_REL_ALPHA = XXX / 100`, and `spXXX` means `STOP_AFTER_N_HITS_REL = XXX / 100`.
- Absolute components always use two-digit filename tokens: `cfXX` means `LIMIT_MODE = ":abs_num"` with `LIMIT_ABS_NUM = XX`, and `spXX` means `STOP_AFTER_N_HITS = XX`.
- `abs` components use the dedicated count grid `{1,2,5,10}`. Mixed folders therefore use mixed-width names such as `cf001_sp05` and `cf02_sp005`.
- `5-cutfiltering-rel-spsel-abs/` uses the full relative CF range `{001,005,010,025,050,075,100}` crossed with the absolute SP range `{01,02,05,10}`.
- `6-cutfiltering-abs-spsel-rel/` uses the full absolute CF range `{01,02,05,10}` crossed with the relative SP range `{001,005,010,025,050,075,100}`.
- `1-cf-spsel-hybrid/*.toml`: test cut filtering plus SP selection together: (cut filtering, sp selection)
  (0.01, 1.00), (0.01, 0.75), (0.01, 0.50), (0.01, 0.25), (0.01, 0.10), (0.01, 0.05)
  (0.05, 1.00), (0.05, 0.75), (0.05, 0.50), (0.05, 0.25), (0.05, 0.10)
  (0.10, 1.00), (0.10, 0.75), (0.10, 0.50), (0.10, 0.25)
  (0.25, 1.00), (0.25, 0.75), (0.25, 0.50)
  (0.50, 1.00), (0.50, 0.75)
  (0.75, 1.00), 
  (1.00, 1.00) (baseline)
   -> best combination of cut filtering and SP selection (cutfiltering_sp_selection.toml)
- `2-cf-spsel-hybrid-viol/*.toml` same as above, but instead of diversity hybrid, use violation only

## Full Test Variants
- `cutfiltering.toml`: best cut filtering strategy from grid search 0.1 rel_num (not rel_num_iter)
- `cutfiltering_sp_selection.toml`: best cut filtering strategy + best subproblem selection strategy 
- for now, don't test PB or minprio in combination with cut filtering

## Carrier (carried from prior experiments)
- as always, LEQRF
- best hybrid score and decay setting from 3-SP_selection
- cut filtering settings always:
  - :DIVERSITY_AGG_HYB
  - :rel mode,
  - vary NUM_CLUSTERS_REL_ALPHA for the percentage of cuts to add,
  - do k_plus_1 mode and max_violation representative, and 
  - weighted_agg and weighted_agg_normalize on.