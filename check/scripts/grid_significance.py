"""One-sided Wilcoxon comparisons of paired grid solve times.

SciPy's signed-rank test uses configuration minus control solve times and the
alternative 'less'. Holm adjustment retains the full planned family, as in
R's p.adjust(method='holm', n=...).
"""
import math


def one_sided_wilcoxon_test(candidate, control):
    from scipy.stats import wilcoxon

    if not candidate or len(candidate) != len(control):
        raise ValueError('Require equally sized nonempty paired samples')
    if not all(math.isfinite(t) and t >= 0 for t in candidate + control):
        raise ValueError('Require finite nonnegative solve times')
    differences = [a - b for a, b in zip(candidate, control)]
    nonzero = [d for d in differences if d != 0]
    ties = len(nonzero) - len(set(abs(d) for d in nonzero))
    if not nonzero:
        return dict(p_value=1.0, wilcoxon_statistic=0.0, wilcoxon_method='all pairs equal',
                    nonzero_pairs=0, absolute_ties=0)
    # SciPy's exact distribution requires distinct nonzero absolute differences.
    method = 'exact' if len(nonzero) == len(differences) and not ties else 'approx'
    result = wilcoxon(differences, alternative='less', zero_method='wilcox',
                      correction=False, method=method)
    return dict(p_value=float(result.pvalue), wilcoxon_statistic=float(result.statistic),
                wilcoxon_method=method, nonzero_pairs=len(nonzero), absolute_ties=ties)


def holm_adjust(pvalues, family_size):
    if family_size < len(pvalues):
        raise ValueError('Family size cannot be smaller than available tests')
    if not all(math.isfinite(p) and 0 <= p <= 1 for p in pvalues.values()):
        raise ValueError('Require finite p-values in [0, 1]')
    adjusted, previous = {}, 0.0
    for rank, (name, p) in enumerate(sorted(pvalues.items(), key=lambda item: (item[1], item[0]))):
        previous = max(previous, min(1.0, (family_size - rank) * p))
        adjusted[name] = previous
    return adjusted


def significance_stars(p):
    if p is None:
        return ''
    return '***' if p < .001 else ('**' if p < .01 else ('*' if p < .05 else ''))


def ratio_label(cell, key='ratio'):
    label = f"{cell[key]:.2f}"
    stars = cell.get('significance_stars', '')
    return label + (rf'\textsuperscript{{{stars}}}' if stars else '')
