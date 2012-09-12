# encoding: utf-8

"""Sequence AnaLyzeR.

pochisq & related functions adapted from http://www.fourmilab.ch/random/
"""

from math import fabs, exp, sqrt, log
from collections import defaultdict

__author__ = 'Eugen Dinca <decostin@gmail.com>'
__version__ = '0.1'


BIG_X = 20.0  # max value to represent exp (x)
Z_MAX = 6.0  # maximum meaningful z value
LOG_SQRT_PI = 0.5723649429247000870717135  # log(sqrt(pi))
I_SQRT_PI = 0.5641895835477562869480795  # 1/sqrt(pi)


def _poz(z):
    if z == 0.0:
        x = 0.0
    else:
        y = fabs(z) / 2.0
        if y >= Z_MAX / 2.0:
            x = 1.0
        elif y < 1.0:
            w = y * y
            x = ((((((((0.000124818987 * w
                        - 0.001075204047) * w + 0.005198775019) * w
                      - 0.019198292004) * w + 0.059054035642) * w
                    - 0.151968751364) * w + 0.319152932694) * w
                  - 0.531923007300) * w + 0.797884560593) * y * 2.0
        else:
            y -= 2.0
            x = (((((((((((((-0.000045255659 * y
                             + 0.000152529290) * y - 0.000019538132) * y
                           - 0.000676904986) * y + 0.001390604284) * y
                         - 0.000794620820) * y - 0.002034254874) * y
                       + 0.006549791214) * y - 0.010557625006) * y
                     + 0.011630447319) * y - 0.009279453341) * y
                   + 0.005353579108) * y - 0.002141268741) * y
                 + 0.000535310849) * y + 0.999936657524

    return (x + 1.0) / 2.0 if z > 0.0 else (1.0 - x) / 2.0


def _ex(x):
    return 0.0 if x < -BIG_X else exp(x)


def pochisq(x2, df):
    """
    Calculates the distribution function of the chi-squared distribution.

    >>> '{0:.4f}'.format(pochisq(2.688, 1))
    '0.1011'

    >>> '{0:.4f}'.format(pochisq(5.99, 2))
    '0.0500'

    >>> '{0:.4f}'.format(pochisq(26.76, 11))
    '0.0050'

    >>> '{0:.4f}'.format(pochisq(8, 15))
    '0.9238'
    """
    if x2 <= 0.0 or df < 1:
        return 1.0

    a = x2 / 2.0
    even = not bool(df & 1)

    y = 0.0
    if df > 1:
        y = _ex(-a)
    s = y if even else 2.0 * _poz(-sqrt(x2))

    if df > 2:
        x2 = (df - 1.0) / 2.0
        z = 1.0 if even else 0.5
        if a > BIG_X:
            e = 0.0 if even else LOG_SQRT_PI
            c = log(a)
            while z <= x2:
                e += log(z)
                s += _ex(c * z - a - e)
                z += 1.0

            return s
        else:
            e = 1.0 if even else I_SQRT_PI / sqrt(a)
            c = 0.0
            while z <= x2:
                e *= a / z
                c += e
                z += 1.0

            return c * y + s

    return s


def build_histograms(sequences_list):
    """
    Builds a columnwise histogram for all sequences in sequences_list.

    >>> sequences_list = ['abc', 'bca']
    >>> expected_histograms = [{'a': 1, 'b': 1},
    ...     {'b': 1, 'c': 1}, {'a': 1, 'c': 1}]
    >>> histograms = build_histograms(sequences_list)
    >>> histograms == expected_histograms
    True
    """
    # we're assumming all sequences to be of same length
    sequence_len = len(sequences_list[0])
    # construct the zero count dictionary for all possible values in the range
    histograms = [defaultdict(int) for i in range(sequence_len)]

    # go through all sequences, one vertical slice at a time
    v_slices = zip(*sequences_list)
    for i in range(len(v_slices)):
        for char in v_slices[i]:
            histograms[i][char] += 1

    return histograms


def character_count_tester(sequences_list):
    """
    Returns a (probability observed distribution is uniform, histogram) tuple,
    for each column calculated across all sequences in sequences_list.
    It assumes all sequences have the same length.

    >>> sequences_list = ['abc', 'bca']
    >>> expected_test_results = [(1.0, {'a': 1, 'b': 1}),
    ...     (1.0, {'b': 1, 'c': 1}), (1.0, {'a': 1, 'c': 1})]
    >>> test_results = character_count_tester(sequences_list)
    >>> test_results == expected_test_results
    True
    """
    # build the histograms, one for each position
    histograms = build_histograms(sequences_list)
    nb_sequences = len(sequences_list)

    test_results = []

    # go over each position in the sequence & compute the chi-squared value
    for histogram in histograms:
        # how many chars were found in the all sequences at this position
        char_set_size = len(histogram)

        df = char_set_size - 1  # the degrees of freedom
        expected_freq = 1.0 * nb_sequences / \
            char_set_size  # assume an uniform distribution

        d2e = [(freq - expected_freq) * (freq - expected_freq) /
               expected_freq for freq in histogram.values()]
        x2 = sum(d2e)
        test_results.append((pochisq(x2, df), histogram))

    return test_results

if __name__ == '__main__':
    import doctest
    doctest.testmod()
