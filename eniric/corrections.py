"""Module with functions that have precision corrections."""
from typing import Dict


def correct_artigau_2018(band: str) -> float:
    """Apply Artigau 2018 nominal Barnard's Star corrections.

    Table 3: Multiplicative correction factors to be applied
    on the RV precision derived from stellar models.
    These values correspond to the square-root of
    the flux-weighted mean Q ratio between observation
    and models for each bandpass. The nominal values are
    for a comparison with the default model described here,
    but we also explore the impact of other physical parameter
    choices.

             |Nominal|  Variants
    [Fe/H]   |−0.5   |  0.0
    log g    | 5.0   |        5.5
    Teff (K) | 3200  |              3400    3000
    ----------------------------------------
    g        |  0.66 | 0.63  0.69  0.82   0.51
    r        |  0.82 | 0.76  0.84  1.08   0.60
    i        |  0.94 | 0.81  0.99  1.26   0.73
    z        |  1.27 | 1.09  1.20  1.82   0.96
    Y        |  0.29 | 0.30  0.27  0.30   0.25
    J        |  0.38 | 0.40  0.31  0.54   0.37
    H        |  1.37 | 0.95  1.82  1.42   1.23
    K        |  1.47 | 1.06  2.00  1.66   1.27

    Usage:
      Multiply by model precisions to get to real corrections.
    """

    corrections = {"g": 0.66, "r": 0.82, "i": 0.94, "Z": 1.27,
                   "z": 1.27, "Y": 0.29, "J": 0.38, "H": 1.37, "K": 1.47}  # type: Dict[str, float]
    try:
        return corrections[band]
    except KeyError as e:
        raise KeyError("Band '{}' does not have a valid correction key".format(band))
