"""This testing is to prove that the RV precision for condition 2 is incorrect for the paper.

Just taking the main components out of the code.
"""

import numpy as np


# This has the main bug.
def org_clump_tester(clump):
    """Return False when group of consecutive values in clump has a length >=3."""
    tester = True
    for block in clump:
        if len(clump) >= 3:  # clump should be block!
            tester = False
            break
    return tester


def corrected_clump_tester(clump):
    """Return False if a block in clump has a length >=3.

    The length of clump was used instead ot he length of block.
    """
    tester = True
    for block in clump:
        if len(block) >= 3:  # Fixed block!
            tester = False
            break
    return tester


def test_clump_tester():
    """Test old version and corrected version of clump_tester.

    This confirms the suspicions of the bugs.

    Note the sampling of the atmopsheric spectra is 10 so 1-2 zeros is
    unlikely to be picked up in a resolution element of spectra.
    Hence we choose a minimum of 3 consecutive zeros.
    """
    # Clump 1 does not have a group of consecutive Falses with len>=3 but the number of blocks in clump is >=3.
    clump1 = [np.array([0, 0], dtype=bool), np.array([0, 0], dtype=bool),
              np.array([0, 0], dtype=bool), np.array([0], dtype=bool)]
    # Clump1 should return True!
    assert org_clump_tester(clump1) is not True  # This is wrong
    assert corrected_clump_tester(clump1) is True

    # Clump 2 does have a group of consutive Falses>3 but not #blocks>3.
    clump2 = [np.array([0, 0, 0, 0, 0], dtype=bool), np.array([0, 0], dtype=bool)]
    # Clump2 should return False!
    assert org_clump_tester(clump2) is not False  # This is wrong
    assert corrected_clump_tester(clump2) is False

    # Test an array of zeros
    clump3 = [np.array([0, 0, 0, 0, 0], dtype=bool)]
    assert corrected_clump_tester(clump3) is False
    # Test a short array of ones
    clump4 = [np.array([0, 0], dtype=bool)]
    assert corrected_clump_tester(clump4) is True  # Not Long enough.


# This code snippet also has a bug.
def org_zeros_clumper(mask):
    """Purpose: Split into the groups of consecutive Falses."""
    mask_reversed = [not i for i in mask]
    clumps = np.array_split(mask, np.where(np.diff(mask_reversed))[0] + 1)[::2]
    return clumps


def correct_zeros_clumper(mask):
    """Purpose: Split into the groups of consecutive Falses."""
    mask_reversed = [not i for i in mask]
    # Depends on first value of mask
    if mask[0]:  # == True
        clumps = np.array_split(mask, np.where(np.diff(mask_reversed))[0] + 1)[1::2]
    else:
        clumps = np.array_split(mask, np.where(np.diff(mask_reversed))[0] + 1)[::2]

    return clumps


def test_zero_clumper():
    """Test clumping code. To confirm suspicions of the bugs."""
    mask0 = np.array([0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0], dtype=bool)
    expected_clump0 = [np.array([0], dtype=bool),
                       np.array([0, 0, 0], dtype=bool),
                       np.array([0, 0], dtype=bool),
                       np.array([0, 0, 0], dtype=bool),
                       np.array([0], dtype=bool)]

    clumps0 = org_zeros_clumper(mask0)
    clumps0_corr = correct_zeros_clumper(mask0)
    # For mask0 both methods work.
    for i, __ in enumerate(clumps0):
        assert (np.all(clumps0[i] == expected_clump0[i]))
        assert (np.all(clumps0_corr[i] == expected_clump0[i]))

    clumps0_inverted = org_zeros_clumper(~mask0)
    clumps0_corr_inverted = correct_zeros_clumper(~mask0)
    # Check number of clumps if inverted
    # Should change since mask0 has uneven group numbers. 5 False, 4 True.
    assert not (len(clumps0_inverted) != len(clumps0))  # This should not be the case
    assert (len(clumps0_corr_inverted) != len(clumps0_corr))

    # This mask starts with group of 1s so fails for original code.
    mask1 = np.array([1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0], dtype=bool)
    expected_clump1 = [np.array([0, 0, 0], dtype=bool),
                       np.array([0, 0], dtype=bool),
                       np.array([0, 0, 0], dtype=bool),
                       np.array([0], dtype=bool)]
    clumps1 = org_zeros_clumper(mask1)
    clumps1_corr = correct_zeros_clumper(mask1)

    for i, __ in enumerate(clumps1):
        assert not (np.all(clumps1[i] == expected_clump1[i]))  # Failed orginal case
        assert (np.all(clumps1_corr[i] == expected_clump1[i]))

    # Testing corner cases
    masked_zeros = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype=bool)
    expected_zeros = [masked_zeros]
    all_zeros_clumped = correct_zeros_clumper(masked_zeros)
    assert len(expected_zeros) == len(all_zeros_clumped)
    for i, __ in enumerate(all_zeros_clumped):
        assert np.all(all_zeros_clumped[i] == expected_zeros[i])

    masked_ones = np.array([1, 1, 1, 1, 1, 1], dtype=bool)
    expected_ones = []
    all_ones_clumped = correct_zeros_clumper(masked_ones)
    assert len(expected_ones) == len(all_ones_clumped)
    assert all_ones_clumped == []
