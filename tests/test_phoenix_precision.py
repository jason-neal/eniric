import pytest

from scripts.phoenix_precision import select_csv_columns, strip_whitespace


@pytest.mark.parametrize(
    "initial, expected",
    [
        (" test string goes here ", "teststringgoeshere"),
        (" 5500,  4.50, 0.0, K, 43.2322 , ", "5500,4.50,0.0,K,43.2322,"),
    ],
)
def test_strp_space(initial, expected):
    assert expected == strip_whitespace(initial)


def test_select_csv_columns():
    assert select_csv_columns("1,2,3,4,5,6,7,8,9,10") == "1,2,3,4,5,6,7,8"


@pytest.mark.parametrize(
    "initial, ncols, expected",
    [
        ("1,2,3,4,5,6,7", 4, "1,2,3,4"),
        ("hello,2,3,4,5,6,7,8,there, 10", 9, "hello,2,3,4,5,6,7,8,there"),
    ],
)
def test_select_csv_columns_with_ncols(initial, ncols, expected):
    assert select_csv_columns(initial, ncols=ncols) == expected
