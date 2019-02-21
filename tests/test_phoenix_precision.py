import pytest

from scripts.phoenix_precision import select_csv_columns, strip_whitespace, check_model


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


@pytest.mark.parametrize(
    "model_in, expected_model",
    [
        ("aces", "aces"),
        ("phoenix", "aces"),
        ("btsettl", "btsettl"),
    ],
)
def test_check_model(model_in, expected_model):
    model_out = check_model(model_in)
    assert model_out == expected_model


def test_model_depreciation():
    with pytest.deprecated_call():
        check_model("phoenix")


@pytest.mark.parametrize(
    "invalid_model",
    ["phoneix-aces", "best", 42],
)
def test_check_model_error(invalid_model):
    with pytest.raises(ValueError):
        check_model(invalid_model)
