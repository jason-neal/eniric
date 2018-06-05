import pytest

from eniric.corrections import correct_artigau_2018


@pytest.mark.parametrize("band", ["g", "Z", "z", "K", "H", "J", "Y"])
def test_bands_in_correction(band):
    value = correct_artigau_2018(band)
    assert isinstance(value, float)
    assert value > 0.2
    assert value < 1.5


@pytest.mark.parametrize("band", ["x", "p", "a"])
def test_bands_not_in_correction(band):
    with pytest.raises(KeyError):
        correct_artigau_2018(band)


@pytest.mark.parametrize("band, value", [("Y", 0.29), ("K", 1.47)])
def test_corection_specific_values(band, value):
    assert correct_artigau_2018(band) == value
