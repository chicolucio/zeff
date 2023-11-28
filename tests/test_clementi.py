import numpy as np

from zeff import ClementiInfo, clementi


def test_hydrogen():
    # Test with a valid simple element name
    result = clementi("Hydrogen")
    assert isinstance(result, ClementiInfo)
    assert len(result.effective_nuclear_charge) > 0
    assert len(result.screening_values) > 0
    assert len(result.screening_percentage) > 0
    assert np.allclose(result.effective_nuclear_charge, [1])
    assert np.allclose(result.screening_values, [0])
    assert np.allclose(result.screening_percentage, [0])


def test_element_above_86():
    # Test with an element where atomic number is greater than 86
    result = clementi("Uranium")  # Uranium has an atomic number of 92
    assert np.all(np.isnan(result.effective_nuclear_charge))
    assert np.all(np.isnan(result.screening_values))
    assert np.all(np.isnan(result.screening_percentage))


def test_nitrogen():
    # Test with a complex element that has multiple orbitals
    result = clementi("N")
    assert isinstance(result, ClementiInfo)
    assert len(result.effective_nuclear_charge) == 3
    assert len(result.screening_values) == 3
    assert len(result.screening_percentage) == 3
    assert np.allclose(
        result.effective_nuclear_charge, [6.665, 3.847, 3.834], atol=1e-3
    )
    assert np.allclose(result.screening_values, [0.3349, 3.1526, 3.1660], atol=1e-3)
    assert np.allclose(result.screening_percentage, [4.784, 45.037, 45.228], atol=1e-3)
