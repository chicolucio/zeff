import numpy as np

from zeff import SlaterInfo, slater


def test_hydrogen():
    # Test with a valid simple element name
    result = slater("Hydrogen")
    assert isinstance(result, SlaterInfo)
    assert len(result.effective_nuclear_charge) > 0
    assert len(result.screening_values) > 0
    assert len(result.screening_percentage) > 0
    assert np.allclose(result.effective_nuclear_charge, [1])
    assert np.allclose(result.screening_values, [0])
    assert np.allclose(result.screening_percentage, [0])


def test_nitrogen():
    # Test with a complex element that has multiple orbitals
    result = slater("N")
    assert isinstance(result, SlaterInfo)
    assert len(result.effective_nuclear_charge) == 3
    assert len(result.screening_values) == 3
    assert len(result.screening_percentage) == 3
    assert np.allclose(result.effective_nuclear_charge, [6.7, 3.9, 3.9], atol=1e-2)
    assert np.allclose(result.screening_values, [0.3, 3.1, 3.1], atol=1e-2)
    assert np.allclose(result.screening_percentage, [4.28, 44.28, 44.28], atol=1e-2)
