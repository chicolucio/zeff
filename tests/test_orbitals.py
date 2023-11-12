import pytest
from sqlalchemy.exc import NoResultFound

from zeff import OrbitalsInfo, orbitals


def test_valid_element():
    # Test with a valid element name
    result = orbitals("Hydrogen")
    assert isinstance(result, OrbitalsInfo)
    assert result.orbital == ["1s"]
    assert result.principal_quantum_number == [1]
    assert result.orbital_letter == ["s"]
    assert result.azimuthal_quantum_number == [0]


def test_valid_symbol():
    # Test with a valid element symbol
    result = orbitals("H")
    assert isinstance(result, OrbitalsInfo)
    assert result.orbital == ["1s"]
    assert result.principal_quantum_number == [1]
    assert result.orbital_letter == ["s"]
    assert result.azimuthal_quantum_number == [0]


def test_invalid_element():
    # Test with an invalid element name
    with pytest.raises(NoResultFound):
        orbitals("Unobtainium")


def test_invalid_symbol():
    # Test with an invalid element symbol
    with pytest.raises(NoResultFound):
        orbitals("Uo")


def test_complex_element():
    # Test with a complex element that has multiple orbitals
    result = orbitals("Carbon")
    assert isinstance(result, OrbitalsInfo)
    assert result.orbital == ["1s", "2s", "2p"]
    assert result.principal_quantum_number == [1, 2, 2]
    assert result.orbital_letter == ["s", "s", "p"]
    assert result.azimuthal_quantum_number == [0, 0, 1]
