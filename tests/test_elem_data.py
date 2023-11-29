import pandas as pd

from zeff import elem_data


def test_valid_element():
    # Test with a valid element name
    result = elem_data("Hydrogen")
    expected_columns = [
        "n",
        "l",
        "l_num",
        "Orbital",
        "Zeff Slater",
        "S Slater",
        "% S Slater",
        "Zeff Clementi",
        "S Clementi",
        "% S Clementi",
    ]

    assert isinstance(result, pd.DataFrame)
    assert list(result.columns) == expected_columns
    assert len(result) > 0  # Ensure DataFrame is not empty


def test_dataframe_structure():
    # Test the structure of the DataFrame for a known element
    result = elem_data("Carbon")
    assert all(
        isinstance(col, str) for col in result.columns
    )  # Check if all column names are strings
    assert all(
        result[col].dtype for col in result.columns
    )  # Check if all columns have a dtype
    assert result.shape == (3, 10)  # Check if the DataFrame has the correct shape
