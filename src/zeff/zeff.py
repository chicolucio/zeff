from pathlib import Path
from typing import List, NamedTuple, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.axes import Axes
from mendeleev import element
from sqlalchemy.exc import NoResultFound

script_path = Path(__file__).parent.absolute()

plt.style.use(str("/".join([str(script_path), "plots.mplstyle"])))


class OrbitalsInfo(NamedTuple):
    orbital: List[str]
    principal_quantum_number: List[int]
    orbital_letter: List[str]
    azimuthal_quantum_number: List[int]


class SlaterInfo(NamedTuple):
    effective_nuclear_charge: list
    screening_values: list
    screening_percentage: list


class ClementiInfo(NamedTuple):
    effective_nuclear_charge: list
    screening_values: list
    screening_percentage: list


def orbitals(name: str) -> OrbitalsInfo:
    """Extracts n, l, and nl notation from electronic configuration info of
    the Mendeleev package.

    Parameters
    ----------
    name : str
        Name or symbol for the chemical element.

    Returns
    -------
    OrbitalsInfo
        A NamedTuple containing lists with orbital representation (nl notation),
        principal quantum number (n), orbital letter notation (l), and
        azimuthal quantum number.

    Raises
    ------
    NoResultFound
        If the provided element name or symbol is not valid.
    """
    try:
        elem = element(name)
    except NoResultFound as e:
        raise NoResultFound(f"Invalid element name or symbol: {name}") from e

    dict_l = {"s": 0, "p": 1, "d": 2, "f": 3}
    orbital, orbital_n, orbital_l, orbital_l_num = [], [], [], []

    for (n, l), _ in elem.ec.conf.items():
        orbital_n.append(n)
        orbital_l.append(l)
        orbital_l_num.append(dict_l.get(l, -1))
        orbital.append(f"{n}{l}")

    return OrbitalsInfo(orbital, orbital_n, orbital_l, orbital_l_num)


def slater(name: str) -> SlaterInfo:
    """Calculates Slater's screening values, effective nuclear charge, and
    screening percentages for each orbital of an element's electronic configuration.

    Parameters
    ----------
    name : str
        Name or symbol for the chemical element.

    Returns
    -------
    SlaterInfo
        A NamedTuple containing lists of effective nuclear charge values,
        screening values, and screening percentages for each orbital.
    """
    orbital_info = orbitals(name)
    elem = element(name)

    zeff_slater = []
    slater_s = []

    for n, l in zip(orbital_info.principal_quantum_number, orbital_info.orbital_letter):
        screening = elem.ec.slater_screening(n, l)
        slater_s.append(screening)
        zeff_slater.append(elem.atomic_number - screening)

    slater_s_percent = [(i / elem.atomic_number) * 100 for i in slater_s]

    return SlaterInfo(zeff_slater, slater_s, slater_s_percent)


def clementi(name: str) -> ClementiInfo:
    """Calculates screening and effective nuclear charge based on values
    calculated by Clementi and Raimond. No values for Z > 86.

    Parameters
    ----------
    name : str
        Name or symbol for the chemical element.

    Returns
    -------
    ClementiInfo
        A NamedTuple containing lists of effective nuclear charge values,
        screening values, and screening percentages for each orbital.
    """

    orbital_info = orbitals(name)
    elem = element(name)
    zeff_clementi = []
    clementi_s = []

    if elem.atomic_number < 87:
        for n, l in zip(
            orbital_info.principal_quantum_number, orbital_info.orbital_letter
        ):
            screening = elem.atomic_number - elem.zeff(n, l, method="clementi")
            clementi_s.append(screening)
            zeff_clementi.append(elem.zeff(n, l, method="clementi"))
    else:
        for _ in orbital_info.principal_quantum_number:
            clementi_s.append(np.nan)
            zeff_clementi.append(np.nan)

    clementi_s_percent = [(i / elem.atomic_number) * 100 for i in clementi_s]

    return ClementiInfo(zeff_clementi, clementi_s, clementi_s_percent)


def elem_data(name: str) -> pd.DataFrame:
    """Data summary for a given element.

    Parameters
    ----------
    name : str
        Name or symbol for the chemical element.

    Returns
    -------
    pd.DataFrame
        Summary of Zeff and S for a given element with Slater and Clementi
        screening values.
    """
    slater_info = slater(name)
    clementi_info = clementi(name)
    orbital_info = orbitals(name)

    data = pd.DataFrame(
        {
            "n": orbital_info.principal_quantum_number,
            "l": orbital_info.orbital_letter,
            "l_num": orbital_info.azimuthal_quantum_number,
            "Orbital": orbital_info.orbital,
            "Zeff Slater": slater_info.effective_nuclear_charge,
            "S Slater": slater_info.screening_values,
            "% S Slater": slater_info.screening_percentage,
            "Zeff Clementi": clementi_info.effective_nuclear_charge,
            "S Clementi": clementi_info.screening_values,
            "% S Clementi": clementi_info.screening_percentage,
        }
    )

    data_org = data.sort_values(by=["n", "l_num"]).reset_index(drop=True)

    return data_org


def plot_slater_zeff(name: str, ax=None, **kwargs) -> Axes:
    """Plots effective nuclear charge (Slater) per orbital.

    Parameters
    ----------
    name : str
        String with the element symbol.
    ax : matplotlib.axes._subplots.AxesSubplot, optional
        Axes for the plot. If None, a new figure and axes object is created.
    **kwargs
        Additional keyword arguments for the plot.

    Returns
    -------
    matplotlib.axes._subplots.AxesSubplot
        The Axes object with the plot.
    """

    data = elem_data(name)

    if ax is None:
        fig, ax = plt.subplots()

    x = data["Orbital"]
    y = data["Zeff Slater"]

    ax.set_ylabel("Effective nuclear charge")
    ax.set_xlabel("Orbitals")
    ax.set_xticks(range(len(data)))
    ax.set_yticks(np.arange(0, round(max(y)) + 5, 5.0))

    ax.plot(x, y, label=f"Zeff Slater {element(name).symbol}", **kwargs)
    ax.legend(
        loc="best",
        shadow=False,
        fancybox=True,
        bbox_to_anchor=(1, 1),
    )

    return ax


def plot_slater_screening(name: str, ax=None, **kwargs) -> Axes:
    """Plots Shielding /% (Slater) per orbital.

    Parameters
    ----------
    name : str
        String with the element symbol.
    ax : matplotlib.axes._subplots.AxesSubplot, optional
        Axes for the plot. If None, a new figure and axes object is created.
    **kwargs
        Additional keyword arguments for the plot.

    Returns
    -------
    matplotlib.axes._subplots.AxesSubplot
        The Axes object with the plot.
    """

    data = elem_data(name)

    if ax is None:
        fig, ax = plt.subplots()

    x = data["Orbital"]
    y = data["% S Slater"]

    ax.set_ylabel("Shielding / %")
    ax.set_xlabel("Orbitals")
    ax.set_xticks(range(len(data)))
    ax.set_yticks(np.arange(0, round(max(y)) + 5, 10.0))

    ax.plot(x, y, label=f"Shielding % Slater {element(name).symbol}", **kwargs)
    ax.legend(loc="best", shadow=True, fancybox=True)

    return ax


def plot_clementi_zeff(name: str, ax=None, **kwargs) -> Axes:
    """Plots effective nuclear charge (Clementi) per orbital.

    Parameters
    ----------
    name : str
        String with the element symbol.
    ax : matplotlib.axes._subplots.AxesSubplot, optional
        Axes for the plot. If None, a new figure and axes object is created.
    **kwargs
        Additional keyword arguments for the plot.

    Returns
    -------
    matplotlib.axes._subplots.AxesSubplot
        The Axes object with the plot.
    """

    data = elem_data(name)

    if ax is None:
        fig, ax = plt.subplots()

    x = data["Orbital"]
    y = data["Zeff Clementi"]

    ax.set_ylabel("Effective nuclear charge")
    ax.set_xlabel("Orbitals")
    ax.set_xticks(range(len(data)))
    ax.set_yticks(np.arange(0, round(max(y)) + 5, 5.0))

    ax.plot(x, y, label=f"Zef Clementi {element(name).symbol}", **kwargs)
    ax.legend(loc="best", shadow=False, fancybox=True, bbox_to_anchor=(1, 1))

    return ax


def plot_clementi_screening(name: str, ax=None, **kwargs) -> Axes:
    """Plots Shielding /% (Clementi) per orbital.

    Parameters
    ----------
    name : str
        String with the element symbol.
    ax : matplotlib.axes._subplots.AxesSubplot, optional
        Axes for the plot. If None, a new figure and axes object is created.
    **kwargs
        Additional keyword arguments for the plot.

    Returns
    -------
    matplotlib.axes._subplots.AxesSubplot
        The Axes object with the plot.
    """

    data = elem_data(name)

    if ax is None:
        fig, ax = plt.subplots()

    x = data["Orbital"]
    y = data["% S Clementi"]

    ax.set_ylabel("Shielding / %")
    ax.set_xlabel("Orbitals")
    ax.set_xticks(range(len(data)))
    ax.set_yticks(np.arange(0, round(max(y)) + 5, 10.0))

    ax.plot(x, y, label=f"Shielding % Clementi {element(name).symbol}", **kwargs)
    ax.legend(loc="best", shadow=True, fancybox=True)

    return ax


def plot_slater_both(name: str, ax=None) -> Tuple[Axes, Axes]:
    """Plots Effective nuclear charge and shielding /% (Slater) per orbital.

    Parameters
    ----------
    name : str
        String with the element symbol.
    ax : matplotlib.axes._subplots.AxesSubplot, optional
        Axes for the plot. If None, a new figure and axes object is created.

    Returns
    -------
    Tuple[matplotlib.axes._subplots.AxesSubplot, matplotlib.axes._subplots.AxesSubplot]
        The Axes objects with the plot.
    """
    data = elem_data(name)

    if ax is None:
        fig, ax = plt.subplots()

    x = data["Orbital"]
    y = data["Zeff Slater"]
    ax.set_ylabel("Effective nuclear charge", color="C0")
    ax.set_xlabel("Orbitals")
    line1 = ax.plot(x, y, label=f"Zef Slater {element(name).symbol}", color="C0")
    for t in ax.get_yticklabels():
        t.set_color("C0")
    ax.set_xticks(range(len(data)))
    ax.set_yticks(np.linspace(0, round(max(y)) + 1, 5))

    ax2 = ax.twinx()
    y2 = data["% S Slater"]
    line2 = ax2.plot(
        x,
        y2,
        label=f"Shielding % Slater {element(name).symbol}",
        color="C1",
    )
    for t in ax2.get_yticklabels():
        t.set_color("C1")
    ax2.set_ylabel("Shielding / %", color="C1")
    ax2.set_yticks(np.linspace(0, round(max(y2)) + 1, 5))

    lines = line1 + line2
    labels = [line.get_label() for line in lines]

    ax2.legend(lines, labels, loc="upper center", shadow=True, fancybox=True)

    # keep only the vertical grids
    ax.grid(visible=False)
    ax2.grid(visible=False)
    ax.grid(visible=True, which="major", axis="x")

    return ax, ax2


def plot_clementi_both(name: str, ax=None) -> Tuple[Axes, Axes]:
    """Plots Effective nuclear charge and shielding /% (Clementi) per orbital.

    Parameters
    ----------
    name : str
        String with the element symbol.
    ax : matplotlib.axes._subplots.AxesSubplot, optional
        Axes for the plot. If None, a new figure and axes object is created.

    Returns
    -------
    Tuple[matplotlib.axes._subplots.AxesSubplot, matplotlib.axes._subplots.AxesSubplot]
        The Axes objects with the plot.
    """
    data = elem_data(name)

    if ax is None:
        fig, ax = plt.subplots()

    x = data["Orbital"]
    y = data["Zeff Clementi"]
    ax.set_ylabel("Effective nuclear charge", color="C0")
    ax.set_xlabel("Orbitals")
    line1 = ax.plot(x, y, label=f"Zef Clementi {element(name).symbol}", color="C0")
    for t in ax.get_yticklabels():
        t.set_color("C0")
    ax.set_xticks(range(len(data)))
    ax.set_yticks(np.linspace(0, round(max(y)) + 1, 5))

    ax2 = ax.twinx()
    y2 = data["% S Clementi"]
    line2 = ax2.plot(
        x,
        y2,
        label=f"Shielding % Clementi {element(name).symbol}",
        color="C1",
    )
    for t in ax2.get_yticklabels():
        t.set_color("C1")
    ax2.set_ylabel("Shielding / %", color="C1")
    ax2.set_yticks(np.linspace(0, round(max(y2)) + 1, 5))

    lines = line1 + line2
    labels = [line.get_label() for line in lines]

    ax2.legend(lines, labels, loc="upper center", shadow=True, fancybox=True)

    # keep only the vertical grids
    ax.grid(visible=False)
    ax2.grid(visible=False)
    ax.grid(visible=True, which="major", axis="x")

    return ax, ax2
