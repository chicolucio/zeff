from typing import List, NamedTuple

import numpy as np
import pandas as pd
from mendeleev import element
from sqlalchemy.exc import NoResultFound


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
            "Zef Slater": slater_info.effective_nuclear_charge,
            "S Slater": slater_info.screening_values,
            "% S Slater": slater_info.screening_percentage,
            "Zef Clementi": clementi_info.effective_nuclear_charge,
            "S Clementi": clementi_info.screening_values,
            "% S Clementi": clementi_info.screening_percentage,
        }
    )

    data_org = data.sort_values(by=["n", "l_num"]).reset_index(drop=True)

    return data_org


def plot_param(ax=None):
    """Common plot parameters.

    Parameters
    ----------
    ax : matplotlib.axes
        Axes for the plot (the default is None).

    Returns
    -------
    None
        Should be called inside a plot function in order to apply the
        parameters.

    """
    ax.axhline(color="gray", zorder=-1)
    ax.axvline(color="gray", zorder=-1)
    # ax.grid(b=True, which="major", linestyle=":", linewidth=2)
    ax.grid(which="major", linestyle=":", linewidth=2)
    ax.minorticks_on()
    # ax.grid(b=True, which="minor", axis="y", linestyle=":", linewidth=1.0)
    ax.grid(which="minor", axis="y", linestyle=":", linewidth=1.0)
    ax.tick_params(which="both", labelsize=14)


def plot_slater_zef(name, ax=None, **kwargs):
    """Plots effective nuclear charge (Slater) per orbital.

    Parameters
    ----------
    name : string
        String with the element symbol.
    ax : matplotlib.axes
        Axes for the plot (the default is None).
    **kwargs : string
        kwargs for plot parameters.

    Returns
    -------
    None
        Plot.

    """
    plot_param(ax)
    fontsize = 15
    linewidth = 4
    x = elem_data(name)["Orbital"]
    y = elem_data(name)["Zef Slater"]
    ax.set_ylabel("Effective nuclear charge", size=fontsize)
    ax.set_xlabel("Orbitals", size=fontsize)
    ax.set_xticks = [i for i in range(len(elem_data(name).index))]
    ax.set_yticks(np.arange(0, round(max(y)) + 5, 5.0))
    ax.plot(
        x,
        y,
        linewidth=linewidth,
        label=f"Zef Slater {element(name).symbol}",
        **kwargs,
    )
    ax.legend(
        fontsize=fontsize - 1,
        loc="best",
        shadow=False,
        fancybox=True,
        bbox_to_anchor=(1, 1),
    )


def plot_slater_screening(name, ax=None, **kwargs):
    """Plots Shielding /% (Slater) per orbital.

    Parameters
    ----------
    name : string
        String with the element symbol.
    ax : matplotlib.axes
        Axes for the plot (the default is None).
    **kwargs : string
        kwargs for plot parameters.

    Returns
    -------
    None
        Plot.

    """
    plot_param(ax)
    fontsize = 15
    linewidth = 4
    x = elem_data(name)["Orbital"]
    y = elem_data(name)["% S Slater"]
    ax.set_ylabel("Shielding / %", size=fontsize)
    ax.set_xlabel("Orbitals", size=fontsize)
    ax.set_xticks = [i for i in range(len(elem_data(name).index))]
    ax.set_yticks(np.arange(0, round(max(y)) + 5, 10.0))
    ax.plot(
        x,
        y,
        linewidth=linewidth,
        label=f"Shielding % Slater {element(name).symbol}",
        **kwargs,
    )
    ax.legend(fontsize=fontsize - 2, loc="best", shadow=True, fancybox=True)


def plot_clementi_zef(name, ax=None, **kwargs):
    """Plots effective nuclear charge (Clementi) per orbital.

    Parameters
    ----------
    name : string
        String with the element symbol.
    ax : matplotlib.axes
        Axes for the plot (the default is None).
    **kwargs : string
        kwargs for plot parameters.

    Returns
    -------
    None
        Plot.

    """
    plot_param(ax)
    fontsize = 15
    linewidth = 4
    x = elem_data(name)["Orbital"]
    y = elem_data(name)["Zef Clementi"]
    ax.set_ylabel("Effective nuclear charge", size=fontsize)
    ax.set_xlabel("Orbitals", size=fontsize)
    ax.set_xticks = [i for i in range(len(elem_data(name).index))]
    ax.set_yticks(np.arange(0, round(max(y)) + 5, 5.0))
    ax.plot(
        x,
        y,
        linewidth=linewidth,
        label=f"Zef Clementi {element(name).symbol}",
        **kwargs,
    )
    ax.legend(
        fontsize=fontsize - 1,
        loc="best",
        shadow=False,
        fancybox=True,
        bbox_to_anchor=(1, 1),
    )


def plot_clementi_screening(name, ax=None, **kwargs):
    """Plots Shielding /% (Clementi) per orbital.

    Parameters
    ----------
    name : string
        String with the element symbol.
    ax : matplotlib.axes
        Axes for the plot (the default is None).
    **kwargs : string
        kwargs for plot parameters.

    Returns
    -------
    None
        Plot.

    """
    plot_param(ax)
    fontsize = 15
    linewidth = 4
    x = elem_data(name)["Orbital"]
    y = elem_data(name)["% S Clementi"]
    ax.set_ylabel("Shielding / %", size=fontsize)
    ax.set_xlabel("Orbitals", size=fontsize)
    ax.set_xticks = [i for i in range(len(elem_data(name).index))]
    ax.set_yticks(np.arange(0, round(max(y)) + 5, 10.0))
    ax.plot(
        x,
        y,
        linewidth=linewidth,
        label=f"Shielding % Clementi {element(name).symbol}",
        **kwargs,
    )
    ax.legend(fontsize=fontsize - 2, loc="best", shadow=True, fancybox=True)


def plot_slater_both(name, ax=None):
    """Plots Effective nuclear charge and shielding /% (Slater) per orbital.

    Parameters
    ----------
    name : string
        String with the element symbol.
    ax : matplotlib.axes
        Axes for the plot (the default is None).

    Returns
    -------
    None
        Plot.

    """
    plot_param(ax)
    fontsize = 15
    linewidth = 4
    x = elem_data(name)["Orbital"]
    y = elem_data(name)["Zef Slater"]
    ax.set_ylabel("Effective nuclear charge", size=fontsize, color="blue")
    ax.set_xlabel("Orbitals", size=fontsize)
    line1 = ax.plot(
        x,
        y,
        linewidth=linewidth,
        label=f"Zef Slater {element(name).symbol}",
        color="blue",
    )
    ax.set_xticks = [i for i in range(len(elem_data(name).index))]
    ax.set_yticks(np.linspace(0, round(ax.get_ybound()[1] + 1), 5))

    ax2 = ax.twinx()
    plot_param(ax2)
    y2 = elem_data(name)["% S Slater"]
    line2 = ax2.plot(
        x,
        y2,
        linewidth=linewidth,
        linestyle="--",
        label=f"Shielding % Slater {element(name).symbol}",
        color="red",
    )
    ax2.set_ylabel("Shielding / %", size=fontsize, color="red")
    ax2.set_yticks(np.linspace(0, round(ax2.get_ybound()[1] + 1), 5))

    lines = line1 + line2
    labels = [line.get_label() for line in lines]

    ax2.legend(
        lines,
        labels,
        fontsize=fontsize - 1,
        loc="upper center",
        shadow=True,
        fancybox=True,
    )


def plot_clementi_both(name, ax=None):
    """Plots Effective nuclear charge and shielding /% (Clementi) per orbital.

    Parameters
    ----------
    name : string
        String with the element symbol.
    ax : matplotlib.axes
        Axes for the plot (the default is None).

    Returns
    -------
    None
        Plot.

    """
    plot_param(ax)
    fontsize = 15
    linewidth = 4
    x = elem_data(name)["Orbital"]
    y = elem_data(name)["Zef Clementi"]
    ax.set_ylabel("Effective nuclear charge", size=fontsize, color="blue")
    ax.set_xlabel("Orbitals", size=fontsize)
    line1 = ax.plot(
        x,
        y,
        linewidth=linewidth,
        label=f"Zef Clementi {element(name).symbol}",
        color="blue",
    )
    ax.set_xticks = [i for i in range(len(elem_data(name).index))]
    ax.set_yticks(np.linspace(0, round(ax.get_ybound()[1] + 1), 5))

    ax2 = ax.twinx()
    plot_param(ax2)
    y2 = elem_data(name)["% S Clementi"]
    line2 = ax2.plot(
        x,
        y2,
        linewidth=linewidth,
        linestyle="--",
        label=f"Shielding % Clementi {element(name).symbol}",
        color="red",
    )
    ax2.set_ylabel("Shielding / %", size=fontsize, color="red")
    ax2.set_yticks(np.linspace(0, round(ax2.get_ybound()[1] + 1), 5))

    lines = line1 + line2
    labels = [line.get_label() for line in lines]

    ax2.legend(
        lines,
        labels,
        fontsize=fontsize - 1,
        loc="upper center",
        shadow=True,
        fancybox=True,
    )
