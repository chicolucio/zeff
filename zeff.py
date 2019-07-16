from mendeleev import element
import numpy as np
import pandas as pd


def orbitals(name):
    """Extracts n, l, and nl notation from electronic configuration info of
    Mendeleev package.

    Parameters
    ----------
    name : string
        Name or symbol for the chemical element.

    Returns
    -------
    lists
        Lists with orbital representation; quantrum principal number; orbital
        letter notation and quantum azimutal number.
    """
    elem = element(name)
    orbital_n = []
    orbital_l = []

    for k, v in elem.ec.conf.items():
        orbital_n.append(k[0])
        orbital_l.append(k[1])

    dict_l = {'s': 0, 'p': 1, 'd': 2, 'f': 3}
    orbital_l_num = []
    for i in orbital_l:
        if i in dict_l:
            orbital_l_num.append(dict_l[i])

    orbital = []
    for i in range(len(orbital_n)):
        orbital.append(str(orbital_n[i]) + orbital_l[i])

    return orbital, orbital_n, orbital_l, orbital_l_num


def slater(name):
    """In quantum chemistry, Slater's rules provide numerical values for the
    effective nuclear charge concept. In a many-electron atom, each electron is
    said to experience less than the actual nuclear charge owing to shielding
    or screening by the other electrons.

    Parameters
    ----------
    name : string
        Name or symbol for the chemical element.

    Returns
    -------
    type list
        Three lists. The first one with the screening values, the second with
        the effective nuclear charge values and the third one with the
        screening percentage for each orbital of the electronic configuration
        of the element.

    """
    elem = element(name)
    slater_s = []
    zeff_slater = []

    for k, v in elem.ec.conf.items():
        slater_s.append(elem.ec.slater_screening(k[0], k[1]))
        zeff_slater.append(elem.atomic_number -
                           elem.ec.slater_screening(k[0], k[1]))

    slater_s_percent = [(i / elem.atomic_number) * 100 for i in slater_s]

    return slater_s, zeff_slater, slater_s_percent


def clementi(name):
    """In a many-electron atom, each electron is said to experience less than
    the actual nuclear charge owing to shielding or screening by the other
    electrons. This functions calculates screening and effective nuclear charge
    based on screening values calculate by Clementi and Raimond with SCF
    algorithms. No values for Z > 86.

    Parameters
    ----------
    name : string
        Name or symbol for the chemical element.

    Returns
    -------
    lists
        Three lists. The first one with the screening values, the second with
        the effective nuclear charge values and the third one with the
        screening percentage for each orbital of the electronic configuration
        of the element.
    """
    elem = element(name)
    clementi_s = []
    zeff_clementi = []

    if elem.atomic_number < 87:
        for k, v in elem.ec.conf.items():
            clementi_s.append(elem.atomic_number -
                              elem.zeff(k[0], k[1], method='clementi'))
            zeff_clementi.append(elem.zeff(k[0], k[1], method='clementi'))
    else:
        for i in range(len(orbitals.orbital_n)):
            clementi_s.append(np.nan)
            zeff_clementi.append(np.nan)

    clementi_s_percent = [(i / elem.atomic_number) * 100 for i in clementi_s]

    return clementi_s, zeff_clementi, clementi_s_percent


def elem_data(name):
    """Data summary for a given element.

    Parameters
    ----------
    name : string
        Name or symbol for the chemical element.

    Returns
    -------
    Pandas DataFrame.
        Summary of Zeff and S for a given element with Slater and Clementi
        screening values.
    """

    slater_s, zeff_slater, slater_s_percent = slater(name)
    clementi_s, zeff_clementi, clementi_s_percent = clementi(name)
    orbital, orbital_n, orbital_l, orbital_l_num = orbitals(name)

    data = pd.DataFrame(data={'n': orbital_n,
                              'l': orbital_l,
                              'l_num': orbital_l_num,
                              'Orbital': orbital,
                              'S Slater': slater_s,
                              '% S Slater': slater_s_percent,
                              'Zef Slater': zeff_slater,
                              'S Clementi': clementi_s,
                              '% S Clementi': clementi_s_percent,
                              'Zef Clementi': zeff_clementi})

    data_org = data.sort_values(by=['n', 'l_num']).reset_index(drop=True)

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
    ax.axhline(color='gray', zorder=-1)
    ax.axvline(color='gray', zorder=-1)
    ax.grid(b=True, which='major', linestyle=':', linewidth=2)
    ax.minorticks_on()
    ax.grid(b=True, which='minor', axis='y', linestyle=':', linewidth=1.0)
    ax.tick_params(which='both', labelsize=14)


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
    x = elem_data(name)['Orbital']
    y = elem_data(name)['Zef Slater']
    ax.set_ylabel('Effective nuclear charge', size=fontsize)
    ax.set_xlabel('Orbitals', size=fontsize)
    ax.set_xticks = ([i for i in range(len(elem_data(name).index))])
    ax.set_yticks(np.arange(0, round(max(y)) + 5, 5.0))
    ax.plot(x, y, linewidth=linewidth, label='Zef Slater {0}'.format(
        element(name).symbol), **kwargs)
    ax.legend(fontsize=fontsize - 1, loc='best', shadow=False, fancybox=True,
              bbox_to_anchor=(1, 1))


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
    x = elem_data(name)['Orbital']
    y = elem_data(name)['% S Slater']
    ax.set_ylabel('Shielding / %', size=fontsize)
    ax.set_xlabel('Orbitals', size=fontsize)
    ax.set_xticks = ([i for i in range(len(elem_data(name).index))])
    ax.set_yticks(np.arange(0, round(max(y)) + 5, 10.0))
    ax.plot(x, y, linewidth=linewidth, label='Shielding % Slater {0}'.format(
        element(name).symbol), **kwargs)
    ax.legend(fontsize=fontsize - 2, loc='best', shadow=True, fancybox=True)


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
    x = elem_data(name)['Orbital']
    y = elem_data(name)['Zef Clementi']
    ax.set_ylabel('Effective nuclear charge', size=fontsize)
    ax.set_xlabel('Orbitals', size=fontsize)
    ax.set_xticks = ([i for i in range(len(elem_data(name).index))])
    ax.set_yticks(np.arange(0, round(max(y)) + 5, 5.0))
    ax.plot(x, y, linewidth=linewidth, label='Zef Clementi {0}'.format(
        element(name).symbol), **kwargs)
    ax.legend(fontsize=fontsize - 1, loc='upper right',
              shadow=True, fancybox=True)


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
    x = elem_data(name)['Orbital']
    y = elem_data(name)['% S Clementi']
    ax.set_ylabel('Shielding / %', size=fontsize)
    ax.set_xlabel('Orbitals', size=fontsize)
    ax.set_xticks = ([i for i in range(len(elem_data(name).index))])
    ax.set_yticks(np.arange(0, round(max(y)) + 5, 10.0))
    ax.plot(x, y, linewidth=linewidth, label='Shielding % Clementi {0}'.format(
        element(name).symbol), **kwargs)
    ax.legend(fontsize=fontsize - 1, loc='lower right',
              shadow=True, fancybox=True)


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
    x = elem_data(name)['Orbital']
    y = elem_data(name)['Zef Slater']
    ax.set_ylabel('Effective nuclear charge', size=fontsize, color='blue')
    ax.set_xlabel('Orbitals', size=fontsize)
    line1 = ax.plot(x, y, linewidth=linewidth, label='Zef Slater {0}'.format(
        element(name).symbol), color='blue')
    ax.set_xticks = ([i for i in range(len(elem_data(name).index))])
    ax.set_yticks(np.linspace(0, round(ax.get_ybound()[1] + 1), 5))

    ax2 = ax.twinx()
    plot_param(ax2)
    y2 = elem_data(name)['% S Slater']
    line2 = ax2.plot(x, y2, linewidth=linewidth, linestyle='--',
                     label='Shielding % Slater {0}'.format(element(name).symbol), color='red')
    ax2.set_ylabel('Shielding / %', size=fontsize, color='red')
    ax2.set_yticks(np.linspace(0, round(ax2.get_ybound()[1] + 1), 5))

    lines = line1 + line2
    labels = [l.get_label() for l in lines]

    ax2.legend(lines, labels, fontsize=fontsize - 1,
               loc='upper center', shadow=True, fancybox=True)


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
    x = elem_data(name)['Orbital']
    y = elem_data(name)['Zef Clementi']
    ax.set_ylabel('Effective nuclear charge', size=fontsize, color='blue')
    ax.set_xlabel('Orbitals', size=fontsize)
    line1 = ax.plot(x, y, linewidth=linewidth,
                    label='Zef Clementi {0}'.format(element(name).symbol),
                    color='blue')
    ax.set_xticks = ([i for i in range(len(elem_data(name).index))])
    ax.set_yticks(np.linspace(0, round(ax.get_ybound()[1] + 1), 5))

    ax2 = ax.twinx()
    plot_param(ax2)
    y2 = elem_data(name)['% S Clementi']
    line2 = ax2.plot(x, y2, linewidth=linewidth, linestyle='--',
                     label='Shielding % Clementi {0}'.format(element(name).symbol), color='red')
    ax2.set_ylabel('Shielding / %', size=fontsize, color='red')
    ax2.set_yticks(np.linspace(0, round(ax2.get_ybound()[1] + 1), 5))

    lines = line1 + line2
    labels = [l.get_label() for l in lines]

    ax2.legend(lines, labels, fontsize=fontsize - 1,
               loc='upper center', shadow=True, fancybox=True)
