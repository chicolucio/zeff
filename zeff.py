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
    Pandas frame.
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
    ax.axhline(color='gray', zorder=-1)
    ax.axvline(color='gray', zorder=-1)
    ax.grid(b=True, which='major', linestyle=':', linewidth=2)
    ax.minorticks_on()
    ax.grid(b=True, which='minor', axis='y', linestyle=':', linewidth=1.0)
    ax.tick_params(which='both', labelsize=14)


def plot_slater_zef(name, ax=None, **kwargs):
    plot_param(ax)
    x = elem_data(name)['Orbital']
    y = elem_data(name)['Zef Slater']
    ax.set_ylabel('Carga nuclear efetiva', size=15)
    ax.set_xlabel('Orbitais', size=15)
    ax.set_xticks = ([i for i in range(len(elem_data(name).index))])
    ax.set_yticks(np.arange(0, round(max(y)) + 5, 5.0))
    ax.plot(x, y, linewidth=4, label='Zef Slater {0}'.format(
        element(name).symbol), **kwargs)
    ax.legend(fontsize=14, loc='upper right', shadow=True, fancybox=True)


def plot_slater_screening(name, ax=None, **kwargs):
    plot_param(ax)
    x = elem_data(name)['Orbital']
    y = elem_data(name)['% S Slater']
    ax.set_ylabel('Blindagem / %', size=15)
    ax.set_xlabel('Orbitais', size=15)
    ax.set_xticks = ([i for i in range(len(elem_data(name).index))])
    ax.set_yticks(np.arange(0, round(max(y)) + 5, 5.0))
    ax.plot(x, y, linewidth=4, label='Blindagem % Slater {0}'.format(
        element(name).symbol), **kwargs)
    ax.legend(fontsize=13, loc='best', shadow=True, fancybox=True)


def plot_clementi_zef(name, ax=None, **kwargs):
    plot_param(ax)
    x = elem_data(name)['Orbital']
    y = elem_data(name)['Zef Clementi']
    ax.set_ylabel('Carga nuclear efetiva', size=15)
    ax.set_xlabel('Orbitais', size=15)
    ax.set_xticks = ([i for i in range(len(elem_data(name).index))])
    ax.set_yticks(np.arange(0, round(max(y)) + 5, 5.0))
    ax.plot(x, y, linewidth=4, label='Zef Clementi {0}'.format(
        element(name).symbol), **kwargs)
    ax.legend(fontsize=14, loc='upper right', shadow=True, fancybox=True)


def plot_clementi_screening(name, ax=None, **kwargs):
    plot_param(ax)
    x = elem_data(name)['Orbital']
    y = elem_data(name)['% S Clementi']
    ax.set_ylabel('Blindagem / %', size=15)
    ax.set_xlabel('Orbitais', size=15)
    ax.set_xticks = ([i for i in range(len(elem_data(name).index))])
    ax.set_yticks(np.arange(0, round(max(y)) + 5, 5.0))
    ax.plot(x, y, linewidth=4, label='Blindagem % Clementi {0}'.format(
        element(name).symbol), **kwargs)
    ax.legend(fontsize=14, loc='lower right', shadow=True, fancybox=True)


'''
def plot_slater_both(name, ax=None):

fazer função para plotar dois gráficos, um eixo Zef e o outro blindagem.
Colorir eixos diferente

def plot_clementi_both(name, ax=None):

fazer função para plotar dois gráficos, um eixo Zef e o outro blindagem

def plot_both

fazer função para traçar 4 gráficos: dois para Zef (Slater e Clementi, linhas
cheias, azuis, eixo esquerda)
e dois para blindagem (idem, linhas vermelhas, tracejadas, eixo direita).
'''
