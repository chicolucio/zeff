.. These are examples of badges you might want to add to your README:
   please update the URLs accordingly

    .. image:: https://readthedocs.org/projects/zeff/badge/?version=latest
        :alt: ReadTheDocs
        :target: https://zeff.readthedocs.io/en/stable/
    .. image:: https://img.shields.io/pypi/v/zeff.svg
        :alt: PyPI-Server
        :target: https://pypi.org/project/zeff/

.. image:: https://img.shields.io/badge/Author-Francisco%20Bustamante-red.svg
    :alt: Francisco Bustamante
    :target: https://www.linkedin.com/in/flsbustamante
.. image:: https://img.shields.io/badge/Python-3.10+-blue.svg
    :alt: Python
    :target: https://www.python.org/
.. image:: https://img.shields.io/badge/License-MIT-blue.svg
    :alt: LICENSE
    :target: LICENSE.txt
.. image:: https://img.shields.io/badge/Contributions-Welcome-brightgreen.svg?style=flat
    :alt: Contributions are welcome
    :target: https://github.com/chicolucio/zeff/issues
.. image:: https://img.shields.io/badge/-PyScaffold-005CA0?logo=pyscaffold
    :alt: Project generated with PyScaffold
    :target: https://pyscaffold.org/

===================
Z\ :subscript:`eff`
===================

    An effective nuclear charge calculator.

.. image:: https://raw.githubusercontent.com/chicolucio/zeff/master/images/zeff_plot.png
   :alt: Image Zeff

Introduction
============

If you are a chemistry teacher or student, chances are that you are going to do Z\
:sub:`eff` and S calculations based on `Slater's Rules`_ or on `Clementi *et al* screening
constants`_. Probably while teaching/learning periodic trends.

This simple project automates these calculations returning them as tables (Pandas
DataFrames) so that one can easily plot the data or manipulate them. Indeed, there are
plot functions that can do graphs providing a convenient visualization functionality for
classes.

Installation
============

Just clone or download this repo. This is not a package (yet, maybe someday :-))

Usage
=====

Look at some examples at the `tutorial`_ available as a `Jupyter notebook`_.

Under the hood - Requirements
=============================

This project relies heavily on `mendeleev`_, `pandas`_ and `matplotlib`_ packages so these must be installed. Check the `requirements.txt`_ file
for requirements and feel free to use it to create a virtual environment.

Contributing
============

All contributions are welcome.

**Issues**

Feel free to submit issues regarding:

- recommendations
- more examples for the tutorial
- enhancement requests and new useful features
- code bugs

**Pull requests**

- before starting to work on your pull request, please submit an issue first
- fork the repo
- clone the project to your own machine
- commit changes to your own branch
- push your work back up to your fork
- submit a pull request so that your changes can be reviewed

For full details on how to contribute, check out the `Contributing Guide`_.

License
=======

MIT, see `LICENSE`_.

Citing
======

If you use Z\ :sub:`eff` in a scientific publication or in classes, please consider
citing as:

F. L. S. Bustamante, *Zeff* - An effective nuclear charge (Z\ :sub:`eff`) and
shielding (S) calculator and graphing tool., 2019 - Available at: `GitHub repository`_


.. _Slater's Rules: https://en.wikipedia.org/wiki/Slater%27s_rules
.. _Clementi *et al* screening constants: https://en.wikipedia.org/wiki/Effective_nuclear_charge#Values
.. _tutorial: notebooks/Zeff_tutorial.ipynb
.. _Jupyter notebook: https://jupyter.org/
.. _mendeleev: https://pypi.org/project/mendeleev/
.. _pandas: https://pandas.pydata.org/
.. _matplotlib: https://matplotlib.org/
.. _requirements.txt: requirements.txt
.. _Contributing Guide: CONTRIBUTING.rst
.. _LICENSE: LICENSE.txt
.. _GitHub repository: https://github.com/chicolucio/zeff
