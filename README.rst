===============================
CINE: Comet INfrared Excitation
===============================


cine is a tool for calculating infrared pumping efficiencies.  At large
nucleocentric distances, one of the main mechanisms for molecular excitation in
comets is the fluorescence by the solar radiation followed by radiative decay
to the ground vibrational state.  This code calculates the effective pumping
rates for rotational levels in the ground vibrational state scaled by the
heliocentric distance of the comet.  These coefficients are useful for modeling
rotational emission lines observed in cometary spectra at sub-millimeter
wavelengths.

Code releases are available on `PyPI <https://pypi.python.org/pypi/cine>`_, and
development happens in the `github project page
<https://github.com/migueldvb/cine>`_.


Installation
------------

cine can be installed using `pip <https://pypi.python.org/pypi/pip>`_:

.. code-block:: bash

    $ pip install cine

or by cloning the github repository:

.. code-block:: bash

    $ # If you have a github account:
    $ git clone git@github.com:migueldvb/cine.git
    $ # If you do not:
    $ git clone https://github.com/migueldvb/cine.git
    $ cd cine
    $ python setup.py install
    $ # Or if you do not have root privileges:
    $ python setup.py install --user


Requirements
------------

The code requires the standard scientific Python packages (`numpy
<http://www.numpy.org/>`_, `scipy <https://www.scipy.org/>`_, and `pandas
<http://pandas.pydata.org/>`_) and astropy's affiliated package `astroquery
<https://github.com/astropy/astroquery>`_.  to access the HITRAN and Lamda
databases. Running the tests requires `nose
<https://pypi.python.org/pypi/nose>`_.


Example
-------

``cine`` is a command-line tool that is included in the package to generate
pumping rates for several molecules. For example, to obtain the effective
pumping rates between the seven lowest rotational levels in the ground
vibrational state of HDO you can run the following command once ``CINE`` has
been installed:

.. code-block:: bash

    $ cine --mol HDO --nlevels 7

This should create a file named ``G_HDO.dat`` which contains the pumping rates
G :subscript:`ij` in units of s :superscript:`-1` between the rotational levels
i and j shown in the first two columns. Note that the levels use zero-based
indexing.

.. code-block:: bash

    0 3 2.568872e-05
    0 4 2.570305e-05
    0 5 1.552757e-05
    1 2 6.253229e-05
    1 6 2.987896e-05
    2 1 6.196215e-05
    2 6 4.410062e-05
    3 0 7.547422e-05
    3 4 3.103947e-05
    3 5 5.048423e-05
    4 0 1.253741e-04
    4 3 5.128064e-05
    4 5 4.679292e-05
    5 0 7.481781e-05
    5 3 8.287649e-05
    5 4 4.643613e-05
    6 1 4.820172e-05
    6 2 7.201329e-05

To include more levels in the calculation, change the ``-n/-nlevels`` command-line
option to a larger value.  cine has a ``-h/--help`` argument that presents an
usage explanation describing each optional argument.

These coefficients are useful for deriving molecular production rates from cometary
lines observed at sub-millimeter wavelengths combined with a code that
solves the radiative transfer equations such as `LIME
<https://github.com/lime-rt/lime>`_.


Downloading HITRAN data
-----------------------

To download the molecular data cine uses the ``astroquery.hitran`` and
``astroquery.lamda`` tools.  Set the ``LAMDA_DATA`` and ``HITRAN_DATA``
environment variables (otherwise, the default
``~/.astropy/cache/astroquery/Lamda``  and
``~/.astropy/cache/astroquery/hitran`` will be used),


Tests
-----

If ``nose`` is installed the tests can be run from the root of the repository as:

.. code-block:: bash

    $ python setup.py test


Contributing
------------

Any questions or bug reports can be raised in github's `issue tracker
<https://github.com/migueldvb/cine/issues>`_ or `pull requests
<https://github.com/migueldvb/cine/pulls>`_.


License
-------

Copyright 2017 Miguel de Val-Borro

``CINE`` is free software made available under the MIT License.
For details see the LICENSE file.
