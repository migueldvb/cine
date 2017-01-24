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

Code releases are available on `PyPI <http://pypi.python.org/pypi/cine>`_, and
development happens in the `github project page
<http://github.com/migueldvb/cine>`_.


Installation and Dependencies
-----------------------------

The code requires the standard scientific Python packages (numpy, scipy, and
pandas) and astropy's affiliated package astroquery to access the HITRAN and
Lamda databases. Running the tests requires `nose <https://pypi.python.org/pypi/nose>`_.

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


License
-------

Copyright 2017 Miguel de Val-Borro

cine is free software made available under the MIT License.
For details see the LICENSE file.
