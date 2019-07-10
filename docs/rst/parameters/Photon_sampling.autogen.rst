===============
Photon_sampling
===============

Photon_sampling.approach
========================
Choice of whether and how to use stratified sampling in creating photons during the
ionization stage of the calculation.

**Type:** Enumerator

**Values:**

T_star
  Sets a single band based on the temperature given

cv
  Traditional cv setup

yso
  YSO setup

AGN
  Test for balance matching the bands we have been using for AGN runs

min_max_freq
  Mode 1 sets a single wide band defined by f1 and f2

user_bands
  User-defined bands

cloudy_test
  Set up to compare with cloudy power law table command note
  that this also sets up the weight and photon index for the PL, to ensure a continuous distribution

wide
  Test for balance to have a really wide frequency range

logarithmic
  Generalized method to set up logarithmic bands


**File:** bands.c


Photon_sampling.nbands
----------------------
Python uses stratified samplign to generate photons during the ionization phase.  This
parameter allows the user to define the number of bands for stratified sampling, if s/he
wants to customize the bands used for the generation of photons

**Type:** Integer

**Values:** Greater than 0

**Parent(s):**

* :ref:`Photon_sampling.approach`: user_bands


**File:** bands.c


Photon_sampling.band_min_frac
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
When specifying manually the bands used for generating photons during the ionization phase, this
parameter specifies the The minimum fraction of photons to be generated in this energy band.
The number of times this parameter will be reqested depends upon the number of bands.  The summ
of the fractions need not sum to 1, in which case the remaining photons will be distributed according
to the luminosity in the energy bands

**Type:** Double

**Values:** Greater than 0 and should sum to less than 1 over all bands

**Parent(s):**

* :ref:`Photon_sampling.nbands`: Greater than 0, once per band


**File:** bands.c


Photon_sampling.band_boundary
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
When the user specifies what bands are used for stratfied sampling, this parameter specifies the boundaries
between energy bands in which a minimum fraction of photons will be generated.  The number of times this
parameter is request depends upon the number of energies bands being used.

**Type:** Double

**Unit:** eV

**Values:** Greater than 0, monotonically increasing

**Parent(s):**

* :ref:`Photon_sampling.nbands`: Greater than 0, once per band


**File:** bands.c


Photon_sampling.high_energy_limit
---------------------------------
Stratified sampling is used during ionization cycles to generate photons.  This parameter
specifies the high energy limit for the frequencies of photons to be generated.

**Type:** Double

**Unit:** eV

**Values:** Greater than :ref:`Photon_sampling.low_energy_limit`

**Parent(s):**

* :ref:`Photon_sampling.approach`: user_bands


**File:** bands.c


Photon_sampling.low_energy_limit
--------------------------------
During the ionization phase, stratified sampling is used to provide good coverage of the full ionizing spectrum. This
parameter sets the lowest envergy (frequency) of for phtoons to be generated whne the user wants to customize the
bands.

**Type:** Double

**Unit:** eV

**Values:** Greater than 0

**Parent(s):**

* :ref:`Photon_sampling.approach`: user_bands


**File:** bands.c


