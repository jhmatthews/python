========
Spectrum
========

Spectrum.orbit_phase
====================
For binary systems, the orbital phase at which the spectrum
is to be extracted (so the effects of an eclipse can be taken
into account in creating the spectrum). Phase 0 corresponds to
inferior conjunciton, that is with the secondary in front (or
depending on inclination angle, partially in front of) the
primary

**Type:** Double

**Values:** Between 0 and 1

**Parent(s):**

* :ref:`Spectrum_cycles`: Greater than or equal to 0

* :ref:`System_type`: binary


**File:** setup.c


Spectrum.no_observers
=====================
The number of different inclination angles for which spectra
will be extracted.

**Type:** Integer

**Values:** Greater than 0

**Parent(s):**

* :ref:`Spectrum_cycles`: Greater than or equal to 0


**File:** setup.c


Spectrum.angle
--------------
The inclination angle with respect to the polar axis for
obtaining a spectrum.  This question will be repeated once
for each desired incliniation

**Type:** Double

**Unit:** Degrees

**Values:** 0 to 90 degrees, where 0 is normal to the disk and 90 is on the disk plane

**Parent(s):**

* :ref:`Spectrum.no_observers`: Greater than 0. Once per observer.


**File:** setup.c


Spectrum.wavemin
================
The minimum wavelength of the final spectra in Angstroms

**Type:** Double

**Unit:** Angstroms

**Values:** Greater than 0

**Parent(s):**

* :ref:`Spectrum_cycles`: Greater than or equal to 0


**File:** setup.c


Spectrum.select_photons_by_position
===================================
Advanced command associated with adding conditions for
the detailed spectra that are extracted.  This command simply
asks whether one would like to select photons by position.  If
so one will be asked to define a spheical region in interms of
its cylindrical coordinates.

**Type:** Boolean (yes/no)

**Parent(s):**

* :ref:`Spectrum_cycles`: Greater than or equal to 0


**File:** setup.c


Spectrum.select_location
------------------------
One of several related parameters that permit one to apply
additional conditions on the location of photons extracted in
the detailed spectrum. The location refers here to the either
where the photons was created or where it last scattered

**Type:** Enumerator

**Values:**

all
  Select photons regardless of where they are generated

below_disk
  Select only photons generated from below (-z) the disk

above_disk
  Select only photons orginating above the disk

spherical_region
  Select photons by defining a spherical region


**Parent(s):**

* :ref:`Spectrum.select_photons_by_position`: ``True``


**File:** setup.c


Spectrum.select_r
^^^^^^^^^^^^^^^^^
Part of a set of parameters which define a spherical region of space from which
photons are to be extracted. select_r defines the radius of the spherical region

**Type:** Double

**Unit:** cm

**Values:** Greater than 0

**Parent(s):**

* :ref:`Spectrum.select_location`: spherical_region


**File:** setup.c


Spectrum.select_rho
^^^^^^^^^^^^^^^^^^^
Advanced command which defines a spherical  region of
space from which photons are to be extracted in constructing a detailed
spectrum.  The region is defined by a cylindrical distance, and z height
and an aximuth, and a radius r.  This parameter defines the rho coordiante
of the region.

**Type:** Double

**Unit:** cm

**Values:** Condition e.g. greater than 0 or list e.g. [1, 2, 5]

**Parent(s):**

* :ref:`Spectrum.select_location`: spherical_region


**File:** setup.c


Spectrum.select_azimuth
^^^^^^^^^^^^^^^^^^^^^^^
Advance command which along with several other parameters
specifies a spherical region of space in cylindrical coordinates.
This parameter desribes the azimuth of the region.  When
this general option is used, a detailed spectrum is constructed
just from photons that originate or scatter int he region

**Type:** Double

**Unit:** Degrees

**Values:** Between 0, and 360 or -180 to 180

**Parent(s):**

* :ref:`Spectrum.select_location`: spherical_region


**File:** setup.c


Spectrum.select_z
^^^^^^^^^^^^^^^^^
Advanced command which defines a spherical  region of
space from which photons are to be extracted in constructing a detailed
spectrum.  The region is defined by a cylindrical distance, and z height
and an aximuth, and a radius r.  This parameter defines the z coordiante
of the region.

**Type:** Double

**Unit:** cm

**Values:** Within the z range of the model

**Parent(s):**

* :ref:`Spectrum.select_location`: spherical_region


**File:** setup.c


Spectrum.type
=============
The type of spectra that are produced in the final spectra. The current choices are flambda, fnu, or basic,
where basic implies simply summing up the energy packets that escape within a particularly wavelength/
frequency bin.

**Type:** Enumerator

**Values:**

flambda
  λF(λ)

fnu
  νF(ν)

basic
  F(λ)


**Parent(s):**

* :ref:`Spectrum_cycles`: Greater than or equal to 0


**File:** setup.c


Spectrum.live_or_die
====================
Normally in creating detailed spectrum Python "extracts" photons in a certain
direction reweighting them to account for the fact that they have been extracted
in a certain direction.  It is possible to just count the photons that are emitted
in a single angle range. The two methods should yield the same or very similar results
but the extraction method is much more efficient and live or die is basically a
diagnostic mode.

**Type:** Enumerator

**Values:**

live.or.die
  Count only those photons that escape within a small angle range towards the observer

extract
  Extract a component of all photons that scatter towards the observer


**Parent(s):**

* :ref:`Spectrum_cycles`: Greater than or equal to 0


**File:** setup.c


Spectrum.select_specific_no_of_scatters_in_spectra
==================================================
Advanced command which allows one to place additional
constraints on the detailed spectra which are extract.
This includes selectiong photons from above or below the
disk, only photons which have scttered, etc.

**Type:** Boolean (yes/no)

**Parent(s):**

* :ref:`Spectrum_cycles`: Greater than or equal to 0


**File:** setup.c


Spectrum.select_scatters
------------------------
Advaned command that allows one to extract photons that
have undergone a certain number of scatters.  If n > MAXSCAT,
that is to say a very large number then all scatters are slected.
If lies between 0 and MAXSCAT then photons will be extracted only
at the point a photon has undergone this number of scatters.  If
n is < 0 then photons with n or greater scattters will be extracted.

**Type:** Integer

**Values:** Greater than 0

**Parent(s):**

* :ref:`Spectrum.select_specific_no_of_scatters_in_spectra`: ``True``


**File:** setup.c


Spectrum.wavemax
================
The maximum wavelength of the detailed spectra that are to be produced

**Type:** Double

**Unit:** Angstroms

**Values:**

Spectrum.wavemin
  Greater than


**Parent(s):**

* :ref:`Spectrum_cycles`: Greater than or equal to 0


**File:** setup.c


