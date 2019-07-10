====
Disk
====

Disk.radmax
===========
The outer edge of the disk.  Photons inside this radius are
absorbed or re-radiated.  Photons which are outside this radius
pass through the disk plane.

**Type:** Double

**Unit:** cm

**Values:** Greater than 0

**Parent(s):**

* :ref:`disk.type`: ``flat``, ``vertically.extended``


**File:** setup_disk.c


Disk.type
=========
Parameter defining whether there is a disk in the system

**Type:** Enumerator

**Values:**

none
  No disk

flat
  Standard flat disk

vertically.extended
  Vertically extended disk


**File:** setup_disk.c


Disk.rad_type_to_make_wind
--------------------------
Multi-line description, must keep indentation.

**Type:** Enumerator

**Values:**

bb
  Multi-line description, must keep indentation.

models
  Multi-line description, must keep indentation.


**Parent(s):**

* :ref:`Disk.type`: ``flat``, ``vertically.extended``


**File:** setup_disk.c


Disk.radiation
--------------
Multi-line description, must keep indentation.

**Type:** Boolean(yes/no)

**Parent(s):**

* :ref:`Disk.type`: ``flat``, ``vertically.extended``


**File:** setup_disk.c


Disk.temperature.profile
^^^^^^^^^^^^^^^^^^^^^^^^
The choice of disk temperature profile

**Type:** Enumerator

**Values:**

standard
  A Shakura - Sunyaev  disk, with a hard inner boundar

readin
  Read the profile in from a file; the user will be queried for the name of the file

yso
  YSO???

analytic
  DEPRECATED??? A profile designed for the situation where the disk is being illuminated by star


**Parent(s):**

* :ref:`Disk.radiation`: ``True``


**File:** setup_disk.c


**Disk.mdot**
"""""""""""""
The mass transfer rate in the disk when considering a standard Shakura-disk.

**Type:** Double

**Unit:** M☉/year

**Parent(s):**

* :ref:`Disk.temperature.profile`: standard


**File:** setup_disk.c


**Disk.T_profile_file**
"""""""""""""""""""""""
When the user chooses to read in the temperature profile as a
function of radius, the user is asked the name of the file that
contains the desired profile.

**Type:** String

**Parent(s):**

* :ref:`Disk.temperature.profile`: readin


**File:** setup_disk.c


Disk.z1
-------
For a vertically extended the disk, the height of the disk is
set to be :ref:`Disk.z0` * :ref:`Disk.radmax` * (r/:ref:`Disk.radmax`)**Disk.z1 where Disk.z1
is the power law index

**Type:** Double

**Values:** Greater than 0

**Parent(s):**

* :ref:`Disk.type`: vertically.extended


**File:** setup_disk.c


Disk.z0
-------
Fractional height at maximum radius.  The physical height at the
outer disk will be this * :ref:`Disk.radmax`.

**Type:** Double

**Values:** Greater than 0

**Parent(s):**

* :ref:`Disk.type`: vertically.extended


**File:** setup_disk.c


Disk.rad_type_in_final_spectrum
===============================
Multi-line description, must keep indentation.

**Type:** Enumerator

**Values:**

bb
  Multi-line description, must keep indentation.

models
  Multi-line description, must keep indentation.

uniform
  Multi-line description, must keep indentation.


**Parent(s):**

* :ref:`Spectrum_cycles`: Greater than or equal to 0


**File:** python.c


