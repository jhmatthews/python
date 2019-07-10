==
BH
==

BH.lum
======
Multi-line description, must keep indentation.

**Type:** Double

**Unit:** ergs/s

**Values:** Condition e.g. greater than 0 or list e.g. [1, 2, 5]

**Parent(s):**

* :ref:`bh.rad_type_to_make_wind`: ``brems``, ``cloudy``, ``model``, ``power``


**File:** setup_star_bh.c


BH.rad_type_to_make_wind
========================
Multi-line description, must keep indentation.

**Type:** Enumerator

**Values:**

bb
  Multi-line description, must keep indentation.

brems
  Multi-line description, must keep indentation.

cloudy
  Multi-line description, must keep indentation.

models
  Multi-line description, must keep indentation.

power
  Multi-line description, must keep indentation.


**File:** setup_star_bh.c


BH.power_law_index
------------------
Multi-line description, must keep indentation.

**Type:** Double

**Values:** Greater than 0

**Parent(s):**

* :ref:`BH.rad_type_to_make_wind`: ``cloudy``, ``power``


**File:** setup_star_bh.c


BH.blackbody_temp
-----------------
Multi-line description, must keep indentation.

**Type:** Double

**Unit:** Kelvin

**Values:** Greater than 0

**Parent(s):**

* :ref:`BH.rad_type_to_make_wind`: bb


**File:** setup_star_bh.c


BH.power_law_cutoff
-------------------
Adds a low-frequency cutoff to the power law spectrum.

**Type:** Double

**Unit:** ??? Hz ???

**Values:** Greater than or equal to 0

**Parent(s):**

* :ref:`BH.rad_type_to_make_wind`: power


**File:** setup_star_bh.c


BH.rad_type_in_final_spectrum
=============================
Multi-line description, must keep indentation.

**Type:** Enumerator

**Values:**

brems
  Multi-line description, must keep indentation.

cloudy
  Multi-line description, must keep indentation.

power
  Multi-line description, must keep indentation.


**Parent(s):**

* :ref:`Spectrum_cycles`: Greater than or equal to 0


**File:** python.c


BH.radiation
============
Whether or not the BH/AGN should radiate.

**Type:** Boolean (yes/no)

**Parent(s):**

* :ref:`System_type`: agn


**File:** setup_star_bh.c


BH.geometry_for_pl_source
=========================
Multi-line description, must keep indentation.

**Type:** Enumerator

**Values:**

lamp_post
  Multi-line description, must keep indentation.

sphere
  Multi-line description, must keep indentation.


**Parent(s):**

* :ref:`bh.radiation`: ``True``


**File:** setup_star_bh.c


