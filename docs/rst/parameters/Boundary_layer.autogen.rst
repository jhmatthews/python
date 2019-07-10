==============
Boundary_layer
==============

Boundary_layer.radiation
========================
Says whether the boundary layer will radiate.

**Type:** Boolean (yes/no)

**Parent(s):**

* :ref:`System_type`: ``star``, ``binary``


**File:** setup_star_bh.c


Boundary_layer.luminosity
-------------------------
The luminosity of the boundary layer

**Type:** Double

**Unit:** ergs/s

**Values:** Greater than 0

**Parent(s):**

* :ref:`Boundary_layer.radiation`: ``True``

* :ref:`Boundary_layer.rad_type_to_make_wind`: power


**File:** setup_star_bh.c


Boundary_layer.temp
-------------------
The temperature of the boundary layer in situations where temperature
is meaningful in generating the spectrum

**Type:** Double

**Unit:** Kelvin

**Values:** Greater than 0

**Parent(s):**

* :ref:`Boundary_layer.radiation`: ``True``

* :ref:`Boundary_layer.rad_type_to_make_wind`: power


**File:** setup.c


Boundary_layer.rad_type_to_make_wind
====================================
Multi-line description, must keep indentation.

**Type:** Enumerator

**Values:**

bb
  Multi-line description, must keep indentation.

models
  Multi-line description, must keep indentation.

power
  Multi-line description, must keep indentation.


**File:** setup_star_bh.c


Boundary_layer.rad_type_in_final_spectrum
=========================================
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


