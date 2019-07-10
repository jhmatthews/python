==============
Central_object
==============

Central_object.radiation
========================
A booliean variable stating whether of not the central object should radiate from its
survace as a star would.

**Type:** Boolean (yes/no)

**Parent(s):**

* :ref:`System_type`: ``star``, ``binary``, ``previous``


**File:** setup_star_bh.c


Central_object.temp
===================
Multi-line description, must keep indentation.

**Type:** Double

**Unit:** Kelvin

**Values:** Greater than zero

**Parent(s):**

* :ref:`System_type`: ``star``, ``binary``, ``previous``


**File:** setup_star_bh.c


Central_object.rad_type_to_make_wind
====================================
Multi-line description, must keep indentation.

**Type:** Enumerator

**Values:**

bb
  Multi-line description, must keep indentation.

models
  Multi-line description, must keep indentation.


**Parent(s):**

* :ref:`System_type`: ``star``, ``binary``, ``previous``


**File:** setup_star_bh.c


Central_object.mass
===================
Mass of the central object

**Type:** Double

**Unit:** M☉

**Values:** Greater than 0

**File:** setup_star_bh.c


Central_object.radius
=====================
Radius of the central object in the system, e.g the white dwarf or black hole

**Type:** Double

**Unit:** cm

**Values:** Greater than 0

**File:** setup_star_bh.c


Central_object.rad_type_in_final_spectrum
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


