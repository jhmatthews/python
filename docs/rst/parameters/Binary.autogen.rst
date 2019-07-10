======
Binary
======

Binary.mass_sec
===============
In binary systems the mass of the secondary. This is used along
with the period to establish the Roche lobes, so that one can
see the effects of eclipses on the system

**Type:** Double

**Unit:** M☉/year

**Values:** Greater than 0

**Parent(s):**

* :ref:`System_type`: binary


**File:** setup_star_bh.c


Binary.period
=============
The perids of a binary system. Along with a mass, the binary period is
used to define the Roche lobe of the system, which in turn can be used
to see the effect of eclipses on the spectrum.  Defining the system as
a secondary also initializes the outer radius of the disk.

**Type:** Double

**Unit:** Hours

**Values:** Greater than 0

**Parent(s):**

* :ref:`System_type`: binary


**File:** setup_star_bh.c


