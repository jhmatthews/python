=====
Hydro
=====

Hydro.thetamax
==============
The maximum theta value to be read in from a hydrodynamic snapshot.
This is typically used to excise a dense disk from the midplane of
such a snapshot. Use a negative value to tell the code to use all
the data.

**Type:** Double

**Unit:** Degrees

**Values:**

-1
  use all data

X
  use up to that angle (typically less than 90)


**Parent(s):**

* :ref:`Wind.type`: hydro


**File:** hydro_import.c


Hydro.file
==========
Relative path to a hydrodynamic snapshot file to be imported.

**Type:** String

**Parent(s):**

* :ref:`Wind.type`: hydro


**File:** hydro_import.c


