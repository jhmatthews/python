============
Stellar_wind
============

Stellar_wind.mdot
=================
Mass loss rate for a wind modelled in terms of the
Caster and Lamemers prescription for a stellar wind.

**Type:** Double

**Unit:** M☉/year

**Values:** Greater than 0

**Parent(s):**

* :ref:`Wind.type`: star


**File:** stellar_wind.c


Stellar_wind.radmin
===================
Inner edge in cm for a stellar wind, normally the
radius of the star.

**Type:** Double

**Unit:** cm

**Values:** Greater than or equal to :ref:`Central_object.radius`

**Parent(s):**

* :ref:`Wind.type`: star


**File:** stellar_wind.c


Stellar_wind.radmax
===================
Multi-line description, must keep indentation.

**Type:** Double

**Unit:** cm

**Values:** Greater than or equal to :ref:`Stellar_wind.radmin`

**Parent(s):**

* :ref:`Wind.type`: star


**File:** stellar_wind.c


Stellar_wind.acceleration_exponent
==================================
Exponent beta for the Caster and Lamers description of a stellar wind
v(r)=v_o + (v_inf - v_o) (1+R_s/r)**beta

**Type:** Double

**Values:** Greater than or equal to 0

**Parent(s):**

* :ref:`Wind.type`: star


**File:** stellar_wind.c


Stellar_wind.v_infinity
=======================
The velocity at large distance of a stellar wind described in terms
of the Casters and Larmers equation
v(r)=v_o + (v_inf - v_o) (1+R_s/r)**beta

**Type:** Double

**Unit:** cm/s

**Values:** Greater than 0

**Parent(s):**

* :ref:`Wind.type`: star


**File:** stellar_wind.c


Stellar_wind.vbase
==================
Multi-line description, must keep indentation.

**Type:** Double

**Unit:** cm/s

**Values:** Condition e.g. greater than 0 or list e.g. [1, 2, 5]

**Parent(s):**

* :ref:`Wind.type`: star


**File:** stellar_wind.c


