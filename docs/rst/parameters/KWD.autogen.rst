===
KWD
===

KWD.mdot_r_exponent
===================
The exponent for the mass loss rate as defined in the KWD model,
m_dot(r) = F(r) ** alpha = T(r) ** (4 * alpha).
F is the local luminous flux and T is the local temperature at a radius R. A
value of 0 sets a uniform mass loss rate.

**Type:** Double

**Values:** Greater than or equal to 0

**Parent(s):**

* :ref:`Wind.type`: kwd


**File:** knigge.c


KWD.v_zero
==========
Multiple of the local sound speed at the base of the wind, this results in
the initial velocity of the wind being able to be greater or less than the
local sound speed.

**Type:** Double

**Unit:** Sound speed at wind base

**Values:** Greater than 0

**Parent(s):**

* :ref:`Wind.type`: kwd


**File:** knigge.c


KWD.acceleration_length
=======================
The size of the acceleration length scale for a disk wind described by the
KWD model.

**Type:** Double

**Unit:** cm

**Values:** Greater than 0

**Parent(s):**

* :ref:`Wind.type`: kwd


**File:** knigge.c


KWD.rmax
========
The radius at which the disk wind terminates, in units of central object
radii. This has to be greater than rmin.

**Type:** Double

**Unit:** :ref:`Central_object.radius`

**Values:** Greater than :ref:`KWD.rmin`

**Parent(s):**

* :ref:`Wind.type`: kwd


**File:** knigge.c


KWD.acceleration_exponent
=========================
Sets the length scale over which the accleration to v_inf is accomplished.
It is the value of the exponent beta for the Caster & Lamers equation of a
stellar wind,
v(r) = v_0 + (v_inf - v_0) * (1 - R_s/r) ** beta.

**Type:** Double

**Values:** Greater than 0

**Parent(s):**

* :ref:`Wind.type`: kwd


**File:** knigge.c


KWD.rmin
========
The radius at which the disk wind begins, in units of central object radii.

**Type:** Double

**Unit:** :ref:`Central_object.radius`

**Values:** Greater than 1

**Parent(s):**

* :ref:`Wind.type`: kwd


**File:** knigge.c


KWD.d
=====
The ratio d/d_min is used to describe the degree of geometric collimation of
the disk wind in the KWD model. However, d (the distance to the focal point in
central object radii) is used as this provides a more natural parameter.

**Type:** Double

**Unit:** :ref:`Central_object.radius`

**Values:** Greater than 0

**Parent(s):**

* :ref:`Wind.type`: kwd


**File:** knigge.c


KWD.v_infinity
==============
The velocity at large distances of a steller wind described by the KWD model,
in units of escape velocity. Described in terms of Castor & Lamers equation,
v(r) = v_0 + (v_inf - v_0) * (1 - R_s/r) ** beta.

**Type:** Double

**Unit:** Escape velocity

**Values:** Greater than 0

**Parent(s):**

* :ref:`Wind.type`: kwd


**File:** knigge.c


