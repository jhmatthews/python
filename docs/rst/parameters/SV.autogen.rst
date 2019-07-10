==
SV
==

SV.thetamax
===========
The angle at which the wind rises from the outermost launching radius in a Shlossman-Vitello type disk wind.
This angle is measured with respect to the vertical (z) direction i.e. zero describes a vertical wind.
See figure 1 of Shlossman & Vitello 1993, ApJ 409,372.

**Type:** Double

**Unit:** Degrees

**Values:** Greater than sv.thetamin

**Parent(s):**

* :ref:`Wind.type`: SV


**File:** sv.c


SV.acceleration_exponent
========================
Power-law acceleration exponent (i.e. alpha) of a line driven wind in a Shlosman & Vitello (SV) CV disk wind model.
Sets the length scale over which the accleration to v_inf is accomplished.
This value is a constant; when equal to 1 the results resemble those of a linear velocity law.
Typically for an SV type wind this power law exponent is 1.5.
See equation (2) Shlosman & Vitello 1993, ApJ 409, 372.

**Type:** Double

**Values:** Greater than 0

**Parent(s):**

* :ref:`Wind.type`: SV


**File:** sv.c


SV.diskmin
==========
The innermost radius from which the wind rises in a Shlossman-Vitello type disk wind.
This radius is measured along the radial disk (r) direction i.e. zero describes the centre of the central object
(white dwarf)
See figure 1 of Shlosman & Vitello 1993, ApJ 409,372.

**Type:** Double

**Unit:** cm

**Values:** Greater than or equal to :ref:`Central_object.radius`

**Parent(s):**

* :ref:`Wind.type`: SV


**File:** sv.c


SV.diskmax
==========
The outermost radius from which the wind rises in a Shlossman-Vitello type disk wind.
This radius is measured along the radial disk (r) direction i.e. zero describes the centre of the central object
(white dwarf)
See figure 1 of Shlosman & Vitello 1993, ApJ 409,372.

**Type:** Double

**Unit:** cm

**Values:** Greater than or equal to :ref:`sv.diskmin` (inner radius disk wind)

**Parent(s):**

* :ref:`Wind.type`: SV


**File:** sv.c


SV.acceleration_length
======================
The size of the acceleration length scale for a disk wind described by the
Shlosman Vitelo model. See equation (2) Shlosman & Vitelo ApJ (1993),409,372

**Type:** Double

**Unit:** cm

**Values:** Greater than 0

**Parent(s):**

* :ref:`Wind.type`: SV


**File:** sv.c


SV.mdot_r_exponent
==================
The exponent for the mass loss rate as defined in the Shlosman Vitelo model,
See lambda in equation (4) Shlosman & Vitelo,ApJ,1993,409,372.

**Type:** Double

**Values:** Greater than or equal to 0. 0 sets a uniform mass loss rate.

**Parent(s):**

* :ref:`Wind.type`: SV


**File:** sv.c


SV.v_zero_mode
==============
Multi-line description, must keep indentation.

**Type:** Enumerator

**Values:**

fixed
  Multi-line description, must keep indentation.

sound_speed
  Multi-line description, must keep indentation.


**Parent(s):**

* :ref:`Wind.type`: SV


**File:** sv.c


SV.v_zero
---------
The velocity at the wind base.

**Type:** Double

**Unit:** ['Speed of sound in the wind', 'cm/s']

**Values:** Greater than 0

**Parent(s):**

* :ref:`SV.v_zero_mode`: ``sound_speed``, ``fixed``


**File:** sv.c


SV.v_infinity
=============
Asymptotic (i.e. final) velocity of a line driven wind in a Shlosman & Vitello CV disk wind model.
Assumed to scale with the local velocity at the base of the streamline.
See equation (2) Shlosman & Vitello 1993, ApJ 409, 372.

**Type:** Double

**Unit:** Escape velocity

**Values:** Greater than 0

**Parent(s):**

* :ref:`Wind.type`: SV


**File:** sv.c


SV.thetamin
===========
The angle at which the wind rises from the innermost launching radius in a Shlossman-Vitello type disk wind.
This angle is measured with respect to the vertical (z) direction. I.e. zero descirbes a vertical wind.
See figure 1 of Shlossman & Vitello 1993, ApJ, 409, 372.

**Type:** Double

**Unit:** Degrees

**Values:** Greater than 0

**Parent(s):**

* :ref:`Wind.type`: SV


**File:** sv.c


