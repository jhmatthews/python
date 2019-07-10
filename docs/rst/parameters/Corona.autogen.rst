======
Corona
======

Corona.radmax
=============
The corona is a box-shaped region which sits immediately
above the disk.  radmax defines the outer edge of the box.

**Type:** Double

**Unit:** cm

**Values:** Greater than :ref:`Central_object.radius`

**Parent(s):**

* :ref:`Wind.type`: corona


**File:** corona.c


Corona.base_den
===============
The coronal model is defined in terms of a base density
and a scale height

**Type:** Double

**Unit:** number/cm**3

**Values:** Greater than 0

**Parent(s):**

* :ref:`Wind.type`: corona


**File:** corona.c


Corona.radmin
=============
The corona is a box-shaped region which sits immediately
above the disk. radmin defines the inner edge of the box.

**Type:** Double

**Unit:** cm

**Values:** Greater than :ref:`Central_object.radius`

**Parent(s):**

* :ref:`Wind.type`: corona


**File:** corona.c


Corona.zmax
===========
The corona is a box-shaped region which sits immediately
above the disk.  zmax defines the height of the box.

**Type:** Double

**Unit:** cm

**Values:** Greater than that the radius of the central object

**Parent(s):**

* :ref:`Wind.type`: corona


**File:** corona.c


Corona.scale_height
===================
The coronal model is defined in terms of a base density
and a scale height

**Type:** Double

**Unit:** cm

**Values:** Greater than 0

**Parent(s):**

* :ref:`Wind.type`: corona


**File:** corona.c


Corona.vel_frac
===============
For the coronal model, the azimuthal velocity is
given by the velocity of the underlying disk.  One
can also give the corona a radial velocity, which is
a fraction of the disk velocity.  (As coded, if this
number is positive, the velicty is the r direction is
toward the central object).

**Type:** Double

**Unit:** Disk velocity

**Values:** Any, 0 implies no radial velocity.

**Parent(s):**

* :ref:`Wind.type`: corona


**File:** corona.c


