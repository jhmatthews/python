==========
Homologous
==========

Homologous.radmin
=================
The starting point of for madel of a homologous flow, a model in
which the velocity at any radius is proportional to the radius

**Type:** Double

**Unit:** cm

**Values:** Greater than or equal to :ref:`Central_object.radius`

**Parent(s):**

* :ref:`Wind.type`: homologous


**File:** homologous.c


Homologous.radmax
=================
Maximum extent of the homologous wind.

**Type:** Double

**Unit:** cm

**Values:** Greater than :ref:`Homologous.radmin`

**Parent(s):**

* :ref:`Wind.type`: homologous


**File:** homologous.c


Homologous.density_exponent
===========================
The power law exponent which defines the decline in density of
a homologous flow as a function of radious.

**Type:** Double

**Values:** Greater than 0 for a density that declines with radius

**Parent(s):**

* :ref:`Wind.type`: homologous


**File:** homologous.c


Homologous.vbase
================
Velocity at the base of the wind

**Type:** Double

**Unit:** cm

**Values:** Greater than 0

**Parent(s):**

* :ref:`Wind.type`: homologous


**File:** homologous.c


Homologous.boundary_mdot
========================
The mass loss rate at the base of the wind in a homlogous flow model, a flow
in which the velocity is proporional to the radius.  In general, mdot will
decline with radius, depending on the exponent of the power law that describes
the trend in density.

**Type:** Double

**Unit:** M☉/yr

**Values:** Greater than 0

**Parent(s):**

* :ref:`Wind.type`: homologous


**File:** homologous.c


