=====
Shell
=====

Shell.wind_v_at_rmin
====================
The velocity of a shell wind at the inner edge of the
shell - the variation of the velocity in the shell is
set by the velocity law exponent. It allows a gradient
to be enforced.

**Type:** Double

**Unit:** cm/s

**Values:** Greater than or equal to 0

**Parent(s):**

* :ref:`Wind.type`: shell


**File:** shell_wind.c


Shell.wind.radmax
=================
Multi-line description, must keep indentation.

**Type:** Double

**Unit:** cm

**Values:** Greater than :ref:`Shell.wind.radmin`

**Parent(s):**

* :ref:`Wind.type`: shell


**File:** shell_wind.c


Shell.wind_mdot
===============
The mass loss through a diagnostic shell type wind. One normally sets
this experimentally in order to get a required hydrogen density in
the shell

**Type:** Double

**Unit:** M☉/year

**Values:** Greater than 0

**Parent(s):**

* :ref:`Wind.type`: shell


**File:** shell_wind.c


Shell.wind.v_at_rmax
====================
The velocity of a shell wind at the outer edge of the
shell - the variation of the velocity in the shell is
set by the velocity law exponent. It allows a gradient
to be enforced.

**Type:** Double

**Unit:** cm/s

**Values:** Greater than or equal to 0

**Parent(s):**

* :ref:`Wind.type`: shell


**File:** shell_wind.c


Shell.wind.radmin
=================
The innermost edge of a diagnostic type of wind made up of a single
(ideally thin) shell.

**Type:** Double

**Unit:** cm

**Values:** Greater than 0

**Parent(s):**

* :ref:`Wind.type`: shell


**File:** shell_wind.c


Shell.wind.acceleration_exponent
================================
Exponent beta for the Caster and Lamers description of a stellar wind
v(r)=v_o + (v_inf - v_o) (1+R_s/r)**beta for a shell wind.

**Type:** Double

**Values:** Greater than or equal to 0

**Parent(s):**

* :ref:`Wind.type`: shell


**File:** shell_wind.c


