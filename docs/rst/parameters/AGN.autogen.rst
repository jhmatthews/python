===
AGN
===

AGN.power_law_cutoff
====================
An advanced option - this is a low frequency cutoff for an
AGN power law spectrum. It prevents the powerlaw being
applied to low frequencies and giving an odd SED.

**Type:** Double

**Unit:** Hz

**Values:** Greater than 0

**Parent(s):**

* :ref:`BH.radiation`: ``True``


**File:** setup_star_bh.c


AGN.power_law_index
===================
The exponent alpha in a power law SED applied to an AGN
central source of the form L_nu=K nu**alpha

**Type:** Double

**Values:** Any - but sign is not assumed, so for negative index use a negative value

**Parent(s):**

* :ref:`BH.radiation`: ``True``


**File:** setup_star_bh.c


AGN.bremsstrahlung_alpha
========================
The frequency exponent alpha in bremstrahlung SED of the form
L_nu=nu**alpha exp(-hnu/kT)

**Type:** Double

**Values:** Any - sign is not assumed so use negative if you want negative

**Parent(s):**

* :ref:`BH.rad_type_to_make_wind`: brems


**File:** setup_star_bh.c


AGN.bremsstrahlung_temp
=======================
The temperature T in bremstrahlung SED of the form
L_nu=nu**alpha exp(-hnu/kT)

**Type:** Double

**Unit:** K

**Values:** Greater than 0

**Parent(s):**

* :ref:`BH.rad_type_to_make_wind`: brems


**File:** setup_star_bh.c


AGN.lamp_post_height
====================
Multi-line description, must keep indentation.

**Type:** Double

**Unit:** :ref:`Central_object.radius`

**Values:** Greater than 0

**Parent(s):**

* :ref:`BH.geometry_for_pl_source`: lamp_post


**File:** setup_star_bh.c


