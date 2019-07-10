====
Diag
====

Diag.write_atomicdata
=====================
Choose whether to write the atomic data that is being used to
an output file.

**Type:** Boolean (yes/no)

**File:** setup_domains.c


Diag.adjust_grid
================
Choose whether or not you would like to adjust the scale length
for the logarithmic grid. Advanced command.

**Type:** Boolean (yes/no)

**Parent(s):**

* :ref:`Wind.number_of_components`: Greater than 0. Once per domain.


**File:** setup_domains.c


Diag.extra
==========
Decide whether or not to use extra diagnostics in advanced mode.
If set, this triggers a many extra questions that allow one to investigate
things such as photon cell statistics, the velocity gradients in cells and
the resonant scatters in the wind

**Type:** Boolean (yes/no)

**File:** python.c


Diag.save_cell_statistics
-------------------------
Choose whether to save the statistics for a selection of save_cell_statistics.
If yes, it looks for a file called "diag_cells.dat" which contains the cells to track,
and saves the photon details (weights, frequencies) for those that interact in
the cell. Useful for checking the detailed MC radiation field in a cell.

**Type:** Boolean (yes/no)

**Parent(s):**

* :ref:`Diag.extra`: ``True``


**File:** diag.c


Diag.make_ioncycle_tables
-------------------------
Multi-line description, must keep indentation.

**Type:** Boolean (yes/no)

**Parent(s):**

* :ref:`Diag.extra`: ``True``


**File:** diag.c


Diag.print_dvds_info
--------------------
Print out information about the velocity gradients in the
cells to a file root.dvds.diag.

**Type:** Boolean (yes/no)

**Parent(s):**

* :ref:`Diag.extra`: ``True``


**File:** diag.c


Diag.save_photons
-----------------
Multi-line description, must keep indentation.

**Type:** Boolean (yes/no)

**Parent(s):**

* :ref:`Diag.extra`: ``True``


**File:** diag.c


Diag.save_extract_photons
-------------------------
Multi-line description, must keep indentation.

**Type:** Boolean (yes/no)

**Parent(s):**

* :ref:`Diag.extra`: ``True``


**File:** diag.c


Diag.keep_ioncycle_windsaves
----------------------------
Decide whether or not to keep a copy of the windsave file after
each ionization cycle in order to track the changes as the
code converges. Produces files of format python01.wind_save and so
on (02,03...) for subsequent cycles.

**Type:** Boolean(yes/no)

**Parent(s):**

* :ref:`Diag.extra`: ``True``


**File:** diag.c


Diag.track_resonant_scatters
----------------------------
Multi-line description, must keep indentation.

**Type:** Boolean (yes/no)

**Parent(s):**

* :ref:`Diag.extra`: ``True``


**File:** diag.c


Diag.use_standard_care_factors
==============================
Advanced command which allows one to change
various other defaults associated with
radiative transfer, inclusing the fractional distance
in a cell that a photon can travel

**Type:** Boolean (yes/no)

**File:** diag.c


Diag.lowest_ion_density_for_photoabs
------------------------------------
For efficiency reasons, Python does not try to calculate photoabsorption
for an ion with an extremly low density.  This advance parameter changes
this density limit

**Type:** Double

**Unit:** n/cm**3

**Values:** Greater than 0

**Parent(s):**

* :ref:`Diag.use_standard_care_factors`: ``False``


**File:** diag.c


Diag.keep_photoabs_in_final_spectra
-----------------------------------
This advanced options allows you to include or exclude photoabsorpiotn
in calculating the final spectra.  (but ksl does not know what the
default is)

**Type:** Boolean (yes/no)

**Parent(s):**

* :ref:`Diag.use_standard_care_factors`: ``False``


**File:** diag.c


Diag.fractional_distance_photon_may_travel
------------------------------------------
The distance photon may travel in a cell is limited to prevent a photon
from moving such a long path that the velocity may change non-linearly.
This problem arises primarily when the photon is travelling azimuthally
in the grid.  This changes the default for the fraction of the maximum
distance in a cell.

**Type:** Double

**Values:** 0 to 1

**Parent(s):**

* :ref:`Diag.use_standard_care_factors`: ``False``


**File:** diag.c


