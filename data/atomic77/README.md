This data set is the same as standard73 except for the extrapolated cross sections.

The procedure is pretty klugey:

* identify odd XSections by eye
* Check none of them are the ground state
* If the gradient is above ~3 (I used 2.9 to give a slight buffer) and it is in the 'odd looking' list, then fudge it by setting gradient to -3
* if not, use the gradient at the last point in the data
* extrapolate up to 100Kev



It is prepared using extrapolate_log.py and subroutines in extrapolate_sub.py
