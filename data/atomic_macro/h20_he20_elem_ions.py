

# Output format: Name z  Relative_abun atomic_weight


# Source Verner, Barthel & Tytler 1994 A&AS 108, 287
Element    1    H   12.00     1.007940
Element    2   He   10.99     4.002602

# IonM  element.name element.z ionstate g_ground ionizpotential(eV)  maxlevels maxnlte   Config
# maxlevels and maxnlte set at 30 and 0 respectively are simply place holders


IonM    H   1   1   2   13.59900   20  20     1s(2S_{1/2}) 
IonM    H   1   2   1 1.0000e+20    1   1     Bare 

# the 20 is the number of "non-LTE" levels to be used - this is infact the
# number of levels to be treated (for my purposes). 

IonM   He   2   1   1   24.58800   11  11     1s^2(1S_0)
IonM   He   2   2   2   54.41800   10  10     1s(2S_{1/2})
IonM   He   2   3   1 1.0000e+20    1   1     Bare
 


