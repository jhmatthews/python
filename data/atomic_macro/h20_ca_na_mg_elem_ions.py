

# Output format: Name z  Relative_abun atomic_weight


# Source Verner, Barthel & Tytler 1994 A&AS 108, 287
Element    1    H   12.00     1.007940
Element    11    Na   6.31    22.989768
Element   12   Mg    7.58    24.305000
Element    20    Ca   6.34      40.078

# IonM  element.name element.z ionstate g_ground ionizpotential(eV)  maxlevels maxnlte   Config
# maxlevels and maxnlte set at 30 and 0 respectively are simply place holders


IonM    H   1   1   2   13.59900   20  20     1s(2S_{1/2}) 
IonM    H   1   2   1 1.0000e+20    1   1     Bare 

# the 20 is the number of "non-LTE" levels to be used - this is infact the
# number of levels to be treated (for my purposes). 
IonM   Na   11   1   2   5.13907    20  20     3s(2S_{1/2})
IonM   Na   11   2   1 1.0000e+20    1   1     Bare

IonM   Mg   12   1   1   7.64623    15  15     3s2(1S_{0})
IonM   Mg   12   2   2  15.035266   18  18     3s(2S_{1/2})
IonM   Mg   12   3   1 1.0000e+20    1   1     Bare

IonM   Ca   20   1   1   6.11316    15  15     4s2(1S_{0})
IonM   Ca   20   2   2   11.87      24  24     4s(2S_{1/2})
IonM   Ca   20   3   1 1.0000e+20    1   1     Bare
 


