Parameter settings to repeat stochastic results in paper
--------------------------------------------------------

The code was run using the operating system Ubuntu 18.04.6 LTS. The code 
is written in C and was compiled with gcc supplied with the operating system
gcc (Ubuntu 7.5.0-3ubuntu1~18.04) 7.5.0.  To reproduce the results in the 
figures, the code needs to be compiled with the settings below.

Figure 2b  This uses the default settings of the code.
           Output: frequency distrbution of wsym in output.freq.dist.dat

Figure 2c  Line 30: set method_vertical = 0
           Output: frequency distrbution of wsym in output.freq.dist.dat

Figure 2d,e
           Line   30: for biparental transmission, set method_vertical = 2
                      for maternal   transmission, set method_vertical = 0
           Line   42: set tmax = 10
           Line   53: set niterate = 1
           Line  150: set flag_timeseries = 1
           Line 1265 to 1267: comment out these lines
           Line 1269: remove comment brackets and set wsym to values in the figure
           Output: timeseries in output.timeseries.dat

Figure 3   See C_code_20spp.c with random seed 7597.
           Line   41: for biparental transmission, set method_vertical = 2
                      for maternal   transmission, set method_vertical = 0

Figure 4   Line   30: set method_vertical = 3
           Line   42: set tmax = 20
           Line   44: set Nsymstart = 100
           Line   46: set fmodstart = 0.2
           Line   53: set niterate = 1
           Line  150: set flag_timeseries = 1
           Line 1265 to 1267: comment out these lines
           Line 1269: remove comment brackets and set wsym = 0.6
           Random number seed: 725
           Output: timeseries in output.timeseries.dat
           More than one realization may be needed to get an instance of
           complete fixaton of the the modifier gene M^-

Figure 5c  Line   30: set method_vertical = 3
           Line   36: set method_horizontal = 1
           Line   42: set tmax = 100
           Line   43: set Nstart = 656
           Line   44: set Nsymstart = 80
           Line   46: set fmodstart = 0.02
           Line   53: set niterate = 1
           Line  150: set flag_timeseries = 1
           Line 1265 to 1267: comment out these lines
           Line 1269: remove comment brackets and set wsym = 0.37
           Random number seed: 470
           Output: timeseries in output.timeseries.dat
           More than one realization may be needed to get an instance of
           complete fixaton of the the modifier gene M^-


