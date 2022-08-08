# Maternal-Transmission
Supplemental Code to accompany "Maternal transmission as a symbiont sieve, and the absence of lactation in male mammals", by Brennen Fagan, George W. A. Constable, and Richard Law.

Four files are included:

1. C-README.txt, by Richard Law
2. C-Simulations.c, by Richard Law
3. Mathematica-DynSys.nb, by Brennen T. Fagan
4. Mathematica-Figure5.nb, by Brennen T. Fagan

C-README.txt includes the descriptions of the parameters required to use C-Simulations.c to recreate the images in the main text. After making any required modifications, you will need to compile C-Simulations.c and run the resulting output file. Note that demographic stochasticity can be a problem, which is why we provide random seeds in some locations.

Mathematica code was last edited using Mathematica 12.1.1.0. We adopt the convention therein of referring to the populations by their (mean-field) densities. Variable names are generally the same as in the main text except as follows:

1. x+ -> x
2. x- -> y
3. e0 -> m

The file Mathematica-DynSys.nb is the tool we used to analyse and robustness check the system. It comprises an implementation of the dynamical system as a manipulate object, which allows the user to dynamically change the parameters used. Note that there are three methods to dynamically change the parameters of the system. Initial conditions can be specified by clicking on the phase portait. Individual parameters can be changed by moving the sliders as well as by pressing the "+" button on the right of a slider for more precise manipulation. Additionally, we display all detected fixed points above the phase portrait.

The file Mathematica-Figure5.nb is the tool we used to categorise fixed points as a function of the parameters. Helper functions are defined at the top of the file. The main function is below and recreates a Mathematica version of Figure 5 in the main text. The Mathematica version can be manipulated by re-running the module with different parameter values. The resultant chart will report the birth and death parameters used to construct it (above), as well as all fixed point combinations identified (by color, below). Note also that the chart is plotted "pixel-by-pixel"; zooming in will separate out the pixels. The number of pixels plotted can be increased by changing the step sizes of m and w (in the module declaration).

