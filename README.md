# Maternal-Transmission
Supplemental Code to accompany "Maternal transmission as a symbiont sieve, and the absence of lactation in male mammals", by Brennen Fagan, George W. A. Constable, and Richard Law.

Four files are included:

1. C-README.txt, by Richard Law
2. C-Simulations.c, by Richard Law
3. C_code_20spp.c, By Richard Law
4. Mathematica-DynSys.nb, by Brennen T. Fagan
5. Mathematica-Figure5.nb, by Brennen T. Fagan
6. Maternal_transmission_SI_code.nb, by George W. A. Constable

C-README.txt includes the descriptions of the parameters required to use C-Simulations.c to recreate the images in the main text. After making any required modifications (e.g. changing parameter values), you will need to compile C-Simulations.c and run the resulting output file. Note that demographic stochasticity will in theory lead to different results with each realisation, which is why we provide random seeds in some locations (i.e. to retain consistency with the specific results plotted in the paper).

C_code_20spp.c is a modification of C-Simulations.c that expands the effective number of taxa in the microbiome. The taxa are again assumed to be independent of one another, so that the host fitness is 1 if there are no symbionts and the average of the symbiont induced fitnesses if any symbionts are present. To recreate Figure 3 in the paper, we use random seed 7597.

Mathematica code by Brennen T. Fagan was last edited using Mathematica 12.1.1.0. We adopt the convention therein of referring to the populations by their (mean-field) densities. Variable names are generally the same as in the main text except as follows:

1. x+ -> x
2. x- -> y
3. e0 -> m

The file Mathematica-DynSys.nb is the tool we used to analyse and robustness check the system. It comprises an implementation of the dynamical system as a manipulate object, which allows the user to dynamically change the parameters used. Note that there are three methods to dynamically change the parameters of the system. Initial conditions can be specified by clicking on the phase portait. Individual parameters can be changed by moving the sliders as well as by pressing the "+" button on the right of a slider for more precise manipulation. Additionally, we display all detected fixed points above the phase portrait.

The file Mathematica-Figure5.nb is the tool we used to categorise fixed points as a function of the parameters. Helper functions are defined at the top of the file. The main function is below and recreates a Mathematica version of Figure 5 in the main text. The Mathematica version can be manipulated by re-running the module with different parameter values. The resultant chart will report the birth and death parameters used to construct it (above), as well as all fixed point combinations identified (by color, below). Note also that the chart is plotted "pixel-by-pixel"; zooming in will separate out the pixels. The number of pixels plotted can be increased by changing the step sizes of m and w (in the module declaration).

The file Maternal_transmission_SI_code.nb contains the code that we used to generate the figures in our supplemental information. This supplemental information examines the implications of the benefits of lactation on the population by separately considering the benefits as offspring survival probability and risks as symbiont transmission probability. In the first section, a population with uniparental lactation is subject to invasion by a small biparental invasive population. In the second section, these roles are reversed. In both sections, the benefits and risks are varied between subsections and the volume of milk of the invasive population and danger of an environmentally acquired symbiont (that is nonetheless transmissable via lactation) are plotted against the resulting fraction of males that are uniparental or experience reduced lactation.
