/****************************************************************************************************************
* Updated C code from symbiont-sieve paper									*
*	                                                                                                        *
* This allows the introduction ofmultiple symbiont taxa                                                         *
*                                                                                                               *
* For maternal   transmission, set method_vertical: 0                                                           *
* For biparental transmission, set method_vertical: 2                                                           *
*                                                                                                               *
* Output used in Figure 3 is in output.timeseries.dat                                                           *
*                                                                                                               *
* 14.05.23                                                                                                      *
****************************************************************************************************************/


#include    <stdlib.h>
#include    <math.h>
#include    <stdio.h>
#include    <string.h>


/****************************************************************************************************************
* Constants                                                                                                     *
****************************************************************************************************************/

#define     error    1e-9          /* rounding error allowed for type double                         */
#define     pi       3.1415927

#define     method_microbiome 2     /* switch for:
                                       1: single   symbiont
                                       2: multiple symbiont taxa                                     */

#define     method_fitness   1      /* effect of microbiome taxa on host fitness
                                       1: additive:         mean of single-taxon w's                 */

#define     method_sym       4      /* symbiont effect on host
                                       1: whost acts on intrinsic birth rate: b0*w -- removed
                                       2: whost acts on intrinsic death rate: d0*w -- removed
                                       3: whost acts on intrinsic death rate: d0/w -- removed
                                       4: whost acts as a factor on total     death rate: d/whost    */

#define     method_vertical  0      /* method of vertical transmission of symbionts
                                       0: maternal
                                       1: paternal    (not implemented)
                                       2: biparental: transmitted if one or both parents have it
                                       3: evolving:   modifier gene on Y chromosome: 0 mat; 1 bip    */

#define     method_horiz     0      /* method of horizontal transmission of symbionts
                                       0: off
                                       1: from environment
                                       2: from another host (as in epidemics -- not implemented)     */

#define     nspp             2      /* number of host species + 1                                    */
#define     Nmic	    21      /* maximum number of taxa in microbiome + 1                      */

#define     tmax          10.0      /* maximum time for simulations                                  */
#define     Nstart        1000      /* initial population size of host species (must be > Nsymstart) */
#define     Nsymstart       50      /* initial number of sym+ hosts (use: 0 to check a sym- popn)    */

#define     fmodstart     0.02      /* Starting frequency for gene stopping male transmission applies
                                       only if there is evolution, i.e. method_vertical: 3:
                                       fmodstart ~ 0.0: popn starts close to bip transmission
                                       fmodstart ~ 1.0: popn starts close to mat transmission.
                                       Modifier gene assumed to be on Y chromosome, passed on to and
                                       expressed only in sons: carried only in males, not in females */

#define     niterate        1       /* number of iterations of stochastic process                    */
#define     nbin            50      /* size of array with Wsym frequency distribution                */
#define     binwidth        0.05    /* binwidth for Wsym frequency distribution                      */
#define     fsuccess        0.1     /* threshold freq of sym+ hosts to count as successful invasion;
                                       distinct from the stopping condition in simulate()            */


/****************************************************************************************************************
* Global types                                                                                                  *
****************************************************************************************************************/

/***A structure to hold the basic information on a single individual one host species at a time***/
/***Designed for a bidirectional list in each host population***/
struct  rep
        {   int     sp;                  /* host species                                             */
            int     sex;                 /* host sex                0: female  1: male               */
            int     sym;                 /* host symbiont           0: absent  1: present            */
            int     mic  [Nmic];         /* host symbiont vector    0: sym abs 1: sym present        */
            double  wmic [Nmic];         /* symbiont fitness vector 1; sym abs w: sym pres           */
            double  w1mic[Nmic];         /* symbiont fitness vector 1; sym abs w1:sym pres (singles) */
            double  w2mic[Nmic];         /* symbiont fitness vector 0; sym abs w2:sym pres (pairs)   */
                                         /* this intrinsic effect on host needs to be combined for symbiont interactions on host   */
            int     mod;                 /* modifier gene male transm  0: off(mat) 1: on(bip)        */
            double  whost;               /* fitness of host allowing for symbiont                    */
            double  e;                   /* host prob p.u.t of encounter with symbiont in environment*/
            double  b;                   /* host prob p.u.t of birth                                 */
            double  d;                   /* host prob p.u.t of death                                 */
            struct  rep  *nextP, *prevP; /* pointers for bidirectional list                          */
        };

typedef struct rep  REP;            /* REP is a new data type for declaring structure of type rep    */
typedef        REP  *repP;          /* *repP is pointer to the data type REP.                        */


/***A structure to hold information on properties of a host species***/
struct  species_properties
        {   int     n;              /* number of hosts (total)	                          */
            int     nf;             /* number of hosts (female)                           */
            int     nm;             /* number of hosts (male)                             */

            int     nsym;           /* number of hosts with symbiont                      */
            double  fsym;           /* freqcy of hosts with symbiont                      */

            int     nmic [Nmic];    /* number of hosts with each symbiont (microbiome)    */
            double  fmic [Nmic];    /* freqcy of hosts with each symbiont (microbiome)    */
            double  w1mic[Nmic];    /* fitnes from symbiont singletons    (microbiome)    */
            double  w2mic[Nmic];    /* fitnes from symbiont pairs         (microbiome)    */

            int     nfsym0;         /* number of hosts (female sym-)                      */
            int     nmsym0;         /* number of hosts (male   sym-)                      */
            int     nfsym1;         /* number of hosts (female sym+)                      */
            int     nmsym1;         /* number of hosts (male   sym+)                      */

            int     nmat;           /* number of transmit- males (trans maternal)         */
            int     nbip;           /* number of transmit+ males (trans biparental)       */
            double  fmat;           /* frequency of transmit- gene                        */

            int     nsym0mat;       /* number of sym- transmit- males                     */
            int     nsym0bip;       /* number of sym- transmit+ males                     */
            int     nsym1mat;       /* number of sym+ transmit- males                     */
            int     nsym1bip;       /* number of sym+ transmit+ males                     */
            double  dis;            /* coefficient of disequilibrium                      */

            double  e;              /* prob p.u.t of encounter with symbiont in environment*/

            double  E0;             /* prob p.u.t of symbiont encounter in environment     */
            double  B0;             /* prob p.u.t.of birth (intrinsic)                     */
            double  D0;             /* prob p.u.t of death (intrinsic)                     */
            double  DD;             /* prob p.u.t of death (density dependent component)   */
                                    /* (no L^2: using number of individuals)               */
            repP    firstP;         /* pointer to first individual                         */
            repP    lastP;          /* pointer to last  individual                         */
        };


/***Holds a pointer to the first element of lists of each species***/
struct  rep state[nspp];

/***Holds properties of the species***/
struct  species_properties	param[nspp];


/****************************************************************************************************************
* Global variables                                                                                              *
****************************************************************************************************************/

double  t;                          /* time                                                */
int     nevent;                     /* counter for number of events                        */

/***Variables for random number generator***/
double  drand48();                  /* Function for random number generator                */
int     seed;                       /* Seed for random number generator                    */

double  Wsym;                       /* effect of symbiont on host fitness; implemented as
                                       a factor by which intrinsic rate is multiplied
                                         < 1:  host rate decreases
                                         = 1:  host rate unchanged
                                         > 1:  host rate increases                         */
double  freq[2][nbin];              /* array with freq distrib of Wsym's allowing invasion */

/***Flags to switch on detailed output for checking code and running single realization***/
int     flag_state      = 0;        /* 1: output for system state on;              0: off  */
int     flag_choose     = 0;        /* 1: output for choose event calculations on  0: off  */
int     flag_choosedad  = 0;        /* 1: output for choose event calculations on  0: off  */
int     flag_vertical   = 0;        /* 1: output for symbiont vertical transm  on  0: off
                                       removed this output -- it's the same as in birth()  */
int     flag_horizontal = 0;        /* 1: output for symbiont horizontl transm on  0: off  */
int     flag_birth      = 0;        /* 1: output for birth calculations on         0: off  */
int     flag_death      = 0;        /* 1: output for death calculations on         0: off  */
int     flag_timeseries = 1;        /* 1: output for timeseries of realization on  0: off  */


/****************************************************************************************************************
*  Function for generating random numbers on a normal distribution N(0,1)                                       *
****************************************************************************************************************/
double  random_normal()
{
        int i;
        double theta, y, value1, value2;

        y     = drand48();
        y     = -2.0 * log(y);
        theta = drand48();
        theta = theta * 2.0 * pi;

        value1 = sqrt(y) * cos(theta);        /*only this one being used*/
        value2 = sqrt(y) * sin(theta);

        return(value1);
}


/****************************************************************************************************************
* Input parameter values                                                                                        *
****************************************************************************************************************/

void    initialize_parameters()
{
        FILE    *out;
        int     sp;

    /***Enter host species population parameters***/
        sp = 1;
        param[sp].n  = Nstart;
        param[sp].B0 = 4.0;
        param[sp].D0 = 1.0;
        param[sp].DD = 0.001;

    /***Rule for horizontal transmission of symbionts***/
        switch(method_horiz)
        {       case 0: /*horizontal transmission of symbionts switched off*/
                          param[sp].E0 = 0.0;
                          break;
                case 1: /*horizontal transmission of symbionts from environment*/
                          param[sp].E0 = 0.1;
                          break;
                case 2: /*horizontal transmission of symbionts from other hosts (as in epidemics)*/
                          printf("method_horiz:%d epidemic model still TO DO -- program halted\n", method_horiz);
                          exit(1);
                          break;
        }

    /***Record of parameter values***/
        out = fopen("output.detail.dat", "a");
        fprintf(out, "HOST POPULATION PARAMETERS\n");
        sp = 1;
        fprintf(out, "  species: %d\n",         sp);
        fprintf(out, "  param[sp].n:    %5d\n", param[sp].n);
        fprintf(out, "  param[sp].E0: %7.3f\n", param[sp].E0);
        fprintf(out, "  param[sp].B0: %7.3f\n", param[sp].B0);
        fprintf(out, "  param[sp].D0: %7.3f\n", param[sp].D0);
        fprintf(out, "  param[sp].DD: %7.3f\n", param[sp].DD);
        fprintf(out, "  Nhat = (B0*0.5 - D0) / DD  (0.5 because only females give birth)\n");

        fprintf(out, "\nREALISATIONS\n");
        fprintf(out, "  tmax:            %4.0f    time over which stochastic realisation runs\n", tmax);
        fprintf(out, "  nspp-1:          %4d	   number of host species\n",                     nspp-1);
        fprintf(out, "  Nstart:          %4d	   initial population size of host spp\n",        Nstart);

        fprintf(out, "\nSYMBIONTS\n");
        fprintf(out, "  method_microbio: %4d    # symbiont taxa: 1 single taxon; 2 multiple taxa\n",               method_microbiome);
        fprintf(out, "  method_fitness:  %4d    host fitness:    1 additive;     2 pairs          3 ring of pairs\n", method_fitness);
        fprintf(out, "  method_vertical: %4d    vertical transm: 0 maternal;     2 biparental;    3 evolving\n",     method_vertical);
        fprintf(out, "  method_horiz:    %4d    horizont transm: 0 absent;       1 environmental; 2 from hosts\n",      method_horiz);
        fprintf(out, "  method_sym:      %4d    symbionts affect host as:        4 D/w\n",                                method_sym);
        fprintf(out, "  Nmic:            %4d    number of symbiont taxa in microbiome\n",         Nmic-1);
        fprintf(out, "  Nsymstart:       %4d    initial number of hosts with each symbiont\n", Nsymstart);
/*      fprintf(out, "  fsuccess:        %4.2f    threshold freq for successful invasion\n",    fsuccess);
*/
        if(method_vertical==3)
        {
        fprintf(out, "\nMALE TRANSMISSION MODIFIER\n");
        fprintf(out, "  fmodstart:     %6.4f	initial frequency of no-male-transmit modifier\n",	fmodstart);
        }

        fprintf(out, "\nSTOCHASTIC REALISATIONS\n");
        fprintf(out, "  niterate:        %4d    number of iterations of stochastic process\n", niterate);
        fprintf(out, "  nbin:            %4d    number of bins for Wsym frequency distribution\n", nbin);
        fprintf(out, "  binwidth:        %4.2f    binwidth for Wsym frequency distribution\n", binwidth);
        fprintf(out, "\n");
        fflush(out);
        fclose(out);
}


/****************************************************************************************************************
* This procedure outputs the current state of the system                                                        *
****************************************************************************************************************/

void    output_system_state()
{
        FILE    *out;
        int     sp, sm, i, j;
        repP    nowP;


    /***Output bidirectional list of host individuals by species***/
        out = fopen("output.detail.dat", "a");
        fprintf(out, "\n");
        for (sp=1; sp<nspp; ++sp)
        {       nowP = param[sp].firstP;
                fprintf(out, "Time:%9.6f nevent:%d sp:%d   n:%d   nf:%d   nm:%d\n", 
                              t, nevent, sp, param[sp].n, param[sp].nf, param[sp].nm);
                fprintf(out, "Species:%d   n:%d   nf:%d   nm:%d\n", sp, param[sp].n, param[sp].nf, param[sp].nm);
                fprintf(out, "   i       prevP        nowP       nextP  sex    b      d      e      modifier  whost   mic\n");
                i = 1;
                while (nowP != NULL)
                {       fprintf(out, "%4d %11d %11d %11d    %d   %6.3f %6.3f %6.3f     %2d    %6.3f   ",
                                i, nowP->prevP, nowP, nowP->nextP, 
                                nowP->sex, nowP->b, nowP->d, nowP->e, nowP->mod, nowP->whost);
                        for (sm=1; sm<Nmic; ++sm)  fprintf(out, "%d",     nowP->mic  [sm]);   fprintf(out, "  w1: ");
                        for (sm=1; sm<Nmic; ++sm)  fprintf(out, "%4.2f ", nowP->w1mic[sm]);
                        if (method_fitness ==2 || method_fitness ==3)
                        {       fprintf(out, "  w2: ");
                                for (sm=1; sm<Nmic; ++sm)  fprintf(out, "%4.2f ", nowP->w2mic[sm]);
                        }
                        fprintf(out, "\n");
                        i = i + 1;
                        nowP = nowP->nextP;
                }
        }
        fflush(out);
        fclose(out);
}


/****************************************************************************************************************
*  Procedure to construct fitness pickP->whost of a single host pickP, given its microbiome                     *
****************************************************************************************************************/
void    host_fitness_microbiome(pickP)
repP	pickP;
{
        int     countsym, sm;
        double  toadd;

    /***Methods to construct the effect of symbiont tax on host fitness***/
    /***The number ofsymbiont taxa can change over time (biparental), so # taxa averaged changes as well***/
	switch(method_fitness)
	{
                case 1: /*additive: mean of single-taxon w's*/
                        pickP->whost = 0.0;
                        countsym     = 0;
                        for (sm=1; sm<Nmic; ++sm)
                        if  (pickP->mic[sm] == 1)
                        {       pickP->whost += pickP->w1mic[sm];
                                countsym    += 1;
	                }
                        if (countsym == 0)   pickP->whost = 1;    
                        if (countsym >  0)   pickP->whost = pickP->whost/countsym;
                        break;

        }
}


/****************************************************************************************************************
* This procedure initializes properties of individual hosts at start of simulate				*
* Microbiome version: allows single and multiple symbiont taxa per host                                         *
****************************************************************************************************************/

void    initialize_hosts_microbiome()
{
        int     sp, sm, i, j;
        double  rn, add, W1sym, W2sym;
        repP    firstP, lastP, nowP, prevP;


        sp = 1;

    /***Set the size of host population***/
        param[sp].n = Nstart;
        if (param[sp].n < 1)
        {       printf("n:%d  Must have individuals in the population -- program halted\n", param[sp].n);
                exit(1);
        }

    /***Check number of hosts carrying symbiont***/
        if (Nstart < Nsymstart)
        {       printf("Nstart:%d  Nsymstart:%d not a good starting point -- program halted\n", Nstart, Nsymstart);
                exit(1);
        }


/***BIDIRECTIONAL LIST OF ALL HOSTS USING POINTERS***/

    /***Set the first and last pointers to null***/
        param[sp].firstP = NULL;
        param[sp].lastP  = NULL;

    /***Set the start pointer as long as individuals actually exist***/
        param[sp].firstP = &state[sp];
        prevP = NULL;
        nowP  = param[sp].firstP;

    /***An extra individual '0' is created and removed below, so all memory allocation is by malloc***/
        for (i=0;   i<=param[sp].n;  ++i)
        {
            /***Allocate the species of the individual***/
                nowP->sp = sp;

            /***Sex of the individual***/
                rn = drand48();
                if(rn <  0.5) nowP->sex = 0; /* female */
                if(rn >= 0.5) nowP->sex = 1; /* male   */

            /***Initialise the microbiome with all symbiont taxa absent (0), and effect on host fitness at 1.0***/
                for (sm=1; sm<Nmic; ++sm)   nowP->  mic[sm] = 0;
                for (sm=1; sm<Nmic; ++sm)   nowP->w1mic[sm] = 1.0;
                for (sm=1; sm<Nmic; ++sm)   nowP->w2mic[sm] = 1.0;

            /***Evolution of transmission (case 3): males have modifier gene (dummy value for females: gene on Y chr***/
                switch (method_vertical)
                {       case 0: /*maternal: male transmission always off*/
                                nowP->mod = 99;
                                break;
                        case 1: /*paternal (not implemented)*/
                                break;
                        case 2: /*biparental: male transmission always on*/
                                nowP->mod = 99;
                                break;
                        case 3: /*modifier gene for  male transmission: 0:off maternal, 1:on biparental*/
                                /*mod not defined for females; dummy value inserted*/
                                if (nowP->sex==0)               nowP->mod = 99;
                                if (nowP->sex==1)
                                {       rn = drand48();
                                        if (rn <  fmodstart)    nowP->mod = 0;
                                        else                    nowP->mod = 1;
                                }
                                break;
                }

            /***Create the pointer links to previous and next individual***/
                nowP->prevP = prevP;
                nowP->nextP = NULL;
                if (i != param[sp].n)
                {       nowP->nextP = (repP) malloc(sizeof(REP));
                        if (!nowP->nextP)
                        {       printf("Unable to allocate memory to pointer: exiting from program\n");
                               exit(1);
                        }
                }

            /***Update the pointers***/
                prevP = nowP;
                nowP  = nowP->nextP;
        }
        param[sp].lastP = prevP;

    /***Remove individual '0' so all individuals used are created by malloc***/
        nowP             = param[sp].firstP;
        nowP             = nowP->nextP;
        param[sp].firstP = nowP;
        nowP->prevP      = NULL;


/***CONSTRUCT MICROBIOME***/

    /***Symbionts independently distributed across hosts***/
    /***# symbiont taxa: Nmic-1;  # hosts: Nstart;  # symbionts of each taxon at start: Nsymstart)***/
    /***All the work with a single symbiont taxon in first submission to Nature Comms inserted as the case where Nmic-1=1***/
        for (sm=1; sm<Nmic; ++sm)
        {
                switch (method_microbiome)
		{       case 1: /*single symbiont: Wsym set in main() (as in first submission to Nature Comms)*/
                                param[sp].w1mic[1] = Wsym;
                                param[sp].w2mic[1] = 0.0;
                                break;
                        case 2: /*multiple symbiont taxa with random W1sym: N(1,0.3): <= W1sym <= 2.5*/
                                W1sym = -1.0;
                                while (W1sym < error || W1sym > 2.5)
                                W1sym = random_normal()*0.3 + 1.0;
                                param[sp].w1mic[sm] = W1sym;
                                break;
		}

            /***Parameter values for pair interactions, W2sym; only used when multiple symbiont taxa are present***/
                if (method_microbiome == 2)
                switch (method_fitness)
                {       case 1: /*additive: w2mic not used*/
                                param[sp].w2mic[sm] = 0.0;
                                break;
                        case 2: /*pair interaction*/
                                W2sym = -1.0;  
                                while (W2sym<error || W2sym>2.0) 
                                W2sym = random_normal()*0.3 + 1.0; 
                                param[sp].w2mic[sm] = W2sym;
                                break;
                        case 3: /*ring of pairs across all taxa uses random w2mic*/
                                W2sym = -1.0;  
                                while (W2sym<error || W2sym>2.0) 
                                W2sym = random_normal()*0.3 + 1.0; 
                                param[sp].w2mic[sm] = W2sym;
                                break;
		}

            /***Give a fixed number of hosts a symbiont at location mic[sm] (# hosts: Nstart;  # to have symbiont: Nsymstart)***/
                if (Nsymstart > 0)
                {       i = 1;
                        while (i <= Nsymstart)
                        {
                /***Choose a host at random -- the host must previously have no symbiont of this taxon***/
                    label_redo:rn   = drand48();
                               add  = 0.0;
                               nowP = param[sp].firstP;
                               while (nowP != NULL)
                               {       add = add + 1.0/Nstart;
                                       if (rn < add)
                                       {       add = add - 1.0/Nstart;
                                               if (nowP->mic[sm] == 0) goto label_host;  /*continue:  taxon mic[sm]  absent in chosen host*/
                                               if (nowP->mic[sm] == 1) goto label_redo;  /*try again: taxon mic[sm] present in chosen host*/
                                               if (nowP->mic[sm] != 0 && nowP->mic[sm] != 1)
                                               { printf("Error in initialising microbiome -- exiting\n");   exit(1); }
                                        }
                                        nowP  = nowP->nextP;
		               }

                 /***Create a symbiont in chosen host at location mic[sm]***/
                     label_host:nowP-> mic[sm]  = 1;
                                nowP->w1mic[sm] = param[sp].w1mic[sm];
                                nowP->w2mic[sm] = param[sp].w2mic[sm];
                                i = i + 1;
                        }
                }
        }

    /***Call procedure to construct host fitness whost from its whole microbiome***/
        nowP = param[sp].firstP;
        while (nowP != NULL)
        {      host_fitness_microbiome(nowP);
               nowP  = nowP->nextP;
	}

/***INITIAL VALUES FOR BIRTH AND DEATH RATES***/

    /***Unlike the previous code, birth and death rates are computed for all hosts here;  if Nsymstart=0, MUST ensure whost=1***/
        nowP = param[sp].firstP;
        while (nowP != NULL)
        {
            /***Birth, death and encounter rates (death allows for microbiome; only females give birth!)***/
                switch(method_sym)
                {       case 1: printf("method_sym:1 not in use -- exiting\n");   exit(1);   break;
                        case 2: printf("method_sym:2 not in use -- exiting\n");   exit(1);   break;
                        case 3: printf("method_sym:3 not in use -- exiting\n");   exit(1);   break;
                        case 4: /*total host death rate becomes D/whost (applies to females and males*/
                                switch(nowP->sex)
                                {      case 0: /*female*/
                                               nowP->b = param[sp].B0;
                                               nowP->d = (param[sp].D0 + param[sp].DD * param[sp].n) / nowP->whost;
                                               nowP->e = param[sp].E0;
                                               break;
                                       case 1: /*male*/
                                               nowP->b = 0.0;
                                               nowP->d = (param[sp].D0 + param[sp].DD * param[sp].n) / nowP->whost;
                                               nowP->e = param[sp].E0;
                                               break;
                                }
                                nowP  = nowP->nextP;
                }
        }

    /***Output initial state of system***/
/*      if (flag_state == 1)    output_system_state();
*/
}


/****************************************************************************************************************
* Procedure for choosing probabilities of birth and death per unit time of each target host                     *
****************************************************************************************************************/

void    individual_event_probabilities()
{
        int     sp, i;
        repP    targP;         /* pointer to target individual*/

    /***Loop for setting probabilties p.u.t. for individuals***/
        for (sp=1; sp<nspp; ++sp)
        {       targP = param[sp].firstP;
                while (targP != NULL)
                {
                    /***Prob. of encounter p.u.t. (horizontal transfer) held fixed at param[sp].E0***/

                    /***Prob. of birth p.u.t. fixed at param[sp].B0 * whost;  targP->b entered when targP initialized***/

                    /***Prob. of death p.u.t. depends on n and on whether host is sym+ ***/
                        switch(method_sym)
                        {       case 1: printf("method_sym:1 not in use -- exiting\n");   exit(1);   break;
                                case 2: printf("method_sym:2 not in use -- exiting\n");   exit(1);   break;
                                case 3: printf("method_sym:3 not in use -- exiting\n");   exit(1);   break;
                                case 4: /*total host death rate becomes D/whost (applies to females and males)*/
                                        targP->d = (param[sp].D0 + param[sp].DD * param[sp].n) / targP->whost;
                                        break;
                        }

                    /***End simulation if a negative birth or death rate has been encountered***/
                        if (targP->b < 0.0)
                        {       printf("Negative birth rate:  program halted\n");
                                exit(1);
                        }
                        if (targP->d < 0.0)
                        {       printf("Negative death rate:  program halted\n");
                                exit(1);
                        }

                        targP = targP->nextP;
                }
        }
}


/****************************************************************************************************************
* Procedure for choosing births and deaths at random.  Returns:                                                 *
*   dt:    random time step from Gillespie algorithm	                                                           *
*   event: random choice of birth or death                                                                      *
*   pickP: pointer to random individual chosen for the event                                                    *
****************************************************************************************************************/

void    choose_event(dt, event, pickP)
double  *dt;
char    *event;
repP    *pickP;
{
        FILE    *out;
        int     sp, i;
        repP    nowP;
        double  Esum[nspp], Bsum[nspp], Dsum[nspp];
        double  E, B, D, sum, rn1, rn2, add;

    /***Prob. p.u.t. of symbiont encounter event set in initialize_parameters() (E0=0, when no horiz. transmission)***/
        for (sp=1; sp<nspp; ++sp)
        {       Esum[sp] = 0.0;
                nowP  = param[sp].firstP;
                while (nowP != NULL)
                {       Esum[sp] = Esum[sp] + nowP->e;
                        nowP  = nowP->nextP;
                }
        }

    /***Probability per unit time of birth event***/
        for (sp=1; sp<nspp; ++sp)
        {       Bsum[sp] = 0.0;
                nowP  = param[sp].firstP;
                while (nowP != NULL)
                {       Bsum[sp] = Bsum[sp] + nowP->b;
                        nowP  = nowP->nextP;
                }
        }

    /***Probability per unit time of death event***/
        for (sp=1; sp<nspp; ++sp)
        {       Dsum[sp] = 0.0;
                nowP  = param[sp].firstP;
                while (nowP != NULL)
                {      Dsum[sp] = Dsum[sp] + nowP->d;
                       nowP  = nowP->nextP;
                }
        }

    /***Sum of probabilities of all events***/
        E = 0.0;
        B = 0.0;
        D = 0.0;
        for (sp=1; sp<nspp; ++sp)
        {       E = E + Esum[sp];
                B = B + Bsum[sp];
                D = D + Dsum[sp];
        }
        sum = E + B + D;

    /***Random time to next event (exponential distribution of waiting times -- Gillespie)***/
        if (sum > error)        *dt = -log(1.0 - drand48())/sum;

    /***End simulation if sum of event probabilities = 0***/
        else
        {        printf("Sum of event probabilities = 0:	 program halted\n");
                exit(1);
        }

    /***Choose at random an encounter with symbiont, or birth or death event***/
        rn1 = drand48();
        if (rn1< E/sum)                 *event = 'e';
        if (rn1>=E/sum && rn1<(E+B)/sum)*event = 'b';
        if (rn1>=(E+B)/sum)             *event = 'd';

    /***Initialize variables for working out the individual to which the event happens***/
        rn2  = drand48();
        add  = 0.0;

    /***Choose at random the individual to which event happens***/
        switch (*event)
        {
                case 'e': /*symbiont encounter event; horizontal tranmission*/
                        for (sp=1; sp<nspp; ++sp)
                        {       nowP = param[sp].firstP;
                                while (nowP != NULL)
                                {       add = add + nowP->e/E;
                                        if (rn2 < add)
                                        {       add = add - nowP->e/E;
                                                goto label_e;
                                        }
                                        nowP  = nowP->nextP;
                                }
                        }
                label_e:break;

                case 'b': /*birth event; males have b=0 and cannot be chosen*/
                        for (sp=1; sp<nspp; ++sp)
                        {       nowP = param[sp].firstP;
                                while (nowP != NULL)
                                {       add = add + nowP->b/B;
                                        if (rn2 < add)
                                        {       add = add - nowP->b/B;
                                                goto label_b;
                                        }
                                        nowP  = nowP->nextP;
                                }
                        }
                label_b:break;

                case 'd': /*death event*/
                        for (sp=1; sp<nspp; ++sp)
                        {       nowP = param[sp].firstP;
                                while (nowP != NULL)
                                {       add = add + nowP->d/D;
                                        if (rn2 < add)
                                        {       add = add - nowP->d/D;
                                                goto label_d;
                                        }
                                        nowP  = nowP->nextP;
                                }
                        }
                label_d:break;
        }
        *pickP  = nowP;

    /***Record of calculations***/
        if (flag_choose == 1)
        {       out = fopen("output.detail.dat", "a");
                fprintf(out, "\nTime:%9.6f  State at end of choose_event procedure\n", t+*dt);
                fprintf(out, "E/sum:%6.3f  B/sum:%6.3f  D/sum:%6.3f  rn1:%6.4f	event:%c    ", 
                        E/sum, B/sum, D/sum, rn1, *event);
                fprintf(out, "rn2:%6.4f  add:%6.4f  pickP:%p  firstP:%p\n", rn2, add, *pickP, param[nowP->sp].firstP);
                fprintf(out, " pickP:%p  sex:%d  sym:%d  mod:%d\n", *pickP, nowP->sex, nowP->sym, nowP->mod);
                fflush (out);
                fclose (out);
        }
}


/****************************************************************************************************************
* Choose random father for birth event (needed if transmission is biparental)                                    *
****************************************************************************************************************/

void    choose_dad(sp, dadP)
int     sp;              /* species of mother (same sp for father)     */
repP    *dadP;	        /* pointer to the individual chosen as father	*/
{
        FILE    *out;
        repP    nowP;
        int     ndad;
        double  dadsum;
        double  rn, add;

    /***Number of potential fathers***/
        ndad   = 0;
        dadsum = 0.0;
        nowP   = param[sp].firstP;
        while (nowP != NULL)
        {       if (nowP->sex == 1)
                {       ndad   += 1;    
                        dadsum += 1.0;
                }
                nowP  = nowP->nextP;
        }
        if (ndad == 0)
        {       printf("no males in population -- program halted\n");
                exit(1);
        }

    /***Choose father at random***/
        rn   = drand48();
        add  = 0.0;
        nowP = param[sp].firstP;
        while (nowP != NULL)
        {       if (nowP->sex == 1)
                {       add = add + 1.0/dadsum;
                        if (rn < add)
                        {       add = add - 1.0/dadsum;
                                goto label_dad;
                        }
                }
		nowP  = nowP->nextP;
        }
        label_dad: *dadP = nowP;

    /***Record of calculations***/
        if (flag_choosedad == 1)
        {       out = fopen("output.detail.dat", "a");
                fprintf(out, "\nTime:%9.6f  State at end of choose_dad procedure\n", t);
                fprintf(out, "firstP:%11d rn:%6.4f  add:%6.4f  dadsum:%6.4f\n", param[nowP->sp].firstP, rn, add, dadsum);
                fprintf(out, "  dadP:%11d  sex:%d  mod:%2d  whost:%5.3f\n", *dadP, nowP->sex, nowP->mod, nowP->whost);
                fflush (out);
                fclose (out);
        }
}


/****************************************************************************************************************
* Construct microbiome in newborn host individual by vertical transmission (single and multiple symbiont taxa)  *
****************************************************************************************************************/
void    transmission_microbiome(pickP, dadP, birthP)
repP    pickP;          /* pointer to mother    */
repP    dadP;           /* pointer to father    */
repP    birthP;         /* pointer to offspring */
{
        FILE    *out;
        int     i, sm;
        double  rn, store;

    /***Symbiont status of host offspring depends on mode of transmission; now a whole microbiome***/
        switch(method_vertical)
        {
                case 0: /*maternal: symbionts all from mother*/
                        for (sm=1; sm<Nmic; ++sm)    birthP->  mic[sm] = pickP->  mic[sm];
                        for (sm=1; sm<Nmic; ++sm)    birthP->w1mic[sm] = pickP->w1mic[sm];
                        for (sm=1; sm<Nmic; ++sm)    birthP->w2mic[sm] = pickP->w2mic[sm];
                        break;

                case 1: /*paternal*/
                        printf("Paternal transmission not implemented: exiting from program\n");
                        exit(1);
                        break;

                case 2: /*biparental: host offspring gets the symbiont if one or both parents contain it*/
                        for (sm=1; sm<Nmic; ++sm)
                        {       if (pickP->mic[sm]==1 || dadP->mic[sm]==1)
                                {       
                                        birthP-> mic[sm] = 1;
                                        if (pickP->mic[sm]==1 && dadP->mic[sm]==0) {  birthP->w1mic[sm] = pickP->w1mic[sm];
                                                                                      birthP->w2mic[sm] = pickP->w2mic[sm];  }
                                        if (pickP->mic[sm]==0 && dadP->mic[sm]==1) {  birthP->w1mic[sm] =  dadP->w1mic[sm];
                                                                                      birthP->w2mic[sm] =  dadP->w2mic[sm];  }
                                        if (pickP->mic[sm]==1 && dadP->mic[sm]==1) {  birthP->w1mic[sm] = pickP->w1mic[sm];
                                                                                      birthP->w2mic[sm] = pickP->w2mic[sm];  }
                                }
                                else
                                {       birthP->  mic[sm] = 0;
                                        birthP->w1mic[sm] = 1.0;
                                        birthP->w2mic[sm] = 0.0;
                                }
                                /*Check to make sure both parents have the same symbiont at element sm? (wmic's should be equal)*/
			}
                        break;

                case 3: /*evolving: host offspring gets the symbiont if mother has it, OR if father has */
                        /*it and also has modfier gene mod:1; this allows transmission from father to evolve*/
                        for (sm=1; sm<Nmic; ++sm)
                        {       if (pickP->mic[sm]==1 || (dadP->mic[sm]==1 && dadP->mod==1))
                                {       birthP-> mic[sm] = 1;
                                        if (pickP->mic[sm]==1 && dadP->mic[sm]==0) {  birthP->w1mic[sm] = pickP->w1mic[sm];
                                                                                      birthP->w2mic[sm] = pickP->w2mic[sm];  }
                                        if (pickP->mic[sm]==0 && dadP->mic[sm]==1) {  birthP->w1mic[sm] =  dadP->w1mic[sm];
                                                                                      birthP->w2mic[sm] =  dadP->w2mic[sm];  }
                                        if (pickP->mic[sm]==1 && dadP->mic[sm]==1) {  birthP->w1mic[sm] = pickP->w1mic[sm];
                                                                                      birthP->w2mic[sm] = pickP->w2mic[sm];  }
                                }
                                else
                                {       birthP->  mic[sm] = 0;
                                        birthP->w1mic[sm] = 1.0;
                                        birthP->w2mic[sm] = 0.0;
                                }
                        }
                        break;
        }

    /***Call procedure to construct fitness whost of newborn host from its whole microbiome***/
        host_fitness_microbiome(birthP);

    /***Birth and death rate of newborn host individual computed in procedure birth()***/

}


/****************************************************************************************************************
* Birth of individual.                                                                                          *
*   pickP:  pointer to mother                                                                                   *
*   dadP:   pointer to father                                                                                   *
*   birthP: pointer to newborn individual                                                                       *
*                                                                                                               *
* Calls: choose_dad(), needed for symbiotic status of newborn individual, if transmission biparental            *
* Calls: symbiont_transmission() to determine whether newborn individual is sym+ or sym-                        *
*                                                                                                               *
* The new individual is placed immediately after the mother in the list;  its properties are defined in this    *
* procedure. This includes whether a newborn male will be able to transmit a symbiont to its offspring.         *
****************************************************************************************************************/

void    birth(pickP)
repP    pickP;         /* pointer to mother, already chosen to give birth */
{
        FILE    *out;
        int     sp, sm;
        double  rn;
        repP    birthP;       /* pointer to newborn individual            */
        repP    dadP;         /* pointer to individual chosen as father   */


/***CHOOSE FATHER OF NEWBORN HOST***/

    /***A random father from choose_dad(): needed for symbiont transmission and for evolution of modifier gene***/
    /***Must be done before newborn individual added to the list -- otherwise neonate could become its own father!***/
        choose_dad(pickP->sp, &dadP);


/***CREATE NEWBORN HOST***/

    /***Allocate memory for newborn individual***/
        birthP = NULL;
        birthP = (repP) malloc(sizeof(REP));
        if (!birthP)
        {       printf("Unable to allocate memory to pointer: exiting from program\n");
                exit(1);
        }

    /***WARNING!  It is essential to alter the pointers in the correct sequence below***/

    /***If the new individual is at the end of the list, the pointer to the next individual is set to NULL***/
    /*** and a record of the tail pointer is made***/
        if (pickP->nextP == NULL)
        {       birthP->nextP          = NULL;
                param[pickP->sp].lastP = birthP;
        }

    /***If the new individual is not at the end of the list, it is linked to the next individual ***/
        else
        {       birthP->nextP          = pickP->nextP;
                pickP->nextP->prevP    = birthP;
        }

    /***Lastly the links are made between the parent and the new individual***/
        birthP->prevP = pickP;
        pickP->nextP  = birthP;


/***PROPERTIES OF NEWBORN HOST***/

    /***Neonate: same species as mother***/
        birthP->sp = pickP->sp;

    /***Neonate: sex chosen at random with probability 0.5***/
        rn = drand48();
        if (rn <  0.5)  birthP->sex = 0;        /* female */
        if (rn >= 0.5)  birthP->sex = 1;        /* male   */

    /***Neonate: modifier gene for male transmission (only used in case 3; not needed for case 0 and 2)***/
        switch (method_vertical)
        {       case 0: /*maternal: male transmission always switched off*/
                        birthP->mod = 99;
                        break;
                case 1: /*paternal (not implemented)*/
                        break;
                case 2: /*biparental: male transmission always switched on*/
                        birthP->mod = 99;
                        break;
                case 3: /*evolution: modifier gene inherited from father (on Y chromosome)*/
                        /*daughters do not carry the modifier gene; dummy value inserted*/
                        if (birthP->sex==0)     birthP->mod = 99;
                        if (birthP->sex==1)     birthP->mod = dadP->mod;
                        break;
        }

    /***Neonate: microbiome of newborn host comes in the tramsmission process***/  
    /***Host fitness procedure called from within transmission_microbiome()***/ 
        transmission_microbiome(pickP, dadP, birthP);

    /***Neonate: birth, death and encounter rates (death allows for microbiome; only females give birth!)***/
        sp = birthP->sp;
        switch(method_sym)
        {       case 1: printf("method_sym:1 not in use -- exiting\n");   exit(1);   break;
                case 2: printf("method_sym:2 not in use -- exiting\n");   exit(1);   break;
                case 3: printf("method_sym:3 not in use -- exiting\n");   exit(1);   break;
                case 4: /*symbiont makes total host death rate D/whost (applies to females and males)*/
                        switch(birthP->sex)
                        {       case 0: /*female*/
                                        birthP->b = param[sp].B0;
                                        birthP->d = (param[sp].D0 + param[sp].DD * param[sp].n) / birthP->whost;
                                        birthP->e = param[birthP->sp].E0;
                                        break;
                                case 1: /*male*/
                                        birthP->b = 0.0;
                                        birthP->d = (param[sp].D0 + param[sp].DD * param[sp].n) / birthP->whost;
                                        birthP->e = param[birthP->sp].E0;
                                        break;
                        }
                        break;
        }

    /***The number of individuals of this species is increased by 1***/
        param[pickP->sp].n = param[pickP->sp].n + 1;

    /***Details of calculation***/
        if (flag_birth == 1)
        {       out = fopen("output.detail.dat", "a");
                fprintf(out, "\nTime:%9.6f  State at end of birth procedure\n", t);
                switch(method_vertical)
                {       case 0: fprintf(out, "Maternal transmission: microbiome from mother\n");                   break;
                        case 1: fprintf(out, "Paternal transmission: microbiome from father\n");                   break;
                        case 2: fprintf(out, "Biparental transmission: union of mother and father microbiomes\n"); break;
                        case 3: fprintf(out, "Transmission evolving\n"); break;
                }

                fprintf(out, " pickP:%11d  sex:%d  mod:%2d  whost:%5.3f  mic:", pickP, pickP->sex, pickP->mod, pickP->whost);
                for (sm=1; sm<Nmic; ++sm)  fprintf(out, "%d",     pickP->mic  [sm]);   fprintf(out, "  w1: ");
                for (sm=1; sm<Nmic; ++sm)  fprintf(out, "%4.2f ", pickP->w1mic[sm]);
                if (method_fitness==2 || method_fitness==3)
                {       fprintf(out, "  w2: ");
                        for (sm=1; sm<Nmic; ++sm)  fprintf(out, "%4.2f ", pickP->w2mic[sm]);
		}
                fprintf(out, "\n");

                fprintf(out, "  dadP:%11d  sex:%d  mod:%2d  whost:%5.3f  mic:",  dadP,  dadP->sex,  dadP->mod,  dadP->whost);
                for (sm=1; sm<Nmic; ++sm)  fprintf(out, "%d",      dadP->mic  [sm]);   fprintf(out, "  w1: ");
                for (sm=1; sm<Nmic; ++sm)  fprintf(out, "%4.2f ",  dadP->w1mic[sm]);
                if (method_fitness==2 || method_fitness==3)
                {       fprintf(out, "  w2: ");
                        for (sm=1; sm<Nmic; ++sm)  fprintf(out, "%4.2f ",  dadP->w2mic[sm]);
		}
                fprintf(out, "\n");

                fprintf(out, "birthP:%11d  sex:%d  mod:%2d  whost:%5.3f  mic:",birthP,birthP->sex,birthP->mod,birthP->whost);
                for (sm=1; sm<Nmic; ++sm)  fprintf(out, "%d",    birthP->mic  [sm]);   fprintf(out, "  w1: ");
                for (sm=1; sm<Nmic; ++sm)  fprintf(out, "%4.2f ",birthP->w1mic[sm]);
                if (method_fitness==2 || method_fitness==3)
                {       fprintf(out, "  w2: ");
                        for (sm=1; sm<Nmic; ++sm)  fprintf(out, "%4.2f ",birthP->w2mic[sm]);
		}
                fprintf(out, "\n");

                fflush (out);
                fclose (out);
        }
}


/****************************************************************************************************************
* Death of individual                                                                                           *
****************************************************************************************************************/

void    death(pickP)
repP    pickP;           /* pointer to the individual chosen for event */
{
        FILE    *out;

    /***If the dead individual is in the interior of the list, the the two adjoining individauls are linked***/
        if (pickP->nextP !=NULL && pickP->prevP !=NULL)
        {       pickP->prevP->nextP     = pickP->nextP;
                pickP->nextP->prevP     = pickP->prevP;
        }

    /***If the dead individual is at the head of the list, the next individual becomes the head***/
        if (pickP->prevP == NULL && pickP->nextP != NULL)
        {       pickP->nextP->prevP     = NULL;
                param[pickP->sp].firstP = pickP->nextP;
        }

    /***If the dead individual is at the tail of the list, the previous individual becomes the tail***/
        if (pickP->nextP == NULL && pickP->prevP != NULL)
        {       pickP->prevP->nextP     = NULL;
                param[pickP->sp].lastP  = pickP->prevP;
        }

    /***If the dead individual is at the head AND the tail of the list, the species is extinct***/
        if (pickP->nextP == NULL && pickP->prevP == NULL)
        {       param[pickP->sp].firstP = NULL;
                param[pickP->sp].lastP  = NULL;
        }

    /***The number of individuals of this species is decreased by 1***/
	param[pickP->sp].n = param[pickP->sp].n - 1;

    /***Record of calculations***/
        if (flag_death == 1)
        {       out = fopen("output.detail.dat", "a");
                fprintf(out, "\nTime:%9.6f  State at end of death procedure\n", t);
                fprintf(out, "  pickP: %p  pickP->sex:%d  pickP->sym:%d  pickP->mod:%d	pickP->whost:%6.4f\n", 
                        pickP,  pickP->sex,  pickP->sym,  pickP->mod,  pickP->whost);
                fflush (out);
                fclose (out);
        }
}


/****************************************************************************************************************
* Symbiont horizontal transmission.                                                                             *
*                                                                                                               *
* This is simple: it just changes the host from sym- to sym+, if the host was previously sym-.                  *
* No change to the number of hosts in the population.                                                           *
****************************************************************************************************************/

void    horiz_trans(pickP)
repP    pickP;               /* pointer to host encountering a symbiont */
{
        FILE    *out;
        int     sp, sm;

        sp = 1;
        sm = 1;

    /***Host encountering a symbiont in environment takes it up, changing its death rate***/
        switch (method_microbiome)
        {       case 1: /*single symbiont*/
                        if(pickP->mic[sm] == 0)
                        {       pickP->mic[sm]   = 1;
                                pickP->w1mic[sm] = Wsym;
                                pickP->whost     = Wsym;
                                pickP->d         = (param[sp].D0 + param[sp].DD * param[sp].n) / pickP->whost;
                        }
                        break;
                case 2: /*multiple symbiont taxa: rules for horizonal transmission notyet decided*/
                        printf("Rules for horiz trans of multiple symbiont taxa as yet undecided -- exiting\n");  exit(1);
                        break;
		}

    /***Details of calculation***/
        if (flag_horizontal == 1)
        {       out = fopen("output.detail.dat", "a");
                fprintf(out, "\nTime:%9.6f  State at end of horiz_trans procedure\n", t);
                fprintf(out, " pickP:%11d  sex:%d  mod:%2d  whost:%5.3f  mic:", pickP, pickP->sex, pickP->mod, pickP->whost);
                for (sm=1; sm<Nmic; ++sm)  fprintf(out, "%d",     pickP->mic  [sm]);   fprintf(out, "  w1: ");
                for (sm=1; sm<Nmic; ++sm)  fprintf(out, "%4.2f ", pickP->w1mic[sm]);
                if (method_fitness==2 || method_fitness==3)
                {       fprintf(out, "  w2: ");
                        for (sm=1; sm<Nmic; ++sm)  fprintf(out, "%4.2f ", pickP->w2mic[sm]);
		}
                fprintf(out, "\n");
                fflush (out);
                fclose (out);
        }
}


/****************************************************************************************************************
* Calculates various measures at each time, and puts them into output timeseries file:                          *
* Single symbiont version                                                                                       * 
****************************************************************************************************************/

void    output_timeseries_single(dt, event, pickP)
double  dt;               /* time increment to the event                 */
char    event;            /* event chosen in each iteration: b, d        */
repP    pickP;            /* pointer to the individual chosen for event  */
{
        FILE    *out;
        int     sp, sm;
        repP    nowP;     /* pointer to individuals in list              */

        sp = 1;
        sm = 1;

    /***Breakdown to numbers by sex and symbiont status (just for checking single timeseries)***/
        param[sp].nfsym1 = 0;
        param[sp].nmsym1 = 0;
        param[sp].nfsym0 = 0;
        param[sp].nmsym0 = 0;

        nowP = param[sp].firstP;
        while (nowP != NULL)
        {       if (nowP->sex==0 && nowP->mic[sm]==1)       param[sp].nfsym1 += 1;
                if (nowP->sex==1 && nowP->mic[sm]==1)       param[sp].nmsym1 += 1;
                if (nowP->sex==0 && nowP->mic[sm]==0)       param[sp].nfsym0 += 1;
                if (nowP->sex==1 && nowP->mic[sm]==0)       param[sp].nmsym0 += 1;
                nowP = nowP->nextP;
        }

    /***Check sex ratio (just for output)***/
        param[sp].nf = 0;
        param[sp].nm = 0;
        nowP = param[sp].firstP;
        while (nowP != NULL)
        {       switch(nowP->sex)
                {       case 0: /*female*/
                                param[sp].nf += 1;
                                break;
                        case 1: /*male*/
                                param[sp].nm += 1;
                                break;
                }
                nowP = nowP->nextP;
        }

    /***Frequency of gene switching off male transmission (maternal)***/
        param[sp].nmat = 0;
        param[sp].nbip = 0;
        nowP = param[sp].firstP;
        while (nowP != NULL)
        {       if (nowP->sex==1)
                {       if(nowP->mod==0)    param[sp].nmat += 1;
                        else                param[sp].nbip += 1;
                }
                nowP = nowP->nextP;
        }
        param[sp].fmat = (double)param[sp].nmat / (double)(param[sp].nmat + param[sp].nbip);

    /***Association of sym status and transmission type, including coefficient of disequilibrium***/
        param[sp].nsym0mat = 0;
        param[sp].nsym1mat = 0;
        param[sp].nsym0bip = 0;
        param[sp].nsym1bip = 0;
        nowP = param[sp].firstP;
        while (nowP != NULL)
        {       if (nowP->sex==1)
                {       if(nowP->mic[sm]==0 && nowP->mod==0)     param[sp].nsym0mat += 1;
                        if(nowP->mic[sm]==1 && nowP->mod==0)     param[sp].nsym1mat += 1;
                        if(nowP->mic[sm]==0 && nowP->mod==1)     param[sp].nsym0bip += 1;
                        if(nowP->mic[sm]==1 && nowP->mod==1)     param[sp].nsym1bip += 1;
                }
                nowP = nowP->nextP;
        }
        param[sp].dis = (double)(param[sp].nsym0mat * param[sp].nsym1bip - param[sp].nsym1mat *param[sp].nsym0bip)
                      / (double)(param[sp].nm * param[sp].nm);

    /***Record results on screen***/
        printf("t:%9.6f sp:%d   %6d  %c  ", t, pickP->sp, nevent, event);
        printf("%4d  %4d:%4d  %6.4f    %6.4f  ", 
               param[sp].n, param[sp].nf, param[sp].nm, param[sp].fmic[sm], param[sp].fmat);
        printf("\n");

    /***Record results to file***/
        out = fopen("output.timeseries.dat", "a");
        if (t-dt < error)
        {       fprintf(out, "#time     sp  nevent event   n    nf    nm  fsymbiont  nmat  nbip fmaternal   ");
                fprintf(out, "mat- mat+ bip- bip+      diseq\n");
	}
        fprintf(out, "%9.6f  %d  %6d   %c  ", t, pickP->sp, nevent, event);
        fprintf(out, "%4d  %4d  %4d  %6.4f    ", param[sp].n,	 param[sp].nf, param[sp].nm, param[sp].fmic[sm]);
        fprintf(out, " %4d  %4d  %6.4f    ", param[sp].nmat, param[sp].nbip, param[sp].fmat);
        fprintf(out, " %4d %4d %4d %4d    ", param[sp].nsym0mat, param[sp].nsym1mat, param[sp].nsym0bip, param[sp].nsym1bip);
        fprintf(out, " %7.4f  ", param[sp].dis);
        fprintf(out, "\n");
        fflush (out);
        fclose (out);
}


/****************************************************************************************************************
* Calculates various measures at each time, and puts them into output timeseries file:                          *
* Multiple symbiont taxa -- microbiome version                                                                  *
****************************************************************************************************************/

void    output_timeseries_microbiome(dt, event, pickP)
double  dt;               /* time increment to the event                 */
char    event;            /* event chosen in each iteration: b, d        */
repP    pickP;            /* pointer to the individual chosen for event  */
{
        FILE    *out;
        int     sp, sm;
        double  sum;
        repP    nowP;     /* pointer to individuals in list              */

        sp = 1;

    /***Record results on screen***/
        printf("t:%9.6f sp:%d   %6d  %c  ", t, pickP->sp, nevent, event);
        printf("%4d  %4d:%4d  %6.4f    %6.4f  ", 
               param[sp].n, param[sp].nf, param[sp].nm, param[sp].fsym, param[sp].fmat);
/*      printf("  address:%d  ", pickP);
*/      printf("\n");

    /***Record results to file***/
        out = fopen("output.timeseries.dat", "a");
        if (t-dt < error)
        {       fprintf(out, "#time     sp  nevent event n w1:");
                for(sm=1; sm<Nmic; ++sm)    fprintf(out, "%5.2f ",   param[sp].w1mic[sm]);    fprintf(out, "\n");
                fprintf(out, "#                            w2:");
                for(sm=1; sm<Nmic; ++sm)    fprintf(out, "%5.2f ",   param[sp].w2mic[sm]);    fprintf(out, "\n");
                fprintf(out, "\n");
	}

        fprintf(out, "%9.6f  %d  %6d   %c  %5d ", t, pickP->sp, nevent, event, param[sp].n);
        for(sm=1; sm<Nmic; ++sm)    fprintf(out, "%5.3f ", param[sp].fmic[sm]);    fprintf(out, "   ");
        for(sm=1; sm<Nmic; ++sm)    fprintf(out, "%5d ",   param[sp].nmic[sm]);    fprintf(out, "\n");
        fflush (out);
        fclose (out);

    /***Display the microbiomes of all hosts at the end (recorded in output.detail.dat)***/
    /***They ore ordered because offspring inserted immediately after the mother***/ 
/*        if (t > tmax)    output_system_state();
*/
}


/****************************************************************************************************************
* Runs a realization of the stochastic process. At each time step:                                              *
* Assigns to each individual a probability p.u.t. of a birth, death or horizontal transmission event in:        *
*     individual_event_probabilities()                                                                          *
* Chooses a random time, event and individual in:                                                               *
*     choose_event()                                                                                            *
* Updates the state of the system, calling one of:                                                              *
*     birth()                                                                                                   *
*     death()                                                                                                   *
*     horiz_trans()                                                                                             *
****************************************************************************************************************/

void    simulate()
{
        FILE    *out;
        int     i, j, sp, sm;

        double  dt;             /* time increment to the event	                   */
        char    event;          /* event chosen in each iteration: b, d, e       */
        repP    pickP;          /* pointer to the individual chosen for event    */
        repP    nowP;

        sp = 1;

    /***Initialise hosts and their microbiomes: includes singletons, and mutliple symbiont taxa***/
        initialize_hosts_microbiome();

    /***Loop through births and deaths until sym+ hosts reach a frequency ~0.0 or ~1.0 (tmax is a backstop)***/
        t = 0.0;
        nevent = 0;

    /***Use one of the loops below; comment out the other one***/

    /***Loop to run for a fixed time tmax -- needed for displaying a single time series***/
        while (t<=tmax)

    /***Loop to run until loss or fixation of the symbiont***/
/*      param[sp].fsym = 0.1;
        if (method_microbiome == 2) { printf("time loop not checked for microbiome -- exiting\n"); exit(1); }
        while (t<=tmax && param[sp].fsym>error && param[sp].fsym<1.0-error && param[sp].n>1)
*/      {
            /***Output state of system***/
                if (flag_state == 1)    output_system_state();

            /***Calculate the birth and death probabilities per unit time of each individual***/
                individual_event_probabilities();

            /***Choose a birth or death event at random***/
                choose_event(&dt, &event, &pickP);

            /***Update the time and event number***/
                t = t + dt;
                nevent += 1;

            /***Update the state of the system***/
                switch  (event)
                {       case 'e':       horiz_trans(pickP);  break;
                        case 'b':       birth(pickP);        break;
                        case 'd':       death(pickP);        break;
                }

            /***Number and frequency of hosts with symbiont (needed in main() below)***/
                for(sm=1; sm<Nmic; ++sm)   param[sp].nmic[sm] = 0;
                nowP = param[sp].firstP;
                while (nowP != NULL)
                {       for(sm=1; sm<Nmic; ++sm)   if (nowP->mic[sm]==1)    param[sp].nmic[sm] += 1;
                        nowP = nowP->nextP;
                }
                for(sm=1; sm<Nmic; ++sm)   param[sp].fmic[sm] = (double)param[sp].nmic[sm] / (double)param[sp].n;

            /***Output to file output.timeseries.dat if required (calculates various extra measures along the way)***/
                if (flag_timeseries == 1)
                switch (method_microbiome)
                {       case 1: output_timeseries_single    (dt, event, pickP);  break;
                        case 2: output_timeseries_microbiome(dt, event, pickP);  break;
		}

            /***If event was death, free the memory previously given to the individual***/
                if (event == 'd')       free(pickP);
        }
}


/****************************************************************************************************************
* Procedure to accumulate frequency distribution of Wsym's that allowing invasion of symbiont                   *
****************************************************************************************************************/

void    fitness_freq(iter)
int     iter;
{
        int     sp, bin;
        double  binmin, binmax;

    /***At start set the bins***/
        if (iter == 0)
        {       for (bin=0; bin<nbin; ++bin)
                {       freq[0][bin] = bin * binwidth;
                        freq[1][bin] = 0.0;
                }
        }

    /***Thereafter, accumulate frequencies of Wsym***/
        else
        {       sp  = 1;
                bin = (int) (Wsym / binwidth);
                freq[1][bin] += 1.0;
                printf("Wsym:%6.4f   freq[%2d]:%4.0f  ", Wsym, bin, freq[1][bin]);
        }
}


/****************************************************************************************************************
* Main body of program                                                                                          *
*   Loops through a specified number of realizations (niterate) of the stochastic process.                      *
*   Calls procedure simulate(), which runs the realization.                                                     *
*   Checks whether the realization meets a condition for successful invasion by the symbiont.                   *
*   If successful, records in output.invasions.dat, and adds w into a frequency distribution in fitness_freq()  *
****************************************************************************************************************/

int     main()
{
        FILE    *out;
        int     iter, bin, sp, count_fsuccess, count_fixed, count_lost, count_undecided;

    /***Initialize output files (ensures they are empty at the start)***/
        fclose(fopen("output.detail.dat",              "w"));
        fclose(fopen("output.timeseries.dat",          "w"));
        fclose(fopen("output.invasions.dat",           "w"));
        fclose(fopen("output.invasions.fixed.dat",     "w"));
        fclose(fopen("output.invasions.lost.dat",      "w"));
        fclose(fopen("output.invasions.undecided.dat", "w"));
        fclose(fopen("output.freq.distrib.dat",        "w"));

    /***Seed the random number generator***/
/*        printf("\nSeed for random number generator (up to 6 digits):      ");
        scanf("%6d", &seed);
*/
seed = 7597;
        srand48(seed);
        out = fopen("output.detail.dat", "w");
        fprintf(out, "Seed for random number generator:%d\n\n", seed);
        fflush (out);
        fclose (out);

    /***Set the birth and death rate parameters***/
        initialize_parameters();

    /***Set header line for output file for list of symbiont invasions***/
        out = fopen("output.invasions.dat", "a");
        fprintf(out, "# iter       Wsym      n   nsym     fsym        t   nevent\n");
        fflush (out);
        fclose (out);
        out = fopen("output.invasions.fixed.dat", "a");
        fprintf(out, "# iter  c_fixed       Wsym      n   nsym     fsym       t   nevent\n");
        fflush (out);
        fclose (out);
        out = fopen("output.invasions.lost.dat", "a");
        fprintf(out, "# iter   c_lost       Wsym      n   nsym     fsym       t   nevent\n");
        fflush (out);
        fclose (out);
        out = fopen("output.invasions.undecided.dat", "a");
        fprintf(out, "# iter c_undecided    Wsym      n   nsym     fsym       t   nevent\n");
        fflush (out);
        fclose (out);

    /***A check to ensure method_micobiome:1 for Nmic:2 (single symbiont), and method_micobiome:2 for Nmic:>2***/
        if ((method_microbiome==1 && Nmic>2) || (method_microbiome==2 && Nmic==2)) 
        {       printf("method_microbiome:%d is incompatible with Nmic:%d -- exiting\n", method_microbiome, Nmic);
                exit(1);
	}

    /***A check to ensure method_fitness:1 for Nmic:2 (no pairwise interactions among symbiont taxa if only 1 taxon present)***/
        if ((method_fitness==2 && Nmic==2) || (method_fitness==3 && Nmic==2)) 
        {       printf("Pairwise interactions not possible with a single symbiont taxon -- exiting\n");
                exit(1);
	}

    /***A check to avoid running > 1 iteration when examining time series***/
        if (flag_timeseries == 1 && niterate > 1)       
        {       printf("Set niterate = 1 to examine a single timeseries: exiting from program\n");
                exit(1);
        }

    /***Set up output file for bin_invasion()***/
        iter = 0;
        fitness_freq(iter);

    /***Loop for iterating stochastic process***/
        sp = 1;
        count_fsuccess  = 0;
        count_fixed     = 0;
        count_lost      = 0;
        count_undecided = 0;
        for (iter=1; iter<=niterate; ++iter)
        {
                if (Nmic > 2 && niterate > 1)
                {       printf("Multiple iterations implemented only for a single symbiont -- exiting\n");
		        exit(1);
		}

            /***Choose an effect of symbiont on host fitness***/

            /***A single fixed value for one iteration***/
/*                if (niterate == 1)  Wsym = 0.37;
*/              if (niterate == 1)  Wsym = 0.4;

            /***Random values for multiple iterations***/
                else
                {       Wsym = -1.0;
                        while (Wsym < error || Wsym > 2.5)
                        Wsym = random_normal()*0.30 + 1.0;
		}

            /***Realization of stochastic birth-death process***/
                simulate();

                printf("\niteration:%3d   fsuccess:%3d   fixed:%3d   lost:%3d   undecided:%3d   ",
                       iter, count_fsuccess, count_fixed, count_lost, count_undecided);
                printf("\n");

            /***Add to list of symbiont invasions satisfying condtion: fsuccess***/
                if (param[sp].fsym >= fsuccess)
                {       count_fsuccess += 1;
                        out = fopen("output.invasions.dat", "a");
                        fprintf(out, "%6d   %8.6f   %4d   %4d   %6.4f   %6.2f  %6d", 
                                iter, Wsym, param[sp].n, param[sp].nsym, param[sp].fsym, t, nevent);
                        fprintf(out, "\n");
                        fflush (out);
                        fclose (out);

                    /***Accumulate frequency distribution of Wsym values allowing symbiont invasion***/
                    /***Half-bin shift for gnuplot box plot***/
                        fitness_freq(iter);
                        out = fopen("output.freq.distrib.dat", "w");
                        for (bin=0; bin<nbin; ++bin)
                        {       fprintf(out, "%6.4f   %6.4f   %4.0f", freq[0][bin], freq[0][bin]+binwidth/2.0, freq[1][bin]);
                                fprintf(out, "\n");
                        }
                        fflush (out);
                        fclose (out);
                }

            /***Add to list of iterations in which symbionts: fixed***/
                if (param[sp].fsym >= 1.0-error)
                {       count_fixed += 1;
                        out = fopen("output.invasions.fixed.dat", "a");
                        fprintf(out, "%6d   %6d   %8.6f   %4d   %4d   %6.4f   %6.2f  %6d", 
                                iter, count_fixed, Wsym, param[sp].n, param[sp].nsym, param[sp].fsym, t, nevent);
                        fprintf(out, "\n");
                        fflush (out);
                        fclose (out);
                }

            /***Add to list of iterations in which symbionts: lost***/
                if (param[sp].fsym <= error)
                {       count_lost += 1;
                        out = fopen("output.invasions.lost.dat", "a");
                        fprintf(out, "%6d   %6d   %8.6f   %4d   %4d   %6.4f   %6.2f  %6d", 
                                iter, count_lost, Wsym, param[sp].n, param[sp].nsym, param[sp].fsym, t, nevent);
                        fprintf(out, "\n");
                        fflush (out);
                        fclose (out);
                }

            /***Add to list of iterations in which symbionts: fate undecided***/
                if (t >= tmax)
                {       count_undecided += 1;
                        out = fopen("output.invasions.undecided.dat", "a");
                        fprintf(out, "%6d   %6d   %8.6f   %4d   %4d   %6.4f   %6.2f  %6d", 
                                iter, count_undecided, Wsym, param[sp].n, param[sp].nsym, param[sp].fsym, t, nevent);
                        fprintf(out, "\n");
                        fflush (out);
                        fclose (out);
                }
        }

        out = fopen("output.detail.dat", "a");
        fprintf(out, "\nNumber of invasions fsuccess:  %6d",   count_fsuccess);
        fprintf(out, "\nNumber of invasions fixed:     %6d",   count_fixed);
        fprintf(out, "\nNumber of invasions lost:      %6d",   count_lost);
        fprintf(out, "\nNumber of invasions undecided: %6d\n", count_undecided);
        fflush (out);
        fclose (out);
}



