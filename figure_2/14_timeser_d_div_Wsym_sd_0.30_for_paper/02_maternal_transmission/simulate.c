/****************************************************************************************************************
* Host stochastic birth-death process with vertical transmission of symbionts					*
*														*
* To remove segmentation faults, had to change type float to type double.  Think that's because the function	*
* drand48 is type double, so the variable rn also needed to be type double.											*
*														*
* 25.11.21.													*
****************************************************************************************************************/

#include	<stdlib.h>
#include	<math.h>
#include	<stdio.h>
#include	<string.h>


/****************************************************************************************************************
* Constants													*
****************************************************************************************************************/

#define		error	1e-9		/* rounding error allowed for type double			*/
#define		pi	3.1415927

#define		method_transmit	0	/* transmission of symbionts from host parents to host offspring
					   0: maternal
					   1: paternal    (not implemented)
					   2: biparental: transmitted if one or both parents have it
                                           3: evolving:	  modifier gene on Y chromosome: 0 mat; 1 bip	*/

#define		method_sym	4	/* symbiont effect on host
					   1: wsym acts as a factor on intrinsic birth rate: b0*w   
					   2: wsym acts as a factor on intrinsic death rate: d0*w
					   3: wsym acts as a factor on intrinsic death rate: d0/w
					   4: wsym acts as a factor on total     death rate: d/w	*/

#define		nspp		2	/* number of host species + 1					*/
#define		tmax		20.0	/* maximum time for simulations					*/
#define		Nstart		1000	/* initial population size of host species (must be > Nsymstart)*/
#define         Nsymstart       100     /* initial number of sym+ hosts (use: 0 to check a sym- popn)	*/

#define		fmodstart	0.5     /* Starting frequency for gene stopping male transmission:
					   fmodstart ~ 0.0: popn starts close to biparental transmission
					   fmodstart ~ 1.0: popn starts close to maternal transmission.
					   Modifier gene assumed to be on Y chromosome, passed on to and 
					   expressed only in sons: carried only in males, not in females*/

#define		niterate	1	/* number of iterations of stochastic process			*/
#define		nbin		50	/* size of array with Wsym frequency distribution		*/
#define		binwidth	0.05	/* binwidth for Wsym frequency distribution			*/
#define		fsuccess	0.1	/* threshold freq of  sym+ hosts to count as successful invasion;	
					   distinct from the stopping condition in simulate()		*/


/****************************************************************************************************************
* Global types													*
****************************************************************************************************************/

/***A structure to hold the basic information on a single individual one host species at a time***/
/***Designed for a bidirectional list in each host population***/
struct	rep
	{	int	sp;			/* host species						*/
		int	sex;			/* host sex                   0: female;  1: male	*/
		int	sym;			/* host symbiont              0: absent; 1: present	*/
                int     mod;                    /* modifier male transmission 0: off(mat);  1: on(bip)	*/
		double	wsym;			/* fitness of host allowing for symbiont		*/
		double	b;			/* probability of birth	per unit time			*/
		double	d;			/* probability of death	per unit time			*/
                struct	rep	*nextP, *prevP;	/* pointers for bidirectional list			*/
	};

typedef struct	rep	REP;			/*  REP	 is a new data type for declaring structure of type rep */
typedef		REP	*repP;			/* *repP is pointer to the data type REP.  (C book p613)	*/

/***A structure to hold information on properties of a host species***/
struct	species_properties
	{	int	n;			/* number of hosts (total)		*/
		int	nf;			/* number of hosts (female)		*/
		int	nm;			/* number of hosts (male)		*/
		int	nsym;			/* number of hosts with symbiont	*/
		double	fsym;			/* frequency of hosts with symbiont	*/

		int	nfsym0;			/* number of hosts (female sym-)	*/
		int	nmsym0;			/* number of hosts (male   sym-)	*/
		int	nfsym1;			/* number of hosts (female sym+)	*/
		int	nmsym1;			/* number of hosts (male   sym+)	*/

		int	nmat;			/* number of transmit- males (mat)	*/
		int	nbip;			/* number of transmit+ males (bip)	*/
                double  fmat;                   /* freq   of transmit- gene		*/

		int	nsym0mat;		/* number of sym- transmit- males	*/
		int	nsym0bip;		/* number of sym- transmit+ males	*/
		int	nsym1mat;		/* number of sym+ transmit- males	*/
		int	nsym1bip;		/* number of sym+ transmit+ males	*/
		double	dis;			/* coefficient of disequilibrium	*/

		double	B0;			/* intrinsic prob of birth p.u.t.	*/
		double	D0;			/* intrinsic prob of death p.u.t.	*/
		double	DD;			/* DD component to death p.u.t.		*/
						/* (no L^2: using # of individuals)	*/
		double	BD;			/* total difference: B-D		*/
		repP	firstP;			/* pointer to first individual		*/
		repP	lastP;			/* pointer to last  individual		*/
	};


/***Holds a pointer to the first element of lists of each species***/
struct	rep			state[nspp];

/***Holds properties of the species***/
struct	species_properties	param[nspp];


/****************************************************************************************************************
* Global variables                                                                                              *
****************************************************************************************************************/

double	t;					/* time					*/
int	nevent;					/* counter for number of events		*/

/***Variables for random number generator***/
double  drand48();				/* Function for random number generator */
int     seed;					/* Seed for random number generator     */

double	Wsym;					/* effect of symbiont on host fitness; implemented as  a factor 
						   by which intrinsic rate is multiplied
						     < 1:  host rate decreases
						     = 1:  host rate unchanged
						     > 1:  host rate increases					*/
double	freq[2][nbin];				/* array with frequency distribution of Wsym allowing invasion	*/

/***Flags to switch on detailed output for checking code***/
int	flag_state      = 0;			/* 1: detailed output for system state on;             0: off	*/
int	flag_choose     = 0;			/* 1: detailed output for choose event calculations on 0: off	*/
int	flag_choosedad  = 0;			/* 1: detailed output for choose event calculations on 0: off	*/
int	flag_transmit   = 0;			/* 1: detailed output for symbiont transmission on     0: off	*/
int	flag_birth      = 0;			/* 1: detailed output for birth    calculations on     0: off	*/
int	flag_death      = 0;			/* 1: detailed output for death    calculations on     0: off	*/
int	flag_timeseries = 1;			/* 1: detailed output for timeseries of realization on 0: off	*/


/****************************************************************************************************************
*  Function for generating random numbers on a normal distribution on [0,1]					*
****************************************************************************************************************/
double random_normal()
{
        int i;
        double theta, y, value1, value2;

        y     = drand48();
        y     = -2.0 * log(y);
	theta = drand48();
        theta = theta * 2.0 * pi;

        value1 = sqrt(y) * cos(theta);	/*only this one being used*/
        value2 = sqrt(y) * sin(theta);

        return(value1);
}


/****************************************************************************************************************
* Input parameter values											*
****************************************************************************************************************/

void initialize_parameters()
{
	FILE	*out;
	int	sp;

    /***Enter host species population parameters***/
	sp = 1;
	param[sp].n  = Nstart;
	param[sp].B0 = 4.0;  
	param[sp].D0 = 1.0;
	param[sp].DD = 0.001;

    /***Record of parameter values***/
	out = fopen("output.detail.dat", "a");
	fprintf(out, "HOST POPULATION PARAMETERS\n");
	sp = 1;
	fprintf(out, "  species: %d\n",         sp);
	fprintf(out, "  param[sp].n:    %5d\n", param[sp].n);
	fprintf(out, "  param[sp].B0: %7.3f\n", param[sp].B0);
	fprintf(out, "  param[sp].D0: %7.3f\n", param[sp].D0);
	fprintf(out, "  param[sp].DD: %7.3f\n", param[sp].DD);
	fprintf(out, "  Nhat = (B0*0.5 - D0) / DD  (0.5 because only females give birth)\n");

	fprintf(out, "\nREALISATIONS\n");
	fprintf(out, "  tmax:            %4.0f    time over which stochastic realisation runs\n",  tmax);
	fprintf(out, "  nspp-1:          %4d    number of host species\n",                       nspp-1);
	fprintf(out, "  Nstart:          %4d    initial population size of host spp\n",          Nstart);

	fprintf(out, "\nSYMBIONTS\n");
	fprintf(out, "  Nsymstart:       %4d    initial   number of hosts with symbiont\n",        Nsymstart);
	fprintf(out, "  fsuccess:        %4.2f    threshold freq for successful invasion\n",        fsuccess);
	fprintf(out, "  method_transmit: %4d    symbionts: 0 maternal; 2 biparental\n",      method_transmit);
	fprintf(out, "  method_sym:      %4d    symbionts affect host: 1 B0*w; 2 D0*w; 3 D0/w; 4 D/w\n", method_sym);

	if(method_transmit==3)
	{
	fprintf(out, "\nMALE TRANSMISSION MODIFIER\n");
	fprintf(out, "	fmodstart:     %6.4f	initial frequency of no-male-transmit modifier\n",	fmodstart);
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
* This procedure outputs the current state of the system							*
****************************************************************************************************************/

void	output_system_state()
{
	FILE	*out;
	int	sp, i, j;
	repP    nowP;

    /***Output bidirectional list of host individuals by species***/
	out = fopen("output.detail.dat", "a");
	fprintf(out, "\n");
	for (sp=1; sp<nspp; ++sp)
	{	nowP = param[sp].firstP;
		fprintf(out, "Time:%9.6f nevent:%d sp:%d   n:%d   nf:%d   nm:%d\n", t, nevent, sp, param[sp].n, param[sp].nf, param[sp].nm);
		fprintf(out, "Species:%d   n:%d   nf:%d   nm:%d\n", sp, param[sp].n, param[sp].nf, param[sp].nm);
		fprintf(out, "   i       prevP        nowP       nextP  sex     birth  death   sym  wsym     mod\n");
		i = 1;
		while (nowP != NULL)
		{	fprintf(out, "%4d %11d %11d %11d    %d    %6.3f %6.3f     %d %6.3f     %2d  \n", 
				i, nowP->prevP, nowP, nowP->nextP, nowP->sex, nowP->b, nowP->d, nowP->sym, nowP->wsym, nowP->mod);
			i = i + 1;
			nowP = nowP->nextP;
		}
	}
	fflush(out);
	fclose(out);
}


/****************************************************************************************************************
* This procedure initializes properties of individual hosts at start of simulate				*
****************************************************************************************************************/

void	initialize_hosts()
{
	int	sp, i, j;
	double	rn;
	repP	firstP, lastP, nowP, prevP;

	sp = 1;

    /***Set the size of host population***/
	param[sp].n = Nstart;
	if (param[sp].n < 1)
	{	printf("n:%d  Must have individuals in the population -- program halted\n", param[sp].n);
		exit(1);
	}

    /***Set the number of hosts carrying symbiont***/
	if (Nstart < Nsymstart)
	{	printf("Nstart:%d  Nsymstart:%d not a good starting point -- program halted\n", Nstart, Nsymstart);
		exit(1);
	}

    /***Set the first and last pointers to null***/
	param[sp].firstP = NULL;
	param[sp].lastP	 = NULL;

    /***Set the start pointer as long as individuals actually exist***/
	param[sp].firstP = &state[sp];
	prevP = NULL;
	nowP  = param[sp].firstP;

    /***Set up bidirectional list for all individuals by means of pointers***/
    /***An extra individual '0' is created and removed below, so all memory allocation is by malloc***/
	for (i=0;   i<=param[sp].n;  ++i)
	{
	    /***Allocate the species of the individual***/
		nowP->sp = sp;

	    /***Sex of the individual***/
		rn = drand48();
		if(rn <  0.5) nowP->sex = 0; /* female */
		if(rn >= 0.5) nowP->sex = 1; /* male   */

	    /***Initial values for birth and death rates (only females give birth!)***/
		switch(nowP->sex)
		{	case 0: /*female*/
				nowP->b = param[sp].B0;
				nowP->d = param[sp].D0 + param[sp].DD * param[sp].n;
				break;
			case 1: /*male*/
				nowP->b = 0.0;
				nowP->d = param[sp].D0 + param[sp].DD * param[sp].n;
				break;
		}

	    /***Set status to symbiont absent and wsym=1 (can be altered below)***/
		nowP->sym  = 0;
		nowP->wsym = 1.0;
				
	    /***For evolution of transmission, define modifier gene of males (nonsense value for females as gene on Y chromosome)***/
		switch (method_transmit)
		{	case 0: /*maternal: male transmission always off*/
				nowP->mod = 0;
				break;
			case 1: /*paternal (not implemented)*/
				nowP->mod = 0;
				break;
			case 2: /*biparental: male transmission always on*/
				nowP->mod = 1;
				break;
			case 3: /*modifier gene for  male transmission: 0:off maternal, 1:on biparental*/
				/*mod not defined for females; dummy value inserted*/
				if (nowP->sex==0)		nowP->mod = 99;
				if (nowP->sex==1)
				{	rn = drand48();
					if (rn <  fmodstart)	nowP->mod = 0;
					else			nowP->mod = 1;
				}
				break;
		}

	    /***Create the pointer links to previous and next individual***/
		nowP->prevP = prevP;
		nowP->nextP = NULL;
		if (i != param[sp].n)
		{	nowP->nextP = (repP) malloc(sizeof(REP));
			if (!nowP->nextP)
			{	printf("Unable to allocate memory to pointer: exiting from program\n");
				exit(1);
			}
		}

	    /***Update the pointers***/
		prevP = nowP;
		nowP  = nowP->nextP;
	}
	param[sp].lastP = prevP;

    /***Remove individual '0' so all individuals used are created by malloc***/
	nowP		 = param[sp].firstP;
	nowP		 = nowP->nextP;
	param[sp].firstP = nowP;
	nowP->prevP	 = NULL;

    /***Introduce a symbiont to hosts and modify host fitness accordingly (number of hosts with symbiont is Nsymstart)***/
    /***These hosts are at the top of the list***/
	if (Nsymstart > 0)
	{	i = 1;
		nowP = param[sp].firstP;
		while (i <= Nsymstart)
		{	nowP->sym  = 1;
			nowP->wsym = Wsym;
			switch(method_sym)
			{
				case 1:	/*symbiont makes intrinsic host birth rate B0*wsym (applies only to females)*/
					switch(nowP->sex)
					{	case 0: /*female*/
							nowP->b = param[sp].B0 * nowP->wsym;
							break;
						case 1: /*male*/
							nowP->b = 0.0;
							break;
					}
					break;

				case 2:	/*symbiont makes intrinsic host death rate D0*wsym (applies to females and males*/
					nowP->d = param[sp].D0 * nowP->wsym + param[sp].DD * param[sp].n;
					break;

				case 3:	/*symbiont makes intrinsic host death rate D0/wsym (applies to females and males*/
					nowP->d = param[sp].D0 / nowP->wsym + param[sp].DD * param[sp].n;
					break;

				case 4:	/*symbiont makes total host death rate D/wsym (applies to females and males*/
					nowP->d = (param[sp].D0 + param[sp].DD * param[sp].n) / nowP->wsym;
					break;
			}
			i = i + 1;
			nowP  = nowP->nextP;
		}
	}

    /***Output initial state of system***/
/*	if (flag_state == 1)	output_system_state();
*/
}


/****************************************************************************************************************
* Procedure for choosing probabilities of birth and death per unit time of each target host			*
****************************************************************************************************************/

void	individual_event_probabilities()
{
	int     sp, i;
	repP    targP;		/* pointer to target individual*/

    /***Loop for setting probabilties p.u.t. for individuals***/
	for (sp=1; sp<nspp; ++sp)
	{	targP = param[sp].firstP;
		while (targP != NULL)
		{	
		    /***Probability of birth p.u.t. fixed at param[sp].B0 * wsym;  targP->b is entered when targP is initialized***/

		    /***Probability of death p.u.t. depends on n and on whether host os sym+ ***/
			switch(method_sym)
			{
				case 1:	/*symbiont makes intrinsic host birth rate B0*wsym (applies only to females)*/
					targP->d = param[sp].D0 + param[sp].DD * param[sp].n;

				case 2:	/*symbiont makes intrinsic host death rate D0*wsym (applies to females and males)*/
					targP->d = param[sp].D0 * targP->wsym + param[sp].DD * param[sp].n;
					break;

				case 3:	/*symbiont makes intrinsic host death rate D0/wsym (applies to females and males)*/
					targP->d = param[sp].D0 / targP->wsym + param[sp].DD * param[sp].n;
					break;

				case 4:	/*symbiont makes total host death rate D/wsym (applies to females and males)*/
					targP->d = (param[sp].D0 + param[sp].DD * param[sp].n) / targP->wsym;
					break;
			}

		    /***End simulation if a negative birth or death rate has been encountered***/
                        if (targP->b < 0.0)
			{	printf("Negative birth rate encountered:  program halted\n");
				exit(1);
			}
                        if (targP->d < 0.0)
			{	printf("Negative death rate encountered:  program halted\n");
				exit(1);
			}

                        targP = targP->nextP;
		}
	}
}


/****************************************************************************************************************
* Procedure for choosing births and deaths at random.  Returns:							*
*   dt:    random time step from Gillespie algorithm								*
*   event: random choice of birth or death									*
*   pickP: pointer to random individual chosen for the event							*
****************************************************************************************************************/

void	choose_event(dt, event, pickP)
double	*dt;
char	*event;
repP	*pickP;
{
	FILE	*out;
	int	sp, i;
	repP	nowP;
	double	Bsum[nspp], Dsum[nspp];
	double	B, D, sum, rn, add;

    /***Probability per unit time of birth event***/
	for (sp=1; sp<nspp; ++sp)
	{	Bsum[sp] = 0.0;
		nowP  = param[sp].firstP;
		while (nowP != NULL)
		{	Bsum[sp] = Bsum[sp] + nowP->b;
			nowP  = nowP->nextP;
		}
	}

    /***Probability per unit time of death event***/
	for (sp=1; sp<nspp; ++sp)
	{	Dsum[sp] = 0.0;
		nowP  = param[sp].firstP;
		while (nowP != NULL)
		{	Dsum[sp] = Dsum[sp] + nowP->d;
			nowP  = nowP->nextP;
		}
	}

    /***Keep a record of the difference between total birth and death probability p.u.t. (just for interest)***/
	for (sp=1; sp<nspp; ++sp)
	param[sp].BD = Bsum[sp] - Dsum[sp];

    /***Sum of probabilities of all events***/
	B = 0.0;
	D = 0.0;
	for (sp=1; sp<nspp; ++sp)
	{	B = B + Bsum[sp];
		D = D + Dsum[sp];
	}
        sum = B + D;

    /***Random time to next event (exponential distribution of waiting times -- Gillespie)***/
	if (sum > error)	*dt = -log(1.0 - drand48())/sum;

    /***End simulation if sum of event probabilities = 0***/
        else
	{	printf("Sum of event probabilities = 0:	 program halted\n");
		exit(1);
	}

    /***Choose at random a birth or death event***/
	rn = drand48();
	if (rn <  B/sum)	*event = 'b';
	if (rn >= B/sum)	*event = 'd';

    /***Initialize variables for working out the individual to which the event happens***/
	rn   = drand48();
	add  = 0.0;

    /***Choose at random the individual to which event happens***/
	switch (*event)
	{
		case 'b': /*birth event; males have b=0 and cannot be chosen*/
			for (sp=1; sp<nspp; ++sp)
			{	nowP = param[sp].firstP;
				while (nowP != NULL)
				{	add = add + nowP->b/B;
					if (rn < add)
					{	add = add - nowP->b/B;
						goto label_b;
					}
					nowP  = nowP->nextP;
				}
			}
		label_b:break;

		case 'd': /*death event*/
			for (sp=1; sp<nspp; ++sp)
			{	nowP = param[sp].firstP;
				while (nowP != NULL)
				{	add = add + nowP->d/D;
					if (rn < add)
					{	add = add - nowP->d/D;
						goto label_d;
					}
                                        nowP  = nowP->nextP;
				}
			}
		label_d:break;
	}
        *pickP	 = nowP;

    /***Record of calculations***/
	if (flag_choose == 1)
	{	out = fopen("output.detail.dat", "a");
		fprintf(out, "\nTime:%9.6f  State at end of choose_event procedure\n", t+*dt);
		fprintf(out, "B:%6.3f  D:%6.3f  ", B, D);
		fprintf(out, "firstP:%d rn:%6.4f  add:%6.4f  event:%c  pickP:%d\n", param[nowP->sp].firstP, rn, add, *event, *pickP);
		fprintf(out, "  pickP:%d  sex:%d  sym:%d  mod:%d\n", *pickP, nowP->sex, nowP->sym, nowP->mod);
		fflush (out);
		fclose (out);
	}
}


/****************************************************************************************************************
* Choose random father for birth event (needed if transmission is biparental)					*
****************************************************************************************************************/

void	choose_dad(sp, dadP)
int	sp;		/* species of mother (same sp for father)	*/
repP	*dadP;		/* pointer to the individual chosen as father	*/
{
	FILE	*out;
	repP	nowP;
	int	ndad;
	double	dadsum;
	double	rn, add;

    /***Number of potential fathers***/
	dadsum = 0;
	nowP   = param[sp].firstP;
	while (nowP != NULL)
	{	if (nowP->sex == 1)
		{	ndad   += 1;	
			dadsum += 1.0;
		}
		nowP  = nowP->nextP;
	}
	if (ndad == 0)
	{	printf("no males in population -- program halted\n");
		exit(1);
	}

    /***Choose father at random***/
	rn   = drand48();
	add  = 0.0;
	nowP = param[sp].firstP;
	while (nowP != NULL)
	{	if (nowP->sex == 1)
		{	add = add + 1.0/dadsum;
			if (rn < add)
			{	add = add - 1.0/dadsum;
				goto label_dad;
			}
		}
		nowP  = nowP->nextP;
	}
	label_dad: *dadP = nowP;

    /***Record of calculations***/
	if (flag_choosedad == 1)
	{	out = fopen("output.detail.dat", "a");
		fprintf(out, "\nTime:%9.6f  State at end of choose_dad procedure\n", t);
		fprintf(out, "dadsum:%6.3f  ", dadsum);
		fprintf(out, "firstP:%d rn:%6.4f  add:%6.4f  dadP:%d\n", param[nowP->sp].firstP, rn, add, *dadP);
		fprintf(out, "  dadP:%d  sex:%d  sym:%d  mod:%d\n", *dadP, nowP->sex, nowP->sym, nowP->mod);
		fflush (out);
		fclose (out);
	}
}


/****************************************************************************************************************
* Construct symbiont status in newborn host individual								*
****************************************************************************************************************/
void	symbiont_transmission(pickP, dadP, birthP)
repP	pickP;		/* pointer to mother    */
repP	dadP;		/* pointer to father    */
repP	birthP;		/* pointer to offspring */
{
	FILE	*out;
	int	i;
	double	rn, store;

    /***Symbiont status of host offspring depends on mode of transmission***/
    /***For now just assume a single symbiont is transmitted, not a whole microbiome***/
	switch(method_transmit)
	{
		case 0: /*maternal: symbionts all from mother*/
			birthP->wsym = pickP->wsym;
			birthP->sym  = pickP->sym;
			break;

		case 1: /*paternal*/
			printf("Paternal transmission not implemented: exiting from program\n");
                        exit(1);
			break;

		case 2: /*biparental: host offspring gets the symbiont if one or both parents contain it*/
			if (pickP->sym==1 || dadP->sym==1)
			{	birthP->sym  = 1;
				birthP->wsym = Wsym;
			}
			else
			{	birthP->sym  = 0;
				birthP->wsym = 1.0;
			}
			break;

		case 3: /*evolving: host offspring gets the symbiont if mother has it, OR if father has */
			/*it and also has modfier gene mod:1; this allows transmission from father to evolve*/
			if (pickP->sym==1 || (dadP->sym==1 && dadP->mod==1))
                        {	birthP->sym  = 1;
				birthP->wsym = Wsym;
			}
                        else
                        {	birthP->sym  = 0;
				birthP->wsym = 1.0;
			}
                        break;
	}

    /***Details of calculation***/
	if (flag_transmit == 1)
	{	out = fopen("output.detail.dat", "a");
		fprintf(out, "\nTime:%9.6f  State at end of symbiont_transmission procedure\n", t);
		switch(method_transmit)
		{	case 0:	fprintf(out, "Maternal transmission: symbiont from mother\n");			break;
			case 1:	fprintf(out, "Paternal transmission: symbiont from father\n");			break;
			case 2:	fprintf(out, "Biparental transmission: symbiont from mother or father\n");	break;
		}
		fprintf(out, "  pickP: %d  sym:%d  mod:%2d\n",  pickP,  pickP->sym,  pickP->mod);
		fprintf(out, "  dadP:  %d  sym:%d  mod:%2d\n",   dadP,   dadP->sym,   dadP->mod);
		fprintf(out, "  birthP:%d  sym:%d         \n", birthP, birthP->sym             );
		fflush (out);
		fclose (out);
	}
}


/****************************************************************************************************************
* Birth of individual.												*
*   pickP:  pointer to mother											*
*   dadP:   pointer to father											*
*   birthP: pointer to newborn individual									*
*														*
* Calls: choose_dad(), needed for symbiotic status of newborn individual, if transmission biparental		*
* Calls: symbiont_transmission() to determine whether newborn individual is sym+ or sym-			*
*														*
* The new individual is placed immediately after the mother in the list;  its properties are defined in this	*
* procedure. This includes whether a newborn male will be able to transmit a symbiont to its offspring.		*
****************************************************************************************************************/

void	birth(pickP)
repP	pickP;		/* pointer to mother, already chosen to give birth */
{
	FILE	*out;
	int	sp;
	double	rn;
	repP	birthP; /* pointer to newborn individual		*/
	repP	dadP;	/* pointer to individual chosen as father	*/


    /***Allocate memory for newborn individual***/
	birthP = NULL;
	birthP = (repP) malloc(sizeof(REP));
	if (!birthP)
	{	printf("Unable to allocate memory to pointer: exiting from program\n");
		exit(1);
	}

    /***WARNING!  It is essential to alter the pointers in the correct sequence below***/

    /***If the new individual is at the end of the list, the pointer to the next individual is set to NULL***/
    /*** and a record of the tail pointer is made***/
	if (pickP->nextP == NULL)
	{	birthP->nextP = NULL;
			param[pickP->sp].lastP = birthP;
	}

    /***If the new individual is not at the end of the list, it is linked to the next individual ***/
	else
	{	birthP->nextP = pickP->nextP;
			pickP->nextP->prevP = birthP;
	}

    /***Lastly the links are made between the parent and the new individual***/
	birthP->prevP = pickP;
	pickP->nextP  = birthP;


/***PROPERTIES OF NEWBORN INDIVIDUAL***/

    /***Neonate: same species as mother***/
	birthP->sp = pickP->sp;

    /***Neonate: father chosen at random in choose_dad(): needed for symbiont transmission and for evolution of modifier gene***/
	choose_dad(pickP->sp, &dadP);

    /***Neonate: symbiont status computed in symbiont_transmission() (depends on mode of transmission)***/
	symbiont_transmission(pickP, dadP, birthP);

    /***Neonate: sex chosen at random with probability 0.5***/
	rn = drand48();
	if (rn <  0.5)	birthP->sex = 0;	/* female */
	if (rn >= 0.5)	birthP->sex = 1;	/* male   */

    /***Neonate: modifier gene for male transmission (only used in case 3; code should not be needed for case 0 and 2)***/
	switch (method_transmit)
	{	case 0: /*maternal: male transmission always switched off*/
			if (birthP->sex==0)	birthP->mod = 99;
			if (birthP->sex==1)	birthP->mod = 0;
			break;
		case 1: /*paternal (not implemented)*/
			break;
		case 2: /*biparental: male transms=ission always switched on*/
			if (birthP->sex==0)	birthP->mod = 99;
			if (birthP->sex==1)	birthP->mod = 1;
			break;
		case 3: /*evolution: modifier gene inherited from father (on Y chromosome)*/
			/*daughters do not carry the modifier gene; dummy value inserted*/
			if (birthP->sex==0)	birthP->mod = 99;
			if (birthP->sex==1)	birthP->mod = dadP->mod;
			break;
	}

    /***Neonate: birth and death rates depend on symbiont status (only females give birth)***/
	sp = birthP->sp;
	switch(method_sym)
	{
		case 1:	/*symbiont makes intrinsic host birth rate B0*wsym (applies only to females)*/
			switch(birthP->sex)
			{	case 0: /*female*/
					birthP->b = param[sp].B0 * birthP->wsym;
					birthP->d = param[sp].D0 + param[sp].DD * param[sp].n;
					break;
				case 1: /*male*/
					birthP->b = 0.0;
					birthP->d = param[sp].D0 + param[sp].DD * param[sp].n;
					break;
			}
			break;

		case 2:	/*symbiont makes intrinsic host death rate D0*wsym (applies to females and males)*/
			switch(birthP->sex)
			{	case 0: /*female*/
					birthP->b = param[sp].B0;
					birthP->d = param[sp].D0 * birthP->wsym + param[sp].DD * param[sp].n;
					break;
				case 1: /*male*/
					birthP->b = 0.0;
					birthP->d = param[sp].D0 * birthP->wsym + param[sp].DD * param[sp].n;
					break;
			}
			break;

		case 3:	/*symbiont makes intrinsic host death rate D0/wsym (applies to females and males)*/
			switch(birthP->sex)
			{	case 0: /*female*/
					birthP->b = param[sp].B0;
					birthP->d = param[sp].D0 / birthP->wsym + param[sp].DD * param[sp].n;
					break;
				case 1: /*male*/
					birthP->b = 0.0;
					birthP->d = param[sp].D0 / birthP->wsym + param[sp].DD * param[sp].n;
					break;
			}
			break;

		case 4:	/*symbiont makes total host death rate D/wsym (applies to females and males)*/
			switch(birthP->sex)
			{	case 0: /*female*/
					birthP->b = param[sp].B0;
					birthP->d = (param[sp].D0 + param[sp].DD * param[sp].n) / birthP->wsym;
					break;
				case 1: /*male*/
					birthP->b = 0.0;
					birthP->d = (param[sp].D0 + param[sp].DD * param[sp].n) / birthP->wsym;
					break;
			}
			break;
	}

    /***The number of individuals of this species is increased by 1***/
	param[pickP->sp].n = param[pickP->sp].n + 1;

    /***Details of calculation***/
	if (flag_birth == 1)
	{	out = fopen("output.detail.dat", "a");
		fprintf(out, "\nTime:%9.6f  State at end of birth procedure\n", t);
		fprintf(out, "  pickP: %d   pickP->sex:%d   pickP->sym:%d   pickP->mod:%2d   pickP->wsym:%6.4f\n",  pickP,  pickP->sex,  pickP->sym,  pickP->mod,  pickP->wsym);
		fprintf(out, "  dadP:  %d    dadP->sex:%d    dadP->sym:%d    dadP->mod:%2d    dadP->wsym:%6.4f\n",   dadP,   dadP->sex,   dadP->sym,   dadP->mod,   dadP->wsym);
		fprintf(out, "  birthP:%d  birthP->sex:%d  birthP->sym:%d  birthP->mod:%2d  birthP->wsym:%6.4f\n", birthP, birthP->sex, birthP->sym, birthP->mod, birthP->wsym);
		fflush (out);
		fclose (out);
	}
}


/****************************************************************************************************************
* Death of individual                                                                                           *
****************************************************************************************************************/

void	death(pickP)
repP	pickP;		/* pointer to the individual chosen for event */
{
	FILE	*out;

    /***If the dead individual is in the interior of the list, the the two adjoining individauls are linked***/
	if (pickP->nextP !=NULL && pickP->prevP !=NULL)
	{	pickP->prevP->nextP = pickP->nextP;
		pickP->nextP->prevP = pickP->prevP;
	}

    /***If the dead individual is at the head of the list, the next individual becomes the head***/
	if (pickP->prevP == NULL && pickP->nextP != NULL)
	{	pickP->nextP->prevP = NULL;
		param[pickP->sp].firstP = pickP->nextP;
	}

    /***If the dead individual is at the tail of the list, the previous individual becomes the tail***/
	if (pickP->nextP == NULL && pickP->prevP != NULL)
	{	pickP->prevP->nextP = NULL;
		param[pickP->sp].lastP = pickP->prevP;
	}

    /***If the dead individual is at the head AND the tail of the list, the species is extinct***/
	if (pickP->nextP == NULL && pickP->prevP == NULL)
	{	param[pickP->sp].firstP = NULL;
		param[pickP->sp].lastP	= NULL;
	}

    /***The number of individuals of this species is decreased by 1***/
	param[pickP->sp].n = param[pickP->sp].n - 1;

    /***Record of calculations***/
	if (flag_death == 1)
	{	out = fopen("output.detail.dat", "a");
		fprintf(out, "\nTime:%9.6f  State at end of death procedure\n", t);
		fprintf(out, "  pickP: %d  pickP->sex:%d  pickP->sym:%d  pickP->mod:%d  pickP->wsym:%6.4f\n",  pickP,  pickP->sex,  pickP->sym,  pickP->mod,  pickP->wsym);
		fflush (out);
		fclose (out);
	}
}


/****************************************************************************************************************
* Calcualtes various measures at each time, and puts them into output timeseries file				*
****************************************************************************************************************/

void	output_timeseries(dt, event, pickP)
double	dt;		/* time increment to the event			*/
char	event;		/* event chosen in each iteration: b, d		*/
repP	pickP;		/* pointer to the individual chosen for event	*/
{
	FILE	*out;
	int	sp;
	repP	nowP;		/* pointer to individuals in list		*/

	sp = 1;

    /***Breakdown to numbers by sex and symbiont status (just for checking single timeseries)***/
	param[sp].nfsym1 = 0;
	param[sp].nmsym1 = 0;
	param[sp].nfsym0 = 0;
	param[sp].nmsym0 = 0;
	nowP = param[sp].firstP;
	while (nowP != NULL)
	{	if (nowP->sex==0 && nowP->sym==1)	param[sp].nfsym1 += 1;
		if (nowP->sex==1 && nowP->sym==1)	param[sp].nmsym1 += 1;
		if (nowP->sex==0 && nowP->sym==0)	param[sp].nfsym0 += 1;
		if (nowP->sex==1 && nowP->sym==0)	param[sp].nmsym0 += 1;
		nowP = nowP->nextP;
	}
	param[sp].fsym = (double)param[sp].nsym / (double)param[sp].n;

    /***Check sex ratio (just for output)***/
	param[sp].nf = 0;
	param[sp].nm = 0;
	nowP = param[sp].firstP;
	while (nowP != NULL)
	{	switch(nowP->sex)
		{	case 0: /*female*/
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
	{	if (nowP->sex==1)
		{	if(nowP->mod==0)	param[sp].nmat += 1;
			else			param[sp].nbip += 1;
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
	{	if (nowP->sex==1)
		{	if(nowP->sym==0 && nowP->mod==0)	param[sp].nsym0mat += 1;
			if(nowP->sym==1 && nowP->mod==0)	param[sp].nsym1mat += 1;
			if(nowP->sym==0 && nowP->mod==1)	param[sp].nsym0bip += 1;
			if(nowP->sym==1 && nowP->mod==1)	param[sp].nsym1bip += 1;
		}
		nowP = nowP->nextP;
	}
	param[sp].dis = (double)(param[sp].nsym0mat * param[sp].nsym1bip - param[sp].nsym1mat *param[sp].nsym0bip)
		      / (double)(param[sp].nm * param[sp].nm);

    /***Record results on screen***/
	printf("t:%9.6f	sp:%d	%6d  %c  ", t, pickP->sp, nevent, event);
	printf("B-D:%7.3f  ", param[sp].BD);
	printf("%4d  %4d:%4d  %6.4f    %6.4f  ", param[sp].n, param[sp].nf, param[sp].nm, param[sp].fsym, param[sp].fmat);
/*	printf("  address:%d  ", pickP);
*/	printf("\n");

    /***Record results to file***/
	out = fopen("output.timeseries.dat", "a");
	if (t-dt < error)
	fprintf(out, "#time     sp  nevent event   n    nf    nm  fsymbiont  nmat  nbip fmaternal   mat- mat+ bip- bip+      diseq\n");
	fprintf(out, "%9.6f  %d  %6d   %c  ", t, pickP->sp, nevent, event);
	fprintf(out, "%4d  %4d  %4d  %6.4f    ", param[sp].n,    param[sp].nf,   param[sp].nm,   param[sp].fsym);
	fprintf(out, " %4d  %4d  %6.4f    ",     param[sp].nmat, param[sp].nbip, param[sp].fmat);
	fprintf(out, " %4d %4d %4d %4d    ",     param[sp].nsym0mat, param[sp].nsym1mat, param[sp].nsym0bip, param[sp].nsym1bip);
	fprintf(out, " %7.4f  ",                 param[sp].dis);
	fprintf(out, "\n");
	fflush (out);
	fclose (out);
}


/****************************************************************************************************************
* Runs a realization of the stochastic process. At each time step:						*
* Assigns to each individual a probability p.u.t. of a birth and death in: individual_event_probabilities()	*
* Chooses a random time, event and individual in: choose_event()						*
* Updates the state of the system, calling either:								*
*     birth()													*
*     death()													*
****************************************************************************************************************/

void	simulate()
{
	FILE	*out;
	int	i, j, sp;

	double	dt;		/* time increment to the event				*/
	char	event;		/* event chosen in each iteration: b, d			*/
	repP	pickP;		/* pointer to the individual chosen for event		*/
	repP	nowP;

	sp = 1;

    /***Initialisation***/
	initialize_hosts();

    /***Loop through births and deaths until sym+ hosts reach a frequency ~0.0 or ~1.0 (tmax is a backstop)***/
	t = 0.0;
	nevent = 0;
	param[sp].fsym = 0.1;
/*	while (t<=tmax && param[sp].fsym>error && param[sp].fsym<1.0-error && param[sp].n>1)
*/	while (t<=tmax)
	{
            /***Output state of system***/
		if (flag_state == 1)	output_system_state();

            /***Calculate the birth and death probabilities per unit time of each individual***/
		individual_event_probabilities();

            /***Choose a birth or death event at random***/
		choose_event(&dt, &event, &pickP);

            /***Update the time and event number***/
		t = t + dt;
		nevent += 1;

            /***Update the state of the system***/
		switch	(event)
		{	case 'b':	birth(pickP);		break;
			case 'd':	death(pickP);		break;
		}

	    /***Number and frequency of hosts with symbiont (needed in main() below)***/
		param[sp].nsym = 0.0;
		nowP = param[sp].firstP;
		while (nowP != NULL)
		{	if (nowP->sym==1)	param[sp].nsym += 1;
			nowP = nowP->nextP;
		}
		param[sp].fsym = (double)param[sp].nsym / (double)param[sp].n;

	    /***Output to file output.timeseries.dat if required (calculates various extra measures along the way)***/
		if (flag_timeseries == 1)
		output_timeseries(dt, event, pickP);

	    /***If event was death, free the memory previously given to the individual***/
		if (event == 'd')       free(pickP);
	}
}


/****************************************************************************************************************
* Procedure to accumulate frequency distribution of Wsym's that allowing invasion of symbiont			*
****************************************************************************************************************/

void	fitness_freq(iter)
int	iter;
{
	int	sp, bin;
	double	binmin, binmax;

    /***At start set the bins***/
	if (iter == 0)
	{	for (bin=0; bin<nbin; ++bin)
		{	freq[0][bin] = bin * binwidth;
			freq[1][bin] = 0.0;
		}
	}

    /***Thereafter, accumulate frequencies of Wsym***/
	else
	{	sp  = 1;
		bin = (int) (Wsym / binwidth);
		freq[1][bin] += 1.0;
		printf("Wsym:%6.4f  D0/Wsym:%6.4f  ", Wsym, param[sp].D0/Wsym);
		printf("freq[%2d]:%4.0f \n", bin, freq[1][bin]);
	}
}


/****************************************************************************************************************
* Main body of program												*
*   Loops through a specified number of realizations (niterate) of the stochastic process.			*
*   It calls procedure simulate(), which runs the realization.							*
*   Checks whether the realization meets a condition for successful invasion by the symbiont.			*
*   If successful, records in output.invasions.dat, and adds w into a frequency distribution in fitness_freq()	*
****************************************************************************************************************/

int	main()
{
	FILE	*out;
	int	iter, bin, sp, count;

    /***Initialize output files (ensures they are empty at the start)***/
	fclose(fopen("output.detail.dat",       "w"));
	fclose(fopen("output.timeseries.dat",   "w"));
	fclose(fopen("output.invasions.dat",    "w"));
	fclose(fopen("output.freq.distrib.dat", "w"));

    /***Seed the random number generator***/
	printf("\nSeed for random number generator (6 digits):	  ");
	scanf("%6d", &seed);
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

    /***Set up output file for bin_invasion()***/
	iter = 0;
	fitness_freq(iter);

    /***Loop for iterating stochastic process with random Wsym***/
	sp = 1;
	count = 0;
	for (iter=1; iter<=niterate; ++iter)
	{
	    /***Random number from normal distribution with constraints 0: <= Wsym <= 4***/
/*		Wsym = -1.0;
		while (Wsym < error || Wsym > 2.5)
		Wsym = random_normal()*0.30 + 1.0;
		Wsym = 0.2;
		Wsym = 0.6;
		Wsym = 1.0;
*/		Wsym = 1.5;

	    /***Realisation of stochastic birth-death process***/
		simulate();

	    /***Record of Wsym's that allow symbiont invasion***/
		if (param[sp].fsym >= fsuccess)
		{
			count += 1;
			printf("iteration:%3d   count:%3d   ", iter, count);

		    /***Add to list of symbiont invasions***/
			out = fopen("output.invasions.dat", "a");
			fprintf(out, "%6d   %8.6f   %4d   %4d   %6.4f   %6.2f  %6d", iter, Wsym, param[sp].n, param[sp].nsym, param[sp].fsym, t, nevent);
			fprintf(out, "\n");
			fflush (out);
			fclose (out);

		    /***Accumulate frequency distribution of Wsym values allowing symbiont invasion***/
			fitness_freq(iter);
			out = fopen("output.freq.distrib.dat", "w");
			for (bin=0; bin<nbin; ++bin)
			{	fprintf(out, "%6.4f   %4.0f", freq[0][bin], freq[1][bin]);
				fprintf(out, "\n");
			}
			fflush (out);
			fclose (out);
		}
	}

/*	out = fopen("output.detail.dat", "a");
	fprintf(out, "\nNumber of successful invasions: %6d\n", count);
	fflush (out);
	fclose (out);
*/
}
