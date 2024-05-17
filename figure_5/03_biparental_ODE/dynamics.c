/****************************************************************************************
* This program integrates ODEs for host-symbiont systems with maternal and biparental	*
* transmission of the symbiont								*
*                                                                                       *
* Males have two genes at a locus on the Y chromosome giving different transmission	*
* of symbiont from father to infant:							*
*   gene M0: beta_M0 = 0.0 for maternal  transmission					*
*   gene M1: beta_M1 = 1.0 for biparental transmission					*
*                                                                                       *
* 30.12.21.										*
*                                                                                       *
* To match the stochastic realisation, this run uses w=0.37, E0=0.1                     *
* Initial values are taken from the stochastic realisation after the first event        *
* and are given at the start of the procedure integrate().                              *
*                                                                                       *
* 10.03.24                                                                              *
****************************************************************************************/

#include	<stdio.h>
#include	<stdlib.h>
#include	<math.h>

#define         error   1e-9            /* rounding error for type double	*/

#define		B0	4.0		/* birth rate per host individual	*/
#define		D0	1.0		/* death rate intrinsic			*/
#define		DD	0.001		/* death rate density dependent		*/
#define		E0	0.1/*0.50*/	/* encounter rate (horizontal)		*/
#define		V	656.0/*1000.0*/	/* system size				*/

#define		w	0.37/*0.4 0.6*/	/* host sym+ fitness			*/
#define		alpha	1.0		/* maternal transmission proportion	*/
#define		beta_M0	0.0		/* paternal transmission with gene M0	*/
#define		beta_M1	1.0		/* paternal transmission with gene M1	*/

#define		tmax	100.0		/* time period for integration		*/
#define		dt	0.01		/* time step   for integration		*/

struct	density
	{	double	fsym0;		/* host density: female sym-		*/
		double	fsym1;		/*		 female sym+		*/
		double	msym0_M0;	/*		 male   sym-	gene M0	*/
		double	msym1_M0;	/*		 male   sym+	gene M0	*/
		double	msym0_M1;	/*		 male   sym-	gene M1	*/
		double	msym1_M1;	/*		 male   sym+	gene M1	*/
	};

struct	birth
	{	double	fsym0;		/* host birth rate: female sym-		*/
		double	fsym1;		/*		 female sym+		*/
		double	msym0_M0;	/*		 male   sym-	gene M0	*/
		double	msym1_M0;	/*		 male   sym+	gene M0	*/
		double	msym0_M1;	/*		 male   sym-	gene M1	*/
		double	msym1_M1;	/*		 male   sym+	gene M1	*/
	};

struct	death
	{	double	fsym0;		/* host death rate: female sym-		*/
		double	fsym1;		/*		 female sym+		*/
		double	msym0_M0;	/*		 male   sym-	gene M0	*/
		double	msym1_M0;	/*		 male   sym+	gene M0	*/
		double	msym0_M1;	/*		 male   sym-	gene M1	*/
		double	msym1_M1;	/*		 male   sym+	gene M1	*/
	};

struct	encounter			/* host encounter rate (horiz trans)	*/
	{	double	fsym0;		/* host death rate: female sym-		*/
		double	fsym1;		/*		 female sym+		*/
		double	msym0_M0;	/*		 male   sym-	gene M0	*/
		double	msym1_M0;	/*		 male   sym+	gene M0	*/
		double	msym0_M1;	/*		 male   sym-	gene M1	*/
		double	msym1_M1;	/*		 male   sym+	gene M1	*/
	};

/***Host population densities***/
struct	density		x;

/***Host birth rates***/
struct	birth		b;

/***Host death rates***/
struct	death		d;

/***Host encounter rates***/
struct	encounter	e;


/****************************************************************************************
* Output										*
****************************************************************************************/
void	output(t)
double	t; 
{
	FILE    *out;
	double	X, sym1, M0, D, xm;

	out = fopen("output.timeseries.dat", "a");

    /***Host population size and frequency of carriers of symbiont***/
	X    = x.fsym0 + x.fsym1 + x.msym0_M0 + x.msym1_M0 + x.msym0_M1 + x.msym1_M1;
	sym1 = (x.fsym1 + x.msym1_M0 + x.msym1_M1) / X;

    /***Frequency of M0 gene***/
	M0   = (x.msym0_M0+ x.msym1_M0)/ (x.msym0_M0+ x.msym1_M0 + x.msym0_M1+ x.msym1_M1);

    /***Coefficient of disequilibrium***/
	xm = x.msym0_M0 + x.msym1_M0 + x.msym0_M1 + x.msym1_M1;
        D  = (x.msym0_M0 * x.msym1_M1 - x.msym1_M0 * x.msym0_M1) / (xm * xm);

    /***Header for output file***/
	if (t < error)
	fprintf(out, "t         VX       fr_sym1   fr_M0      D          x.fsym0   x.fsym1  x.msym0M0  x.msym1M0  x.msym0M1  x.msym1M1   \n");

    /***Output at each time step***/
	fprintf(out, "%6.3f    %6.1f   %6.4f   %6.4f   %6.4f       ",       t, V*X, sym1, M0, D);  
	fprintf(out, "%6.1f   %6.1f   %6.1f     %6.1f     %6.1f     %6.1f", V*x.fsym0, V*x.fsym1, V*x.msym0_M0, V*x.msym1_M0, V*x.msym0_M1, V*x.msym1_M1);  
	fprintf(out, "\n");
        fflush(out);
	fclose(out);

}


/****************************************************************************************
* Birth rates (total -- not per capita)							*
* Uses the mating table (s x s, etc)							*
****************************************************************************************/
void	birth_rates()
{
											/* f x m */

	b.fsym0 = B0/2.0 * 1.0/(x.msym0_M0 + x.msym0_M1 + x.msym1_M0 + x.msym1_M1) 
		* (  x.fsym0 * x.msym0_M0						/* s x s */
		   + x.fsym0 * x.msym0_M1						/* s x s */

		   + x.fsym1 * x.msym0_M0 * (1.0-alpha  )				/* S x s */
		   + x.fsym1 * x.msym0_M1 * (1.0-alpha  )				/* S x s */

		   + x.fsym0 * x.msym1_M0 * (1.0-beta_M0)				/* s x S */
		   + x.fsym0 * x.msym1_M1 * (1.0-beta_M1)				/* s x S */

		   + x.fsym1 * x.msym1_M0 * (1.0-alpha  ) * (1.0-beta_M0)		/* S x S */
		   + x.fsym1 * x.msym1_M1 * (1.0-alpha  ) * (1.0-beta_M1)		/* S x S */
		  );

	b.fsym1 = B0/2.0 * 1.0/(x.msym0_M0 + x.msym0_M1 + x.msym1_M0 + x.msym1_M1)
		* (  0.0								/* s x s */

		   + x.fsym1 * x.msym0_M0 * alpha					/* S x s */
		   + x.fsym1 * x.msym0_M1 * alpha					/* S x s */

		   + x.fsym0 * x.msym1_M0 * beta_M0					/* s x S */
		   + x.fsym0 * x.msym1_M1 * beta_M1					/* s x S */

		   + x.fsym1 * x.msym1_M0 * (alpha + beta_M0 - alpha*beta_M0)		/* S x S */
		   + x.fsym1 * x.msym1_M1 * (alpha + beta_M1 - alpha*beta_M1)		/* S x S */
		  );

	b.msym0_M0 = B0/2.0 * 1.0/(x.fsym0 + x.fsym1)
		   * (  x.fsym0 * x.msym0_M0						/* s x s */
		      + x.fsym1 * x.msym0_M0 * (1.0-alpha  )				/* S x s */
		      + x.fsym0 * x.msym1_M0 * (1.0-beta_M0)				/* s x S */
		      + x.fsym1 * x.msym1_M0 * (1.0-alpha  ) * (1.0-beta_M0)		/* S x S */
		     );

	b.msym1_M0 = B0/2.0 * 1.0/(x.fsym0 + x.fsym1) 
		   * (  0.0								/* s x s */
		      + x.fsym1 * x.msym0_M0 * alpha					/* S x s */
		      + x.fsym0 * x.msym1_M0 * beta_M0					/* s x S */
		      + x.fsym1 * x.msym1_M0 * (alpha + beta_M0 - alpha*beta_M0)	/* S x S */
		     );

	b.msym0_M1 = B0/2.0 * 1.0/(x.fsym0 + x.fsym1)
		   * (  x.fsym0 * x.msym0_M1						/* s x s */
		      + x.fsym1 * x.msym0_M1 * (1.0-alpha  )				/* S x s */
		      + x.fsym0 * x.msym1_M1 * (1.0-beta_M1)				/* s x S */
		      + x.fsym1 * x.msym1_M1 * (1.0-alpha  ) * (1.0-beta_M1)		/* S x S */
		     );

	b.msym1_M1 = B0/2.0 * 1.0/(x.fsym0 + x.fsym1) 
		   * (  0.0								/* s x s */
		      + x.fsym1 * x.msym0_M1 * alpha					/* S x s */
		      + x.fsym0 * x.msym1_M1 * beta_M1					/* s x S */
		      + x.fsym1 * x.msym1_M1 * (alpha + beta_M1 - alpha*beta_M1)	/* S x S */
		     );
}


/****************************************************************************************
* Death rates (total -- not per capita)							*
****************************************************************************************/
void	death_rates()
{
	double	X;

    /***Total populaton density***/
	X = x.fsym0 + x.fsym1 + x.msym0_M0 + x.msym1_M0 + x.msym0_M1 + x.msym1_M1;

	d.fsym0    = 1.0   * x.fsym0    * (D0 + DD * V * X);
	d.fsym1    = 1.0/w * x.fsym1    * (D0 + DD * V * X);
	d.msym0_M0 = 1.0   * x.msym0_M0 * (D0 + DD * V * X);
	d.msym1_M0 = 1.0/w * x.msym1_M0 * (D0 + DD * V * X);
	d.msym0_M1 = 1.0   * x.msym0_M1 * (D0 + DD * V * X);
	d.msym1_M1 = 1.0/w * x.msym1_M1 * (D0 + DD * V * X);
}


/****************************************************************************************
* Encounter rate with symbiont in the environment -- horizontal (total--not per capita) *
****************************************************************************************/
void	encounter_rates()
{
	e.fsym0    = E0 * x.fsym0;
 	e.msym0_M0 = E0 * x.msym0_M0;
	e.msym0_M1 = E0 * x.msym0_M1;
}


/****************************************************************************************
* Integration loop									*
****************************************************************************************/
void	integrate()
{
	double		t;	/* Current time		*/
	struct	density	xnew;

    /***Initial densties***/
/*	x.fsym0    = 0.45;
	x.fsym1    = 0.05;
	x.msym0_M0 = 0.09;
	x.msym1_M0 = 0.01;
	x.msym0_M1 = 0.36;
	x.msym1_M1 = 0.04;
*/
    /***Initial frequencies:   sym1: 0.01   M0: 0.4/2 (1/2 because it's only in males)***/
/*	x.fsym0    = 0.495;
	x.fsym1    = 0.005;
	x.msym0_M0 = 0.396;		
	x.msym1_M0 = 0.004;
	x.msym0_M1 = 0.099;
	x.msym1_M1 = 0.001;
*/
    /***Initial frequencies:   sym1: 0.01   M0: 0.25/2 (1/2 because it's only in males)***/
/*	x.fsym0    = 0.495;
	x.fsym1    = 0.005;
	x.msym0_M0 = 0.2475;		
	x.msym1_M0 = 0.0025;
	x.msym0_M1 = 0.2475;
	x.msym1_M1 = 0.0025;
*/
    /***Initial frequencies:      sym1: ~0.01  ***/
/*	x.fsym0    = 0.495;
	x.fsym1    = 0.005;
	x.msym0_M0 = 0.37125;
	x.msym1_M0 = 0.00375;
	x.msym0_M1 = 0.12375;		
	x.msym1_M1 = 0.00125;
*/
    /***Initial frequencies:   these values taken from stochastic output after first event ***/
     /* n:        655                                                        */
     /* nfemales: 332                    nfemales(sym-): 292  / 655 = 0.4458 */
     /*                                  nfemales(sym+): 40   / 655 = 0.0611 */
     /* nmales:   323 nmales(sym-): 281  nmales(sym-)M1: 275  / 655 = 0.4198 */
     /*                                  nmales(sym-)M0: 6    / 655 = 0.0092 */
     /*               nmales(sym+): 42   nmales(sym+)M1: 42   / 655 = 0.0641 */
     /*                                  nmales(sym+)M0: 0    / 655 = 0      */
	x.fsym0    = 0.4458;
	x.fsym1    = 0.0611;
	x.msym0_M1 = 0.4198;
	x.msym0_M0 = 0.0092;		
	x.msym1_M1 = 0.0641;
	x.msym1_M0 = 0.0000;

    /***Set time at start***/
        t       = 0.0;
/*	toutput = 0;
*/
    /***Iteration for evolution -- stops at tmax***/
        while (t<tmax)
        {
	    /***Output***/
		output(t);

	    /***Compute birth, death and encounter rates***/
		birth_rates();
		death_rates();
		encounter_rates();

	    /***Integration step***/
		xnew.fsym0    = x.fsym0    + dt * (b.fsym0    - d.fsym0    - e.fsym0); 
		xnew.fsym1    = x.fsym1    + dt * (b.fsym1    - d.fsym1    + e.fsym0); 
		xnew.msym0_M0 = x.msym0_M0 + dt * (b.msym0_M0 - d.msym0_M0 - e.msym0_M0); 
		xnew.msym1_M0 = x.msym1_M0 + dt * (b.msym1_M0 - d.msym1_M0 + e.msym0_M0); 
		xnew.msym0_M1 = x.msym0_M1 + dt * (b.msym0_M1 - d.msym0_M1 - e.msym0_M1); 
		xnew.msym1_M1 = x.msym1_M1 + dt * (b.msym1_M1 - d.msym1_M1 + e.msym0_M1); 

	    /***Update densities***/
		x.fsym0    = xnew.fsym0;
		x.fsym1    = xnew.fsym1;
		x.msym0_M0 = xnew.msym0_M0;
		x.msym1_M0 = xnew.msym1_M0;
		x.msym0_M1 = xnew.msym0_M1;
		x.msym1_M1 = xnew.msym1_M1;

	    /***Update time***/
		t = t + dt;
/*		toutput = toutput + dt;
*/
	}

}


/****************************************************************************************
* Main body of program									*
****************************************************************************************/
int	main()
{

    /***Initialize output file (ensures it is empty at the start)***/
	fclose(fopen("output.timeseries.dat", "w"));

    /***Integration***/
	integrate();
}

