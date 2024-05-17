/****************************************************************************************
* This program integrates ODEs for host-symbiont systems with maternal and biparental	*
* transmission of the symbiont								*
*                                                                                       *
* 30.12.21.										*
****************************************************************************************/

#include	<stdio.h>
#include	<stdlib.h>
#include	<math.h>

#define	mu	1.0
#define	sigma	0.3
#define	dx	0.05
#define	pi	3.1415927


/****************************************************************************************
* Main body of program									*
****************************************************************************************/
int	main()
{
	FILE	*out;
	double	x, y, f, sum;

    /***Initialize output file (ensures it is empty at the start)***/
        out = fopen("output.gaussian.dat", "w");

	x = 0.0;
	sum = 0.0;
	while (x < 2.2)
	{
		f = (x - mu) / sigma;
		y = exp(-0.5*f*f) / (sigma*sqrt(2*pi)) * 5000.0 * dx;

		sum = sum + y;

		fprintf(out, "%6.2f   %10.6f   %10.5f \n", x, y, sum);
		x = x + dx;
	}

	fflush(out);
	fclose(out);
}

