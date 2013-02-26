#include <king/rand.h>

/*
float ran0(long *idum)
{
	long k;
	float ans;
	*idum ^= MASK;						// XORing with MASK allows use of zero and other
	k=(*idum)/IQ;						// simple bit patterns for idum.
	*idum=IA*(*idum-k*IQ)-IR*k;				// Compute idum=(IA*idum) % IM without over-
	if (*idum < 0)
		*idum += IM;					// ﬂows by Schrage’s method.
	ans=AM*(*idum);						// Convert idum to a ﬂoating result.
	*idum ^= MASK;						// Unmask before return.
	return ans;
}

float ran1(long *idum)
{
	int j;
	long k;
	static long iy=0;
	static long iv[NTAB];
	float temp;
	if (*idum <= 0 || !iy)
	{							// Initialize.
		if (-(*idum) < 1)
			*idum=1;				// Be sure to prevent idum = 0.
		else
			*idum = -(*idum);
		for (j=NTAB+7;j>=0;j--)
		{						// Load the shuffle table (after 8 warm-ups).
			k=(*idum)/IQ;
			*idum=IA*(*idum-k*IQ)-IR*k;
			if (*idum < 0)
				*idum += IM;
			if (j < NTAB)
				iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ;						// Start here when not initializing.
	*idum=IA*(*idum-k*IQ)-IR*k;				// Compute idum=(IA*idum) % IM without over-
	if (*idum < 0)
		*idum += IM;					// ﬂows by Schrage’s method.
	j=iy/NDIV;						// Will be in the range 0..NTAB-1.
	iy=iv[j];						// Output previously stored value and reﬁll the
	iv[j] = *idum;						// shuffle table.
	if ((temp=AM*iy) > RNMX)
		return RNMX;					// Because users don’t expect endpoint values.
	else
		return temp;
}
*/

float ran2(long *idum)
{
	int         j;
	long        k;
	static long idum2 = 123456789;
	static long iy    = 0;
	static long iv[NTAB];
	float       temp;

	if( *idum <= 0 )
	{							// Initialize.
		if( -(*idum) < 1 )
			*idum = 1;				// Be sure to prevent idum = 0.
		else
			*idum = -(*idum);

		idum2  = (*idum);

		for(j = NTAB + 7; j >= 0; j--)
		{						// Load the shuﬄe table (after 8 warm-ups).
			k     = (*idum) / IQ1;
			*idum = IA1 * (*idum - k * IQ1) - k * IR1;
			if( *idum < 0 )
				*idum += IM1;
			if( j < NTAB )
				iv[j]  = *idum;
		}
		iy     = iv[0];
	}
	k     = (*idum) / IQ1;					// Start here when not initializing.
	*idum = IA1 * (*idum - k * IQ1) - k * IR1;		// Compute idum=(IA1*idum) % IM1 without
	if(*idum < 0)
		*idum += IM1;					// overﬂows by Schrage’s method.
	k     = idum2 / IQ2;
	idum2 = IA2 * (idum2 - k * IQ2) - k * IR2;		// Compute idum2=(IA2*idum) % IM2 likewise.
	if(idum2 < 0)
		idum2 += IM2;
	j     = iy / NDIV;					// Will be in the range 0..NTAB-1.
	iy    = iv[j] - idum2;					// Here idum is shuﬄed, idum and idum2 are
	iv[j] = *idum;						// combined to generate output.
	if(iy < 1)
		iy    += IMM1;
	if( (temp = AM * iy) > RNMX)
		return RNMX;					// Because users don’t expect endpoint values.
	else
		return temp;
}

float ran3(long *idum)
{
	static int  inext,
		    inextp;
	static long ma[56];					// The value 56 (range ma[1..55]) is special and
	static int  iff = 0; 					// should not be modiﬁed; see Knuth.
	long        mj,
	            mk;
	int         i,
	            ii,
	    	    k;

	if(*idum < 0 || iff == 0)
	{							// Initialization.
		iff     = 1;
		mj      = labs(MSEED - labs(*idum));		// Initialize ma[55] using the seed idum and the
		mj     %= MBIG;					// large number MSEED.
		ma[55]  = mj;
		mk      = 1;
		for(i = 1; i <= 54; i++)
		{						// Now initialize the rest of the table,
			ii     = (21 * i) % 55;			// in a slightly random order,
			ma[ii] = mk;				// with numbers that are not especially random.
			mk     = mj - mk;
			if( mk < MZ )
				mk += MBIG;
			mj     = ma[ii];
		}
		for(k = 1; k <= 4; k++)				// We randomize them by “warming up the gener-
			for(i = 1; i <= 55; i++)
			{					// ator.”
				ma[i] -= ma[1 + (i + 30) % 55];
				if( ma[i] < MZ )
					ma[i] += MBIG;
			}
		inext   = 0;					// Prepare indices forour ﬁrst generated number.
		inextp  = 31;					// The constant 31 is special; see Knuth.
		*idum   = 1;
	}
	if( ++inext == 56 )					// Here is where we start, except on initialization.
		inext   = 1;					// Increment inext and inextp, wrapping around
	if( ++inextp == 56 )
		inextp  = 1;					// 56 to 1.
	mj = ma[inext] - ma[inextp];				// Generate a new random number subtractively.
	if( mj < MZ )
		mj     += MBIG;					// Be sure that it is in range.
	ma[inext] = mj; 					// Store it,
	return mj * FAC;					// and output the derived uniform deviate.
}

