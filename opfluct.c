/* opfluct.c  by Szabo                              2021.08.23.  */
/* donation game on a two-layer square lattice                    */
/* with stochastic imitation of a neighbor within the same layer  */
/* MC simulations for studying spontaneous symmetry breaking      */
/* average C freqs., op fluc., and covariance vs r  */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#define adatnev "opfluct.d2t"      /* file to store data              */
#define rinit       0.0015          /* initial ratio                   */
#define rfin        0.0020          /* final ratio                     */
#define rN          5             /* number of different ratios      */
#define L           800            /* lattice size                    */
#define SIZE        (L*L)          /* number of sites                 */
#define NU_OF_NEIG  4              /* number of neighbors             */
#define MC_STEPS    200000          /* run-time in MCS unit            */
#define TERM        50000          /* thermalization in unit MCS      */
#define temp        0.1            /* noise level for imitation       */

typedef int       tomb1[2][L][L];
typedef double    tomb2[2][2];
typedef int       tomb3[L];
typedef long int  tomb4[3];
typedef double    tomb5[2];

tomb1 player_s;              /* matrix, defining player's strategy    */
tomb2 pm;                    /* payoff matrix                         */
tomb3 csok, nov;             /* neighbor's coordinate involving PBC   */
tomb4 Spopc;		         /* Cs in sublattices                     */

double d, r, prd, PX, PY, x0, xx, xy, xz, xs, z, temperature;
double rat, drat, cav0, cav1, covar, covc, opav, opX, r0, r1, r2;
int SX, SY, rser, strat1, ri, sr, sq, s1, s2, s3, s4, s5, lix, liy;
long int i, j, k, l, ii, iK, s, MCe, step, elestep;
int ix,jx,iy,jy,iz,jz,ip,im,jp,jm;

  FILE * tf;

long randl(long);   /* random number between 0 and num-1 */
double randd(void); /* double random number in  [0, 1)   */

/* Matsumoto Mersenne twister */

/* Period parameters */
#define N 624
#define M 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define LOWER_MASK 0x7fffffffUL /* least significant r bits */

static unsigned long mt[N]; /* the array for the state vector  */
static int mti=N+1; /* mti==N+1 means mt[N] is not initialized */

/* initializes mt[N] with a seed */
void init_genrand(unsigned long s)
{
    mt[0]= s & 0xffffffffUL;
    for (mti=1; mti<N; mti++) {
        mt[mti] =
	    (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti);
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array mt[].                        */
        /* 2002/01/09 modified by Makoto Matsumoto             */
        mt[mti] &= 0xffffffffUL;
        /* for >32 bit machines */
    }
}

/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
void init_by_array(unsigned long init_key[], unsigned long key_length)
// unsigned long init_key[], key_length;
{
    int i, j, k;
    init_genrand(19650218UL);
    i=1; j=0;
    k = (N>key_length ? N : key_length);
    for (; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1664525UL))
          + init_key[j] + j; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++; j++;
        if (i>=N) { mt[0] = mt[N-1]; i=1; }
        if (j>=key_length) j=0;
    }
    for (k=N-1; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1566083941UL))
          - i; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++;
        if (i>=N) { mt[0] = mt[N-1]; i=1; }
    }

    mt[0] = 0x80000000UL; /* MSB is 1; assuring non-zero initial array */
}

/* generates a random number on [0,0xffffffff]-interval */
unsigned long genrand_int32(void)
{
    unsigned long y;
    static unsigned long mag01[2]={0x0UL, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    if (mti >= N) { /* generate N words at one time */
        int kk;

        if (mti == N+1)   /* if init_genrand() has not been called, */
            init_genrand(5489UL); /* a default initial seed is used */

        for (kk=0;kk<N-M;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        for (;kk<N-1;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

        mti = 0;
    }

    y = mt[mti++];

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;
}

/* generates a random number on [0,0x7fffffff]-interval */
long genrand_int31(void)
{
    return (long)(genrand_int32()>>1);
}

/* generates a random number on [0,1]-real-interval */
double genrand_real1(void)
{
    return genrand_int32()*(1.0/4294967295.0);
    /* divided by 2^32-1 */
}

/* generates a random number on [0,1)-real-interval */
double genrand_real2(void)
{
    return genrand_int32()*(1.0/4294967296.0);
    /* divided by 2^32 */
}

/* generates a random number on (0,1)-real-interval */
double genrand_real3(void)
{
    return (((double)genrand_int32()) + 0.5)*(1.0/4294967296.0);
    /* divided by 2^32 */
}

/* generates a random number on [0,1) with 53-bit resolution*/
double genrand_res53(void)
{
    unsigned long a=genrand_int32()>>5, b=genrand_int32()>>6;
    return(a*67108864.0+b)*(1.0/9007199254740992.0);
}
/* These real versions are due to Isaku Wada, 2002/01/09 added */

long randl(long num)      /* random number between 0 and num-1 */
{
  return ((long) genrand_int31() % num);
}

double randd(void)      /* double random number in  [0, 1) */
{
  return((double) genrand_real2());
}

void main()
 {
  init_genrand(5617);    /* Init Matsuomoto random number generator */

  temperature = temp;
  rat = rinit;     drat = (rfin - rinit)/rN;
  for (ix=0; ix<L; ix++) { nov[ix] = ix+1; csok[ix]=ix-1;}
  csok[0]=L-1;  nov[L-1]=0;

 tf=fopen(adatnev,"w+t");
 fprintf(tf,"   r        rhoC0     rhoC1      op        opX       covar ");
 fprintf(tf," opfluct.c, ");
 fprintf(tf,"temp=%3.1f,",temp);
 fprintf(tf," samp=%d, term=%d,",MC_STEPS,TERM);
 fprintf(tf," N=%d\n",SIZE);
 fflush(tf);
 fclose(tf);

 MCe = 2*SIZE;

   for (ix=0; ix<L; ix++) {
    for (jx=0; jx<L; jx++) { strat1 = 1; 
	                         r = randd();
			                 if ( r<0.8) { strat1=0; }
						     player_s[0][ix][jx] = strat1;
							 strat1 = 1; 
	                         r = randd();
			                 if ( r<0.2) { strat1=0; }
						     player_s[1][ix][jx] = strat1;
						  }}   /*   prepared initial state    */

 for (rser=0; rser<=rN; rser++)
  {
    pm[0][0] = 1.0 - rat;    pm[0][1] = - rat;
    pm[1][0] = 1.0;          pm[1][1] = 0.0;

    cav0 = 0.0;   cav1 = 0.0;
	opav = 0.0;   covar = 0.0;    opX = 0.0;

  for(step=1; step<=(MC_STEPS+TERM); step++)
   {
    for(elestep=0; elestep<MCe; elestep++)
     {
      ix = (int) randl(L);
      jx = (int) randl(L);
      lix = (int) randl(2);
      liy = 1 - lix;
      SX = player_s[lix][ix][jx];
      ip = nov[ix]; im = csok[ix];
      jp = nov[jx]; jm = csok[jx];
      ri = (int) randl(4);
	switch (ri) {
		case 0:   iy = ix;   jy = jp;   break;
		case 1:   iy = im;   jy = jx;   break;
		case 2:   iy = ix;   jy = jm;   break;
		case 3:   iy = ip;   jy = jx;   break;
				}
	  SY = player_s[lix][iy][jy];
	if (SX != SY) {
			s = player_s[liy][ix][jx];                  PX=pm[SX][s];
			iz=ip;  jz=jx;    s=player_s[liy][iz][jz];  PX+=pm[SX][s];
			iz=ix;  jz=jp;    s=player_s[liy][iz][jz];  PX+=pm[SX][s];
			iz=im;  jz=jx;    s=player_s[liy][iz][jz];  PX+=pm[SX][s];
			iz=ix;  jz=jm;    s=player_s[liy][iz][jz];  PX+=pm[SX][s];

			ip = nov[iy];      jp = nov[jy];
			im = csok[iy];     jm = csok[jy];

			s = player_s[liy][iy][jy];                 PY=pm[SY][s];
			iz=ip;  jz=jy;   s=player_s[liy][iz][jz];  PY+=pm[SY][s];
			iz=iy;  jz=jp;   s=player_s[liy][iz][jz];  PY+=pm[SY][s];
			iz=im;  jz=jy;   s=player_s[liy][iz][jz];  PY+=pm[SY][s];
			iz=iy;  jz=jm;   s=player_s[liy][iz][jz];  PY+=pm[SY][s];

			prd = 1.0 / (1 + exp(-(PY - PX)/temperature));
			r = randd();
			if ( r<prd) { player_s[lix][ix][jx] = SY;  }
						}

     } /* end of elementary MC step  */

    if(step>TERM) { Spopc[0]=0; Spopc[1]=0; Spopc[2]=0; 

                 for (ix=0; ix<L; ix++) {
                 for (jx=0; jx<L; jx++) {
					sr = player_s[0][ix][jx];
					sq = player_s[1][ix][jx];
                    Spopc[0] += (1 - sr); 
                    Spopc[1] += (1 - sq); 
					Spopc[2] += ((1-sr)*(1-sq));
					}}
				r0=(double) Spopc[0]; 
				r0/=SIZE; cav0 += r0;
                		r1=(double) Spopc[1];
				r1/=SIZE; cav1 += r1;
				opav += (r0-r1);  r2= (r0-r1)*(r0-r1);
				opX += r2;
				r0 =(double) Spopc[2]; 
				r0 /= SIZE; covar += r0;
	   }

   } /*  end of MC steps */

	cav0 /= MC_STEPS;       cav1 /= MC_STEPS; 
	opav /= MC_STEPS; 	opX /= MC_STEPS;
	opX -=(opav*opav);	opX *= SIZE;
	covar /= MC_STEPS;       covar -=(cav0*cav1);

      tf=fopen(adatnev,"a+t");
      fprintf(tf,"%9.6f",rat);
      fprintf(tf,"%10.6f%10.6f",cav0,cav1);
      fprintf(tf,"%10.6f%10.2f%10.6f\n",opav,opX,covar);
      fflush(tf);
      fclose(tf);

   rat += drat;
  } /* rser */

   Vege: ;
 }

/* -------------------------------------------------------- */
