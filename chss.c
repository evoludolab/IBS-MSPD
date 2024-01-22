/* chss.c  by Szabo                                  2021.02.22.  */
/* donation game on a two-layer square lattice                    */
/* with stochastic imitation of a neighbor within the same layer  */
/* MC simulations for studying pattern evolution in a small box   */
/* start from a random initial state                              */
/* snapshots are printed on an eps file                           */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#define epsnev "chss022.eps"      /* eps file to visualize snapshot  */
#define r1          0.0022        /* cost to benefit ratio           */
#define TERM        100           /* thermalization in MCS unit      */
#define L           800           /* linear size of lattice          */
#define Lb          200           /* linear size of box              */
#define Lx          20            /* horizontal location of box      */
#define Ly          20            /* vertical location of box        */
#define SIZE        (L*L)         /* number of sites                 */
#define NU_OF_NEIG  4             /* number of neighbors             */
#define temp        0.1           /* noise level for imitation       */

typedef char      tomb1[L][L][2];
typedef double    tomb2[2][2];
typedef int       tomb3[L];
typedef long int  tomb4[2];
typedef double    tomb5[4];

tomb1 player_s;              /* matrix, defining player's strategy    */
tomb2 pm;                    /* payoff matrix                         */
tomb3 csok, nov;             /* neighbor's coordinate involving PBC   */
tomb5 podo;  		         /* portion of domains within the box     */

double d, r, prd, PX, PY, z, temperature;
double rat, r0;
int SX, SY, tser, strat1, ri, s1, s2, s4, lay1, lay2, Lb2, sk;
long int i, j, k, l, ii, jj, iK, s, MCe, step, elestep, Nsiav;
int ix,jx,iy,jy,iz,jz,ip,im,jp,jm,tMCS;

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
  init_genrand(5391);    /* Init Matsuomoto random number generator */

  temperature = temp;
  rat = r1;
  for (ix=0; ix<L; ix++) { nov[ix] = ix+1; csok[ix]=ix-1;}
  csok[0]=L-1;  nov[L-1]=0;
  MCe = 2*SIZE;

   for (ix=0; ix<L; ix++) {
    for (jx=0; jx<L; jx++) { strat1 = 0;
			     r = randd();
			     if ( r<0.15) { strat1 = 1 - strat1;  }
 			     strat1 = randl(2); 
				 player_s[ix][jx][0] = strat1; 
			     r = randd();
			     if ( r<0.15) { strat1 = 1 - strat1;  }
 			     strat1 = randl(2);
                 player_s[ix][jx][1] = strat1; 
			  }}             /*  prepared random initial state    */

    pm[0][0] = 1.0 - rat;    pm[0][1] = - rat;
    pm[1][0] = 1.0;          pm[1][1] = 0.0;

  for(step=0; step<TERM; step++)
   {
    for(elestep=0; elestep<MCe; elestep++)
     {
      ix = (int) randl(L);
      jx = (int) randl(L);
      lay1 = (int) randl(2);
      lay2 = 1 - lay1;
      SX = player_s[ix][jx][lay1];
      ip = nov[ix]; im = csok[ix];
      jp = nov[jx]; jm = csok[jx];
      ri = (int) randl(4);
	switch (ri) {
		case 0:   iy = ix;   jy = jp;   break;
		case 1:   iy = im;   jy = jx;   break;
		case 2:   iy = ix;   jy = jm;   break;
		case 3:   iy = ip;   jy = jx;   break;
				}
	  SY = player_s[iy][jy][lay1];
	if (SX != SY) {
			s = player_s[ix][jx][lay2];                  PX=pm[SX][s];
			iz=ip;  jz=jx;    s=player_s[iz][jz][lay2];  PX+=pm[SX][s];
			iz=ix;  jz=jp;    s=player_s[iz][jz][lay2];  PX+=pm[SX][s];
			iz=im;  jz=jx;    s=player_s[iz][jz][lay2];  PX+=pm[SX][s];
			iz=ix;  jz=jm;    s=player_s[iz][jz][lay2];  PX+=pm[SX][s];

			ip = nov[iy];      jp = nov[jy];
			im = csok[iy];     jm = csok[jy];

			s = player_s[iy][jy][lay2];                 PY =pm[SY][s];
			iz=ip;  jz=jy;   s=player_s[iz][jz][lay2];  PY+=pm[SY][s];
			iz=iy;  jz=jp;   s=player_s[iz][jz][lay2];  PY+=pm[SY][s];
			iz=im;  jz=jy;   s=player_s[iz][jz][lay2];  PY+=pm[SY][s];
			iz=iy;  jz=jm;   s=player_s[iz][jz][lay2];  PY+=pm[SY][s];

			prd = 1.0 / (1 + exp(-(PY - PX)/temperature));
			r = randd();
			if ( r<prd) { player_s[ix][jx][lay1] = SY;  }
						}
     } /* end of elementary MC step  */
	}  /* end of thermalization      */

                                            /* creation of EPS file */
 tf=fopen(epsnev,"w+t");
 fprintf(tf,"%c!PS-Adobe-2.0\n",'%');  
 fprintf(tf,"%c%cTitle: snapshot for two-layer donation game, L=%d, K=%3.1f, t=%d MCS\n",'%','%',L,temperature,TERM);
 fprintf(tf,"%c%c Lx=%d, Ly=%d, Lb=%d \n",'%','%',Lx,Ly,Lb);  

 fprintf(tf,"%c%cCreator: ch2lsst.c\n",'%','%');  
 fprintf(tf,"%c%cBoundingBox: 15 15 505 505\n",'%','%');
 fprintf(tf,"%c%cOrientation: Portrait\n",'%','%');  
 fprintf(tf,"%c%cEndComments\n",'%','%');  
 fprintf(tf,"/M {rmoveto} bind def\n");  
 fprintf(tf,"/L {lineto} bind def\n");  
 fprintf(tf,"/R {rlineto} bind def\n");  
 fprintf(tf,"/S {moveto} bind def\n");  
 fprintf(tf,"/nr {grestore 0 2 translate gsave}  def\n");  
 fprintf(tf,"/ca {1.0 1.0 1.0 setrgbcolor} bind def\n");  
 fprintf(tf,"/cc {1.0 0.0 0.0 setrgbcolor} bind def\n");  
 fprintf(tf,"/cb {0.0 0.0 1.0 setrgbcolor} bind def\n");
 fprintf(tf,"/cd {0.0 0.0 0.0 setrgbcolor} bind def\n");  
 fprintf(tf,"/a {ca newpath 0 0 S 2 0 R 0 2 R -2 0 R closepath fill 2 0 translate} def\n");  
 fprintf(tf,"/b {cb newpath 0 0 S 2 0 R 0 2 R -2 0 R closepath fill 2 0 translate} def\n");  
 fprintf(tf,"/c {cc newpath 0 0 S 2 0 R 0 2 R -2 0 R closepath fill 2 0 translate} def\n"); 
 fprintf(tf,"/d {cd newpath 0 0 S 2 0 R 0 2 R -2 0 R closepath fill 2 0 translate} def\n");  
 fprintf(tf,"20 20 translate\n");  
 fprintf(tf,"1 1 scale\n");  
 fprintf(tf,"gsave\n");  

   sk = Lb / 40;  Lb2 = 2*Lb;
   for (j = Ly; j < (Ly+Lb); j++) {
    for (ii = 0; ii < sk; ii++)    {
     for (jj = 0; jj < 40; jj++)    {
     i = Lx+40*ii + jj;
     s1 = 2*player_s[i][j][0]+player_s[i][j][1];
      switch (s1) {
        case 0:  fprintf(tf,"a "); break;
        case 1:  fprintf(tf,"b "); break;
        case 2:  fprintf(tf,"c "); break;
		case 3:  fprintf(tf,"d "); break;
                  } /* case s1 */
		                   }
     fprintf(tf,"\n");               }
     fprintf(tf,"nr\n"); }
   fprintf(tf,"0 %d translate\n",(-Lb2));
   fprintf(tf,"0.5 setlinewidth 0 setgray newpath\n");
   fprintf(tf,"0 0 S %d 0 R 0 %d R %d 0 R closepath stroke\n",Lb2,Lb2,(-Lb2));
   fprintf(tf,"showpage\n");  
   fprintf(tf,"%c%cTrailer\n",'%','%');  
   fflush(tf);
   fclose(tf);

	
 }
/* -------------------------------------------------------- */
