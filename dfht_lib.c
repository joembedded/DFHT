/***************************************************************
*
* DFHT_LIB.C
*
* (C)JoEmbedded.de
*
* Discrete Fast Hartley Transform - Digital Signal Processing
* -----------------------------------------------------------
*
* Recipe:
* - Fill an array with periodically equidistant samples
* - As an option weight the data by a 'window' function to 
*   reduce the Leackage Effect (omitted here)
* - Do the transformation 'dfht()'
* - Get the power spectrum with  'power()'
* - Make a decision.
*
* This software is the 'real' variant of the FFT. Because the FFT is
* based on complex numbers, it normally is ideal for stero signals.
* For mono/single signals the DFHT is an easy2use alternative
* (although e.g. on ARM the ARM-Lib-FFT is faster, but less portable).
* The same routine can be used for analysis and synthesizes!
*
* Infos:     Oxford University Press:
*            Ronald N. Bracewell: The Hartley Transformation
*            https://en.wikipedia.org/wiki/Hartley_transform
*
* Discrete  Fourier Transform, N Elements:
* ----------------------------------------
*
* F(v) = 1/N {Sum t=0 to N-1: f(t) * exp(-i * 2*PI * v * t / N)}
* f(t) = {Sum v=0 to N-1: F(v) * exp( i * 2*PI * v * t / N)}
*
* Discrete Hartley Transform, N Elements:
* -----------------------------------------
*
* H(v) = 1/N {Sum t=0 to N-1: f(t) * cas( 2*PI * v * t /N )}
* f(t) = {Sum v=0 to N-1: H(v) * cas ( 2*PI * v * t / N)}
* Note: cas(x) = cos(x) + sin(x)
*
* Power Spectrum:
* for v=0: P(0) = ( H(0)^2 )  Rem:( H(N-0) == H(0) ) // DC-Component
* for (v=1..(N/2 - 1)): P(v) = ( H(v)^2 + H(N-v)^2 ) / 2
* (for v=N/2: P(N/2) = ( H(N/2)^2 )  // unused for Power Spectrum)
*
* Angle:
* for v=0: P(0): (DC component)  Rem:( H(N-0) == H(0) )
* for (v=1..(N/2 - 1)): j(v) = arctan( ( H(v) - H(N-v) ) / ( H(v) + H(N-v) ) )
*
* Mapping the  Hartley Transform <-> Fourier Transform (see Demo 
* to get sin-/cos-Coefficients):
* F(v)=E(v)+ i * O(v) (complex!)
* E(v) = H(v) + H(N - v), O(v) = H(v) - H(N - v)
*
* About Speed:
* - Standard Desktop PC (ca. 2019): 8k DFHT in < 1 msec
* - nRF52@64MHz: 4k DFHT ca. 30 msec ('Release') 
*   (Remark: ARM-Lib-FFT: 4k FFT ca. 10 msec)
* 
******************************************************************/
#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS // For VisualStudio
#define _USE_MATH_DEFINES // for C
#endif

#include <math.h>                               /* Math */
#include <stdio.h>                              /* Standard I/O    */
#include <stdlib.h>		   // exit
#include <stdint.h>

#include "dfht_lib.h"


/* Support */
#if DATA_SIZE == 8
#define DATA_BITS 3
#elif DATA_SIZE == 16
#define DATA_BITS 4
#elif DATA_SIZE == 32
#define DATA_BITS 5
#elif DATA_SIZE == 64
#define DATA_BITS 6
#elif DATA_SIZE == 128
#define DATA_BITS 7
#elif DATA_SIZE == 256
#define DATA_BITS 8
#elif DATA_SIZE == 512
#define DATA_BITS 9
#elif DATA_SIZE == 1024
#define DATA_BITS 10
#elif DATA_SIZE == 2048
#define DATA_BITS 11
#elif DATA_SIZE == 4096
#define DATA_BITS 12
#elif DATA_SIZE == 8192
#define DATA_BITS  13
#elif DATA_SIZE == 16384
#define DATA_BITS 14
#elif DATA_SIZE == 32768
#define DATA_BITS 15
#elif DATA_SIZE == 65536
#define DATA_BITS 16
#else
#error "DATA_SIZE Error"
#endif


// From DATA_SIZE derived constants...
#define NEL 		DATA_SIZE
#define NEL_d2 		(DATA_SIZE/2)
#define NEL_d4 		(DATA_SIZE/4)
#define NEL_d4m1 	(DATA_SIZE/4-1)

static uint8_t init = 0;

/* Modul global variablen */

/***********************************************************
* Some data
***********************************************************/
static float work_buf[DATA_SIZE];                        /* work buffer */
static float sinfunc[3 * NEL_d4 +1 ];

/***********************************************************
* init_dfht(): Initialising
***********************************************************/
void static init_dfht(void) {
	uint16_t i;

	if (init) return;              /* Do it only once */
	init = 1;
	/* init the sine table */
	for (i = 0; i < 3 * NEL_d4 + 1; i++) sinfunc[i] = (float)sin((i * M_TWO_PI) / (float)NEL);

}
/*****************************************************************
* bit_reverse(): permute buffer 'perfect shuffle'
*****************************************************************/
void static bit_reverse(float* pfeld) {             /* pfeld points to the datea */
	uint16_t NEL_cnt;                                 /* Counter */
	uint16_t ind_s;                                 /* Sourceindex */
	uint16_t ind_r;                                 /* Reversed Index */
	uint8_t bits;                                   /* in bits */
	float swap_tmp;

	for (NEL_cnt = 0; NEL_cnt < NEL; NEL_cnt++) {             /* all elements */
		ind_s = NEL_cnt;
		ind_r = 0;
		for (bits = 0; bits < DATA_BITS; bits++) {       /* all bits */
			ind_r <<= 1;                          /* assume 0 */
			if (ind_s & 1) ind_r++;              /* transfer bit */
			ind_s >>= 1;                          /* next */
		}
		if (ind_r > NEL_cnt) {                        /* Swap only if target- */
			swap_tmp = pfeld[ind_r];              /* index > the source! */
			pfeld[ind_r] = pfeld[NEL_cnt];
			pfeld[NEL_cnt] = swap_tmp;
		}
	}
}
/***********************************************************************
* stage_1(), stage_2(), stage_x(): ransformation-Stages
***********************************************************************/
void stage_1(float* sfeld, float* dfeld) {       /* Stage 1        */
	uint16_t NEL_cnt;
	float f0, f1;

	/* Important Stage  1: Choose INPLACE or TRANSPORT because of the Buffer */
	for (NEL_cnt = 0; NEL_cnt < NEL; NEL_cnt += 2) {
		f0 = *sfeld++;                          /* Stage 1: add only */
		f1 = *sfeld++;
		*dfeld++ = f0 + f1;
		*dfeld++ = f0 - f1;
	}
}
void stage_2(float* sfeld, float* dfeld) {       /* Stage 2        */
	uint16_t NEL_cnt;
	float f0, f1, f2, f3;
	for (NEL_cnt = 0; NEL_cnt < NEL; NEL_cnt += 4) {
		f0 = *sfeld++;                          /* Stufe 2: add only too */
		f1 = *sfeld++;
		f2 = *sfeld++;
		f3 = *sfeld++;
		*dfeld++ = f0 + f2;
		*dfeld++ = f1 + f3;
		*dfeld++ = f0 - f2;
		*dfeld++ = f1 - f3;
	}
}

void stage_x(float* sfeld, float* dfeld, uint8_t ist_x) { /* Stage 3,4,5... */
	uint16_t stub, stub_d2;                         /* Step etc. */
	uint16_t n_cnt;                                 /* block count */
	uint16_t i_cnt;                                 /* matrix elements */
	uint16_t tri_s;                                 /* Trigonom step */

	float dt;                                   /* Temp pointer for */
	float* ps1, * pc1;                            /* faster Arrays */
	float* pq1, * pqs, * pqc;
	float* pz1, * pz2;

	stub = 1 << ist_x;                             /* Stage step size */
	stub_d2 = stub >> 1;                            /* Half of it */
	tri_s = NEL / stub;                               /* per Block one repitition */

	/* This code monster does the complete calculation of on block of the stage */
	for (n_cnt = 0; n_cnt < NEL; n_cnt += stub) {

		ps1 = sinfunc;
		pc1 = sinfunc + NEL_d4;

		pz1 = dfeld;                              /* 1-1-Dest. */
        pz2 = dfeld + stub_d2;                      /* 1-2-Dest.*/
		pq1 = sfeld;                              /* 1-Source */
		pqc = sfeld + stub_d2;                      /* Cosine source */
		pqs = sfeld + stub;                         /* Sine source */

		*pz1++ = *pq1 + *pqc;                   /* Exeption 1 */
		*pz2++ = *pq1++ - *pqc++;

		for (i_cnt = 1; i_cnt < stub_d2; i_cnt++) {

			pc1 += tri_s;
			ps1 += tri_s;

			/* Retrograd-Indexing */
			dt = (*pqc++ * *pc1) + (*--pqs * *ps1);

			*pz1++ = *pq1 + dt;
			*pz2++ = *pq1++ - dt;
		}
		sfeld += stub;
		dfeld += stub;
	}

}
/***********************************************************************
*
* void dfht_do(pdaten, stellen, richtung): Die DFHT
*         This DFHT is done INPLACE. Nevertheless we need a temporoary
*         Buffer., because a 'real' inplace solution would have taken
*         a little more code...
*
***********************************************************************/
void dfht_do(float* pdaten, uint8_t richtung) {
	float* z0, * z1, * zt;                  /* Change: Buffer, Data */
	uint8_t i;
	uint16_t j;

#ifdef _VERBOSE
	printf("*** DFHT Initalise: %d Bits (Stages)... ***\n", DATA_BITS);
#endif
	init_dfht();

#ifdef _VERBOSE
	printf("*** DFHT Permute, %d Elements... ***\n", NEL);
#endif
	bit_reverse(pdaten);                        /* permute data */

#ifdef _VERBOSE
	puts("*** DFHT Stage 1... ***");
#endif
#if (DATA_BITS & 1)                           /* End: Data in pdaten[] */
	stage_1(pdaten, pdaten);                /* Odd: stage_1: Inplace */
	z0 = work_buf;
	z1 = pdaten;
#else
	stage_1(pdaten, work_buf);              /* Even: stage_1: Transport */
	z1 = work_buf;
	z0 = pdaten;
#endif
#ifdef _VERBOSE
	puts("*** DFHT Stage 2... ***");
#endif
	stage_2(z1, z0);                            /* Else: Toggle! */
	for (i = 3; i <= DATA_BITS; i++) {
#ifdef _VERBOSE
		printf("*** DFHT Stage %d... ***\n", i);
#endif
		stage_x(z0, z1, i);
		zt = z0; z0 = z1; z1 = zt;
	}
	if (!richtung) {
#ifdef _VERBOSE
		printf("*** DFHT %d normalise elements ... ***\n", NEL);
#endif
		for (j = 0; j < NEL; pdaten[j++] /= NEL);
	}
#ifdef _VERBOSE
	puts("*** DFHT Done! ***");
#endif
}
/***********************************************************************
*
* dfht_power(): Powerspectrum calculation
*          float sfeld[]: The coefficients of the Hartley-Transf.
*          float dfeld[]: Here the result: the spectrum (!!! half size of sfeld !!!)
*          Return: Value of the largest element (use for display, ...)
*
************************************************************************/

float dfht_power(float* sfeld, float* dfeld) {
	float maxp;
	float p;
	uint16_t i;

#ifdef _VERBOSE
	printf("*** DFHT %d Elements powerspectrum... ***\n", NEL);
#endif
	dfeld[0] = 2 * sfeld[0] * sfeld[0];             /* Exeption */
	maxp = dfeld[0];                              /* Maximum init. */
	for (i = 1; i < NEL_d2; i++) {
		p = (sfeld[i] * sfeld[i]) + (sfeld[NEL - i] * sfeld[NEL - i]);
		dfeld[i] = p;
		if (p > maxp) maxp = p;
	}
#ifdef _VERBOSE
	puts("*** DFHT Powerspectrum done! ***");
#endif
	return maxp;
}
/**** END ****/


