/******************************************************************
*
* DFHT.C - Demo of DFHT_LIB.C
*
* (C)JoEmbedded.de
*
******************************************************************/
#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS // For VisualStudio
#define _USE_MATH_DEFINES // for C
#endif


#include <math.h>                               /* Math */
#include <stdio.h>                              /* Standard I/O    */
#include <stdint.h>

#include "dfht_lib.h"   // DATA

/***********************************************************************
* Some DEMO Variables...
************************************************************************/
#define LINES 20	// No of display lines
#define COLUMS 64

static float daten[DATA_SIZE];  // This field holds the data
static float original_daten[DATA_SIZE];  // (Original for verification)
static float power_buf[DATA_SIZE / 2];			/* A buffer for the spectrum */
static char display[LINES][COLUMS + 1]; // A 'virtual' screen

/***********************************************************************
* clr_disp(): Clear 'virtual' screen
************************************************************************/
void clr_disp(void) {
	unsigned int i, j;
	char* pc = &display[0][0];
	for (i = 0; i < LINES; i++) {
		for (j = 0; j < COLUMS; j++) {
			*pc++ = '.';
		}
		*pc++ = '\n';
	}
}

/***********************************************************************
* show_disp(): Move virtual screen to the UART
************************************************************************/
void show_disp(void) {
	unsigned int i;
	char* pc = &display[0][0];
	for (i = 0; i < LINES * (COLUMS + 1); i++) {
		putchar(*pc++);
	}
}
/***********************************************************************
* data_disp(): Fill display for power
************************************************************************/
void data_disp(void) {
	unsigned int i, j;
	for (i = 0; i < DATA_SIZE; i++) {
		j = LINES / 2 - (int)daten[i];
		if (j < LINES) display[j][(i * COLUMS) / DATA_SIZE] = '#';
	}
}

/***********************************************************************
* data_disp(): Fill display for power (only half screen used)
************************************************************************/
void power_disp(float pmax) {
	unsigned int i;
	int yval;
	float fval;
	if (pmax == 0.0) return;	//  /0
	for (i = 0; i < DATA_SIZE / 2; i++) {
		fval = power_buf[i] / pmax;
		yval = LINES - ((int)(fval * (float)LINES));
		if (yval < 0) yval = 0;
		else if (yval >= LINES) yval = LINES - 1;
		display[yval][(i * COLUMS) / DATA_SIZE] = '#';
	}
}

/****************************************************
* ###  M A I N ###
****************************************************/

void main(void) {
	uint16_t err;
	uint16_t i;
	float pmax;
	for (err = 0; err < 5; err++) {	// 5 Runs
		printf("*** DFHT-Demo, Run:%d ****\n",err);
		// Fill field with sample data - 'err' will demonstrate the Leackage Effect 
		for (i = 0; i < DATA_SIZE; i++) {
			float fvt = (float)(0.4 * LINES * sin((2.0 + (float)err / 5.0) * M_TWO_PI * (float)i / (float)DATA_SIZE));
			// fvt += (float) ( (rand() & 255) - 128)/100.0;	// Optinally with Noise
			daten[i] = original_daten[i] = fvt;
		}
#if DATA_SIZE>=32
		printf("Original data with error %d:\n", err);
		clr_disp();
		data_disp();
		show_disp();
		puts("<NL>");
		(void)getchar();
#endif
#if DATA_SIZE<=64
		// Als Zahlen Input
		printf("Original data with error %d:\n", err);
		for (i = 0; i < DATA_SIZE; i++) printf("%d: %f\n", i, daten[i]);
		puts("<NL>");
		(void)getchar();
#endif	
		puts("Wait...");

		dfht_do(daten, 0);	// Do analysis, 0: analysis (1: synthesizes)

		pmax = dfht_power(daten, power_buf);

#if DATA_SIZE>=32
		printf("Power spectrum for error %d (max: %f):\n", err, pmax);
		clr_disp();
		power_disp(pmax);
		show_disp();
		puts("<NL>");
		(void)getchar();
		puts("Wait...");
#endif
		// Als Zahlen Output

#if DATA_SIZE <= 64

#if 0 /* (Enable if interesting) */
		printf("Hartley Coeffs with error %d:\n", err);
		for (i = 0; i < DATA_SIZE; i++) printf("H%d: %f\n", i, daten[i]);

		printf("Fourier Coeffs with error %d:\n", err);
		for (i = 0; i <= DATA_SIZE/2; i++) {
			float k_cos, k_sin;
			if (i == DATA_SIZE / 2) {
				k_cos = daten[DATA_SIZE / 2]; // cos(): +1;-1;+1;...
				k_sin = 0;					  // sin(): 0;0;0...
			}
			else if (i) {
				k_cos = (daten[i] + daten[(DATA_SIZE - i)]);
				k_sin = (daten[i] - daten[(DATA_SIZE - i)]);
			}
			else { // DC Component (i==0)
				k_cos = daten[0]; // cos(): 1;1;1;...
				k_sin = 0;		  // sin(): 0;0;0
			}
			printf("F%d: Kcos: %f Ksin: %f\n", i,k_cos,k_sin);
		}
		puts("<NL>");
		(void)getchar();
#endif // Coeffs

		printf("Power spectrum for error %d:\n", err);
		for (i = 0; i < DATA_SIZE / 2; i++) {
			printf("%d: P:%f\n", i, power_buf[i]);
		}

		printf("Re-SynOriginal (Classical Slow Fourier Syntheseis) data with error %d:\n", err);
		for (unsigned int n = 0; n < DATA_SIZE; n++) {
			float sum = 0;
			for (i = 0; i <= DATA_SIZE / 2; i++) {  
				float frq = (float)n * M_TWO_PI * (float)i / (float)DATA_SIZE;
				float k_cos, k_sin;
				if (i == DATA_SIZE / 2) { 
					k_cos = daten[DATA_SIZE / 2]; // cos(): +1;-1;+1;...
					k_sin = 0;					  // sin(): 0;0;0...
				}
				else if (i) {
					k_cos = (daten[i] + daten[(DATA_SIZE - i)]);
					k_sin = (daten[i] - daten[(DATA_SIZE - i)]);
				}
				else { // DC Component (i==0)
					k_cos = daten[0]; // cos(): 1;1;1;...
					k_sin = 0;		  // sin(): 0;0;0
				}
				sum += k_cos * cos(frq) + k_sin * sin(frq);
			}
			printf("Slow F:%d: %f (delta:%f)\n", n, sum, sum-original_daten[n]);
		}
		puts("<NL>");
		(void)getchar();
#endif	

		dfht_do(daten, 1);	// Do analysis, 0: analysis (1: synthesizes)

#if DATA_SIZE>=32
		printf("Re-SynOriginal data with error %d:\n", err);
		clr_disp();
		data_disp();
		show_disp();
		puts("<NL>");
		(void)getchar();
#endif
#if DATA_SIZE<=64
		// Als Zahlen Syn
		printf("Re-SynOriginal data with error %d:\n", err);
		for (i = 0; i < DATA_SIZE; i++) printf("%d: %f (delta:%f)\n", i, daten[i], daten[i]-original_daten[i]);
		puts("<NL>");
		(void)getchar();
#endif	
	}
}

/**** END MAIN TEST ****/


