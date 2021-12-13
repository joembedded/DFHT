/******************************************************************
*
* DFHT_LIB.H
*
***********************************************************/

#ifndef M_TWO_PI
#define M_TWO_PI		(M_PI*2)		// not ANSII Standard
#endif

/***********************************************************
* Define the size for the transformation
* Must be a power of 2 and size 8..65536
***********************************************************/

//#define DATA_SIZE  8192
#define DATA_SIZE  64


//#define _VERBOSE // If defined: printf(..)

void dfht_do(float* pdaten, uint8_t richtung);
float dfht_power(float* sfeld, float* dfeld);

// END