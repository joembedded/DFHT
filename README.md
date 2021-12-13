# DFHT - Discrete Fast Hartley Transform
_A portable and easy2use alternative for a FFT_

This software is the 'real' variant of the FFT. Because the FFT is
based on complex numbers, it normally is ideal for stero signals.
For mono/single signals the DFHT is an easy2use alternative
(although e.g. on ARM the ARM-Lib-FFT is faster, but less portable).
The same routine can be used for analysis and synthesizes!

## About Speed: ##
- Standard Desktop PC (ca. 2019): 8k DFHT in < 1 msec
- nRF52@64MHz: 4k DFHT ca. 30 msec ('Release') 
  (Remark: ARM-Lib-FFT: 4k FFT ca. 10 msec)

## Installation/Test ##
- Tested with Visual Studio (Windows) and SES (GCC)

## Links ##
- Wikipedia: https://en.wikipedia.org/wiki/Hartley_transform

***
