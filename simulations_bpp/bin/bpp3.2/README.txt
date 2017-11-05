Notes by Ziheng Yang
23 May 2015

(1) To compile, try one of the following

   UNIX gcc compiler:
      cc -o bpp -O3 bpp.c tools.c -lm
      cc -o MCcoal -DSIMULATION bpp.c tools.c -lm

      gcc -o bpp -O3 bpp.c tools.c -lm
      gcc -o MCcoal -DSIMULATION bpp.c tools.c -lm

   INTEL icc compiler:
      icc -o bpp -fast bpp.c tools.c -lm
      icc -o MCcoal -DSIMULATION -fast bpp.c tools.c -lm

   MAC OSX intel:
      cc -o bpp -O3 bpp.c tools.c -lm
      cc -o MCcoal -DSIMULATION -O3 bpp.c tools.c -lm

   Windows Microsoft Visual C++:
      cl -O2 bpp.c tools.c
      cl -O2 -FeMCcoal.exe -DSIMULATION bpp.c tools.c

(2) To run an example analysis, try 

   cd examples

   ../bpp yu2001.bpp.ctl

   ../bpp ChenLi2001.bpp.ctl

   ../bpp lizard.bpp.ctl


Good luck.
