          seed = -123

       seqfile = test3s.txt
      Imapfile = Imap.3s.txt
       outfile = out.txt
      mcmcfile = mcmc.out

* speciesdelimitation = 0 * fixed species tree
* speciesdelimitation = 1 0 2    * speciesdelimitation algorithm0(e)
 speciesdelimitation = 1 1 4 1 * speciesdelimitation algorithm1(a m)
         speciestree = 0 * species-tree by NNI
*        speciestree = 1  0.4 0.2 0.1   * speciestree pSlider ExpandRatio ShrinkRatio

   speciesmodelprior = 1  * 0: uniform LH; 1:uniform rooted trees; 2: uniformSLH; 3: uniformSRooted 

  species&tree = 3  A  B  C
                    2  5  3
                   ((A, B), C);

       usedata = 0    * 0: no data (prior); 1:seq like
         nloci = 2  *100    * number of data sets in seqfile

     cleandata = 0    * remove sites with ambiguity data (1:yes, 0:no)?

    thetaprior = 5 5    # gamma(a, b) for theta
      tauprior = 4 8 1  # gamma(a, b) for root tau & Dirichlet(a) for other tau's

*     locusrate = 0 2.0   # (0: No variation, 1: estimate, 2: from file) & a_Dirichlet (if 1)
*      heredity = 1 4 4   # (0: No variation, 1: estimate, 2: from file) & a_gamma b_gamma (if 1)
* sequenceerror = 0 0 0 0 0 : 0.05 1   # sequencing errors: gamma(a, b) prior

       finetune = 1: 1 1 0.1 0.1 0.2 0.5 1.0  # finetune for GBtj, GBspr, theta, tau, mix, locusrate, seqerr

         print = 1 0 0 0   * MCMC samples, locusrate, heredityscalars Genetrees
        burnin = 10000
      sampfreq = 10
       nsample = 100000

* This is the true tree from which the data test3s.txt are simulated:
* ((A #0.01, B #0.02) :0.01 #0.05, C #0.05) :0.02 #0.04;
*** Note: Make your window wider (144 columns) before running the program.

