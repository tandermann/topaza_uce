          seed = 1

       seqfile = test4s.txt
      Imapfile = Imap.4s.txt
       outfile = out.txt
      mcmcfile = mcmc.txt

* speciesdelimitation = 0 * fixed species tree
*  speciesdelimitation = 1 0 2  * speciesdelimitation algorithm0 and finetune(e)
*  speciesdelimitation = 1 1 2 1 * speciesdelimitation algorithm1 finetune (a m)
           speciestree = 1 * species-tree by NNI; 0 * species tree fixed
*          speciestree = 1  0.4 0.2 0.1   * speciestree pSlider ExpandRatio ShrinkRatio

speciesmodelprior = 1  * 0: uniform LH; 1:uniform rooted trees; 2: uniformSLH; 3: uniformSRooted

  species&tree = 4  A  B  C  D
                    3  3  3  3
                  ((A, B), (C, D));

       usedata = 1    * 0: no data (prior); 1:seq like
         nloci = 2    * number of data sets in seqfile

     cleandata = 0    * remove sites with ambiguity data (1:yes, 0:no)?

    thetaprior = 2 200    # gamma(a, b) for theta
      tauprior = 4 200 1  # gamma(a, b) for root tau & Dirichlet(a) for other tau's

*     locusrate = 0 2.0   # (0: No variation, 1: estimate, 2: from file) & a_Dirichlet (if 1)
*      heredity = 0 4 4   # (0: No variation, 1: estimate, 2: from file) & a_gamma b_gamma (if 1)
* sequenceerror = 0 0 0 0 0 : 0.05 1   # sequencing errors: gamma(a, b) prior

       finetune = 1: 8.42859  0.00499  0.01654  0.01232  0.60284  # finetune for GBtj, GBspr, theta, tau, mix, locusrate, seqerr

         print = 1 0 0 0   * MCMC samples, locusrate, heredityscalars Genetrees
        burnin = 4000
      sampfreq = 1
       nsample = 100000
