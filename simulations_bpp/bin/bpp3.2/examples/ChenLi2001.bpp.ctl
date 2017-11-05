          seed =  -1

       seqfile = ChenLi2001.txt
      Imapfile = ChenLi2001.Imap.txt
       outfile = out.txt
      mcmcfile = mcmc.txt

   speciesdelimitation = 0 * fixed species tree
*  speciesdelimitation = 1 0 2    * speciesdelimitation algorithm0 and finetune(e)
*  speciesdelimitation = 1 1 2 1  * speciesdelimitation algorithm1 finetune (a m)

  species&tree = 4  H  C  G  O
                    1  1  1  1
                 (((H, C), G), O);

       usedata = 1    * 0: no data (prior); 1:seq like
         nloci = 53    * number of data sets in seqfile

     cleandata = 0    * remove sites with ambiguity data (1:yes, 0:no)?

    thetaprior =  2 2000   # gamma(a, b) for theta
      tauprior = 14 1000   # gamma(a, b) for root tau & Dirichlet(a) for other tau's

*     locusrate = 0 2.0   # (0: No variation, 1: estimate, 2: from file) & a_Dirichlet (if 1)
*      heredity = 0 4 4   # (0: No variation, 1: estimate, 2: from file) & a_gamma b_gamma (if 1)
* sequenceerror = 0 0 0 0 : 0.05 1   # sequencing errors: a_gamma, b_gamma

*       finetune = 0: 0.5 0.002 0.0006  0.0004 0.06 0.2 1.0  # auto (0 or 1): finetune for GBtj, GBspr, theta, tau, mix, locusrate, seqerr

       finetune = 1: .01 .01 .01 .01 .01 .01 .01 .01  # auto (0 or 1): finetune for GBtj, GBspr, theta, tau, mix, locusrate, seqerr

         print = 1 0 0 0   * MCMC samples, locusrate, heredityscalars Genetrees
        burnin = 2000
      sampfreq = 2
       nsample = 20000

*** Note: Make your window wider (140 columns) before running the program.
