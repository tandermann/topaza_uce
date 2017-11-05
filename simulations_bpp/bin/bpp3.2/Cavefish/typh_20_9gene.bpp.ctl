          seed =  -1

       seqfile = typh_20_9gene.txt
      Imapfile = typh_20_9gene.Imap.txt
       outfile = out.txt
      mcmcfile = mcmc.txt

 speciesdelimitation = 0 * fixed species tree
* speciesdelimitation = 1 0 2    * speciesdelimitation algorithm0 and finetune(e)
 speciesdelimitation = 1 1 2 1 * speciesdelimitation algorithm1 finetune (a m)
         speciestree = 1 * species-tree by NNI

*   speciesmodelprior = 1  * 0: uniform LH; 1:uniform rooted trees; 2: uniformSLH; 3: uniformSRooted

  species&tree = 7  A  B  C  D  E  F  Sp
                    3  2  3  2  2  8   2
		    ((((B,A),E),((D,F),C)),Sp);


       usedata = 1    * 0: no data (prior); 1:seq like
         nloci = 1 5    * (the last gene is mt NADH2)

     cleandata = 0    * remove sites with ambiguity data (1:yes, 0:no)?

    thetaprior = 2 1000    # gamma(a, b) for theta
      tauprior = 2 1000 1  # gamma(a, b) for root tau & Dirichlet(a) for other tau's

*      heredity = 1 4 4
*     locusrate = 1 5

       finetune = 1: 2 0.0005 0.0020 0.00005 0.05 0.33  # finetune for GBtj, GBspr, theta, tau, mix, locusrate, seqerr

         print = 1 0 0 0   * MCMC samples, locusrate, heredityscalars Genetrees
        burnin = 400
      sampfreq = 2
       nsample = 2000
