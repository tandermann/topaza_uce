          seed =  -1

       seqfile = phryno.txt
      Imapfile = phryno3s.Imap.txt
       outfile = out.txt
      mcmcfile = mcmc.txt

 speciesdelimitation = 0         * fixed species tree
*  speciesdelimitation = 1 0 2   * speciesdelimitation algorithm0(e)
  speciesdelimitation = 1 1 2 1  * speciesdelimitation algorithm1 finetune (a m)
          speciestree = 1        * species-tree by NNI

speciesmodelprior = 1         * 0: uniform labeled histories; 1:uniform rooted trees

  species&tree = 3 B Ce Co
                  82 34 20
                  (B, (Ce, Co));

       usedata = 1    * 0: no data (prior); 1:seq like
         nloci = 2    * number of data sets in seqfile

     cleandata = 0    * remove sites with ambiguity data (1:yes, 0:no)?

    thetaprior = 2 100   # gamma(a, b) for theta
      tauprior = 2 1000   # gamma(a, b) for root tau & Dirichlet(a) for other tau's

       finetune = 1: .05 .0001 .005 .0005 .2 .01 .01 .01  # auto (0 or 1): finetune for GBtj, GBspr, theta, tau, mix, locusrate, seqerr

         print = 1 0 0 1  * print MCMC samples, locusrate, heredityscalars, printGtree
        burnin = 2000
      sampfreq = 2
       nsample = 20000
