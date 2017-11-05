          seed =  -1

       seqfile = phryno.txt
      Imapfile = phryno5s.Imap.txt
       outfile = out.txt
      mcmcfile = mcmc.txt

* speciesdelimitation = 0 * fixed species tree
* speciesdelimitation = 1 0 2 * speciesdelimitation algorithm1 finetune (a m)
  speciesdelimitation = 1 1 2 0.5  * speciesdelimitation algorithm1 finetune (a m)
*          speciestree = 1 * species-tree by NNI

*speciesmodelprior = 1         * 0: uniform labeled histories; 1:uniform rooted trees

  species&tree = 5  NCA SCA NBC CBC SBC
                     18  44  20  34  20
                ((NCA, ((SCA, NBC), CBC)), SBC);

*                (((NCA, SCA), (NBC, CBC)), SBC);


       usedata = 1    * 0: no data (prior); 1:seq like
         nloci = 2    * number of data sets in seqfile

     cleandata = 0    * remove sites with ambiguity data (1:yes, 0:no)?

    thetaprior = 2 1000   # gamma(a, b) for theta
*    thetaprior = 2 100   # BY2013 used this unreasonable prior.  gamma(a, b) for theta
      tauprior = 2 1000   # gamma(a, b) for root tau & Dirichlet(a) for other tau's

      finetune = 1: .01 .0001 .005 .0005 .2 .01 .01 .01  # auto (0 or 1): finetune for GBtj, GBspr, theta, tau, mix, locusrate, seqerr

         print = 1 0 0 0   * MCMC samples, locusrate, heredityscalars Genetrees
        burnin = 4000
      sampfreq = 4
       nsample = 200000
