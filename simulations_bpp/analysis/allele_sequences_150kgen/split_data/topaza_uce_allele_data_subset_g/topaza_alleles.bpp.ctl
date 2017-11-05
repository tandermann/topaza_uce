          seed =  -1

       seqfile = topaza_uce_allele_data_subset_g
      Imapfile = topaza_alleles.Imap.txt
       outfile = out.txt
      mcmcfile = mcmc.txt

   speciesdelimitation = 0
   speciestree = 0

  species&tree = 6  F  D  E  X  Y  Z
                    2  4  4  4  2  4
                 (((D, E),((Y, Z), X)), F);

       usedata = 1
         nloci = 82

     cleandata = 0
*no cleaning of the input alignments (a ‘1’ here would mean that ambiguities and gaps are eliminated form the input alignments

    thetaprior =  2 2000
* specifies the gamma prior G(a, b) for the theta parameters, with the mean to be a/b. In this case, the mean is 2/2000 = 0.001.

      tauprior = 14 1000
* divergence time parameter for the root in the species tree

* no heredity specifies (default=0) since theta is expected to be equal between all UCE loci
* no locusrate specified, because mutation rate is similar between the different UCE loci 

       finetune = 1: .01 .01 .01 .01 .01 .01 .01 .01
* the “1” before the colon means that the MCMC will do automatic adjustments by the program (other option ‘0’ means no adjustments). The following numbers are the initial step sizes (not exactly sure what difference this makes…)

         print = 1
        burnin = 50000
      sampfreq = 10
       nsample = 150000

*** Note: Make your window wider (140 columns) before running the program.
