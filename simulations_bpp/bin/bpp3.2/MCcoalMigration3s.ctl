          seed = 12345

       seqfile = MySeq.txt   * comment out this line if you don't want seqs
      treefile = MyTree.tre   * comment out this line if you don't want trees
      Imapfile = MyImap.txt

  species&tree = 3  A  B  C 
                    1  1  1 
  ((A #1.1, B #1.1) : 0.5 #.12, C #0.3) : 1.0 #.3;

     migration = 5   * number of pops (order fixed by program)

             A     B     C    ABC   AB
      A      0     1.1   0    0     0
      B      1.1   0     0    0     0
      C      0     0     0    0     0
      ABC    0     0     0    0     0
      AB     0     0     0    0     0

   loci&length = 100 1000 * number of loci & number of sites at each locus
*     locusrate = 0.5     * alpha for gamma for locus rate
