  ======  ABOUT sqdifC =====
 Sqdif reads nucleotide sequences of aligned multifasta format
 and calculates the sequence difference of all possible
 sequence pairs based on the method of Miyata and Yasunaga (1980).

 The imput sequence should be nucleotide sequence, 
 and aligned considering Amino acid frame.

Reference:
 Miyata, T., and Yasunaga, T., 1980. Molecular evolution of mRNA: A method
 for estimating evolutionary rates of synonimous and amino acid substitutions
 from homologous nucleotide sequences and its application. J. Mol. Evol.
16: 23-36

 Sqdif is originally written in FORTRAN in 1982 and re-written
 in C in 2008 by Genome Information research Center (this version).

USAGE:
-i Infile. stdin if default.
-o Outfile. stdout if default. (not yet implemented)
-f Output format text (T) or tabular (X). T if default.
-w Window size to. This value must be dividable by 3. Zero if default.
-s Slide size. This value must be dividable by 3. Zero if default.
-r Consider substitution route in counting synonimous sites.
     [T/F] T if default. F if Nei-Gojobori method.
-c [T/F] Calculate difference of all combination [T] or Compare the first sequence with others [F].  T if default.
-W Substitution route is weighted. [T/F] T if default. F if Nei-Gojobori method.
-T Eliminate route of ternilal codon when creating the table of synonimous sites.
     [T/F]. F if default. T if Nei-Gojobori method.
-C Codon table file. Currently 'universal.codon'.
-S Substitution table file. Currently 'MiyataYasunagaTable_real.csv'.
-D show debug information [T/F] T if default.
