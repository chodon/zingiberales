found one frameshift in GSMUA_Achr11T22900_001 for CS54_Musa_coccinea_Musaceae2bcallimusa 
 -- this sequence is probably in error

found by searching initial macse output for internal ! by:
grep -h '!' *NT.fasta | sed 's/^-*\![NACGT!]\{1,2\}\|\!\{1,2\}[ACGTN]\{1,2\}-*$//g' | grep '!' > internalQseqs.txt
grep -B 3 -f internalQseqs.txt 

found 21 seqs with internal Qs in 13 genes, one of the genes is eliminated with LBR already

and looking for frameshifts -- where the sequences do not line up and the amino acid is messed up too
