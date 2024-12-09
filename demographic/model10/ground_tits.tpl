//Parameters for the coalescence simulation program : fsc27
5 samples to simulate :
//Population effective sizes (number of genes)
NSRE1
NSRE2
NSRE3
NSPLT
NSLE
//Samples sizes and samples age
20
20
20
20
20
//Growth rates : negative growth implies population expansion
0
0
0
0
0
//Number of migration matrices : 0 implies no migration between demes
2
//Migration matrix 0
0	0	0	0	0
0	0	0	0	0
0	0	0	0	0
0	0	0	0	M34
0	0	0	M43	0
//Migration matrix 2
0	0	0	0	0
0	0	0	0	0
0	0	0	0	0
0	0	0	0	0
0	0	0	0	0
//historical event: time, source, sink, migrants, new deme size, growth rate, migr mat index
6 historical event
TBOT1   2  2   0 RANCSRE3   0 0
TBOT1   1  1   0 RANCSRE2   0 0
TDIV1   4  3   1 RANC1      0 1
TDIV1   2  3   1 1      0 1
TDIV1   1  3   1 1      0 1
TDIV2   0  3   1 RANC2      0 1
//Number of independent loci [chromosome]
1 0
//Per chromosome: Number of contiguous linkage Block: a block is a set of contiguous loci
1
//per Block:data type, number of loci, per gen recomb and mut rates
FREQ 1 0 3.3e-9


