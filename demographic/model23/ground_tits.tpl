//Parameters for the coalescence simulation program : fsc27
3 samples to simulate :
//Population effective sizes (number of genes)
NNRE
NNPLT
NNLE
//Samples sizes and samples age
18
18
18
//Growth rates : negative growth implies population expansion
0
0
0
//Number of migration matrices : 0 implies no migration between demes
3
//Migration matrix 0
0	0	0
0	0	M12
0	M21	0
//Migration matrix 1
0	M01	0
M10	0	M12
0	M21	0
//Migration matrix 2
0	0	0
0	0	0
0	0	0
//historical event: time, source, sink, migrants, new deme size, growth rate, migr mat index
3 historical event
TBOT1   0  0   0 RANCNRE    0 1
TDIV1   2  1   1 RANC1      0 2
TDIV1   0  1   1 1      0 2
//Number of independent loci [chromosome]
1 0
//Per chromosome: Number of contiguous linkage Block: a block is a set of contiguous loci
1
//per Block:data type, number of loci, per gen recomb and mut rates
FREQ 1 0 3.3e-9


