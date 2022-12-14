This software is used to construct A-Bruijn graph given pair-wise alignment of a reference genome. 

This was developed for a journal article published in 2004: 

Pevzner, Paul A., Haixu Tang, and Glenn Tesler. "De novo repeat classification and fragment assembly." Genome research 14.9 (2004): 1786-1796.

The original software was modified slightly to be used in conjunction with GRASPER ( https://github.com/COL-IU/GRASPER )

************************************
RepGraph is a package to produce 
a repeat graph from a genomic sequence.
************************************

**Compile/Installation
> make install



The user should have at least one of the
following programs to generate alignments
between the genomic sequence and itself:
blast, patternhunter,
or reputer. The following programs can
transform the alignment produced by the
above programs into the alignment file
my program required respectively:
bl2aln -- for blast result
Syntax: bl2aln -i InpFile -o outfile [-l min_leg -d min_id]
-i InpFile: The input file name of reads
-o OutFile: output alignment file
-l min_leg: minimum length of repeats (optional, default 500)
-d min_id: minimum identity of repeats (optional, default 0.95)

ph2aln -- for patternhunter
Syntax: ph2aln -i InpFile -s seqfile -o outfile [-l min_leg -d min_id]
-i InpFile: The input file name of alignments
-s SeqFile: The input file name of genomic sequence
-o OutFile: output alignment file
-l min_leg: minimum length of repeats (optional, default 500)
-d min_id: minimum identity of repeats (optional, default 0.95)

reputer2aln -- for reputer
Syntax: the same as ph2aln.

The output of the above result will be then used by 
program repeat_sin to generate repeat graph from the 
alignment file (see above).
Syntax: repeat_sin -i InpFile -s SeqFile -o outfile [-l shortleg] >tmp.log
-i InpFile: The input file name of alignments
-s SeqFile: The input file name of reads
-o Outfile: The output file name of contigs
-l shortleg (optional): default 50

From the standard output (tmp.log), the user can
see the statistics of the repeat graph and more
important, the embedding of the genomic
sequence on the repeat graph!
The repeat graph is slightly beautified by
remove the short edge (<shortleg).

There are two side programs in this directory
that are requested by Degui and Alkes.
simplifygraph -- simplify a provided repeat graph
	with the algorithm used in repeat_sin;
repeat_sin_read -- construct repeat graph
	based on the user-defined equivalent
	positions (instead of alignments
	produced by one of the above program)
