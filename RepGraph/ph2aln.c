/**
>HEADER
    Copyright (c) 2004 Haixu Tang hatang@indiana.edu

    This file is part of the RepGraph package.

    RepGraph is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    RepGraph is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with RepGraph.  If not, see <http://www.gnu.org/licenses/>.
<HEADER
**/

#include <stdinc.h>
#include <param.h>
#include <extfunc.h>

char	inpfile[100], outfile[100], qualfile[100], overlapfile[100], seqfile[100];
int	qualinp;
int	min_leg;
double	min_id;

void initenv(int argc, char **argv);

main(int argc, char **argv)
{
	int	i, j, k, l, m, n;
	char	**src_seq, **src_name;
	int	*len_seq, num_seq;
	char	temp[100];
	ALIGN	**align, *aln, *aln0;
	FILE	*fp;

	readpar();
	random1(&idum);
	initenv(argc, argv);

/*	Input the length of the reads (required) */

	len_seq = (int *) ckalloc(2 * sizeof(int));
	src_seq = (char **) ckalloc(2 * sizeof(char *));
	src_name = (char **) ckalloc(1 * sizeof(char *));
	src_name[0] = (char *) ckalloc(100 * sizeof(char));

	fp = ckopen(seqfile, "r");
	num_seq = readseq1by1(src_seq, src_name, len_seq, fp);
	fclose(fp);
	printf("Genome length: %d\n", len_seq[0]);

/*	Make reverse complements of input sequences rev(i) --> i + num_seq	*/

	len_seq[1] = len_seq[0];
	src_seq[1] = (char *) ckalloc(len_seq[0] * sizeof(char));
	for(j = 0; j < len_seq[0]; j ++)	{
		src_seq[1][j] = rev(src_seq[0][len_seq[0] - j - 1]);
	}

/*      read in pairwise alignments by Reputer	*/

	align = (ALIGN **) ckalloc(2 * sizeof(ALIGN *));
	fp = ckopen(inpfile, "r");
	n = readph(align, src_seq, len_seq, fp, min_leg, min_id);
	fclose(fp);
	printf("# alignments input: %d.\n", n);

/*	Write alignments	*/

	fp = ckopen(outfile, "w");
	for(m = 0; m < 2; m ++)	{
		n = size_align(align[m]);
		fwrite(&n, sizeof(int), 1, fp);
		aln = align[m];
		while(aln)	{
			fwrite(&(aln -> reads[1]), sizeof(int), 1, fp);
			fwrite(&(aln -> mis_match), sizeof(int), 1, fp);
			fwrite(&(aln -> length), sizeof(int), 1, fp);
			fwrite(aln -> pos[0], sizeof(int), aln -> length, fp);
			fwrite(aln -> pos[1], sizeof(int), aln -> length, fp);
			aln0 = aln -> next;
			free((void *) aln -> pos[0]);
			free((void *) aln -> pos[1]);
			free((void *) aln);
			aln = aln0;
		}
	}
	fclose(fp);
	printf("Done...\n");

	free((void **) align);
	for(i = 0; i < 2 * num_seq; i ++)	{
		free((void *) src_seq[i]);
	}
	for(i = 0; i < num_seq; i ++)	{
		free((void *) src_name[i]);
	}
	free((void **) src_seq);
	free((void **) src_name);
	free((void *) len_seq);
}

void initenv(int argc, char **argv)
{
	int copt;
	int inpseq, outseq;
	extern char *optarg;

	min_leg = 500;
	min_id = 0.99;
	inpseq = qualinp = outseq = 0;

	while ((copt=getopt(argc,argv,"i:o:l:d:s:")) != EOF)	{
		switch(copt) {
			case 'i':
			  inpseq = 1;
			  sscanf(optarg,"%s", inpfile);
			  continue;
			case 'o':
			  outseq = 1;
			  sscanf(optarg,"%s", outfile);
			  continue;
			case 's':
			  sscanf(optarg,"%s", &seqfile);
			  continue;
			case 'l':
			  sscanf(optarg,"%d", &min_leg);
			  continue;
			case 'd':
			  sscanf(optarg,"%lf", &min_id);
			  continue;
			default:
			  printf("ph2aln -i InpFile -s seqfile -o outfile [-l min_leg -d min_id]\n");
			  printf("-i InpFile: The input file name of alignments\n");
			  printf("-s SeqFile: The input file name of genomic sequence\n");
			  printf("-o OutFile: output alignment file\n");
			  printf("-l min_leg: minimum length of repeats.\n");
			  printf("-d min_id: minimum identity of repeats.\n");
			  exit(-1);
		}
		optind--;
	}

	if(inpseq == 0 || outseq == 0)	{
		printf("ph2aln -i InpFile -s seqfile -o outfile [-l min_leg -d min_id]\n");
		printf("-i InpFile: The input file name of alignments\n");
		printf("-s SeqFile: The input file name of genomic sequence\n");
		printf("-o OutFile: output alignment file\n");
		printf("-l min_leg: minimum length of repeats.\n");
		printf("-d min_id: minimum identity of repeats.\n");
		exit(-1);
	}
}
