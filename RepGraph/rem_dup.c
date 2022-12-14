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

#define MAX_NUM_ALN 100000

char	inpfile[100], outfile[100], qualfile[100], overlapfile[100], seqfile[100], lenfile[100], segfile[100];
int	qualinp;
int	min_leg;
double	min_id;

void initenv(int argc, char **argv);
char chk_aln(int *posindex, int index1, int index2);
char chk_ovp(int *pos1, int *pos2);

main(int argc, char **argv)
{
	int	i, j, k, l, m, n;
	char	c;
	int	len_chro[100];
	int	num_chro;
	char	temp[100];
	char	**chrname;
	char	*label;
	char	**seq, **seqname;
	int	*len_seq;
	int	num_seq;
	int	*tmpclass, **alnclass, *size_class, size_tmpclass, num_class;
	int	*tmpaln, *num_aln, **alnindex, size_tmp;
	int	*posindex;
	FILE	*fp, *fp1;

	readpar();
	random1(&idum);
	initenv(argc, argv);

/*	read in chromosome length	*/
	chrname = alloc_name(100, 100);
	fp = ckopen(lenfile, "r");
	num_chro = readlen(fp, len_chro, chrname);
	fclose(fp);
	printf("# chromosome %d\n", num_chro);

/*      read in pairwise alignments 	*/

	len_seq = (int *) ckalloc(MAX_NUM_ALN * sizeof(int));
	seq = (char **) ckalloc(2 * MAX_NUM_ALN * sizeof(char *));
	seqname = (char **) ckalloc(2 * MAX_NUM_ALN * sizeof(char *));
	for(i = 0; i < 2 * MAX_NUM_ALN; i ++)	{
		seqname[i] = (char *) ckalloc(100 * sizeof(char));
	}
	fp = ckopen(inpfile, "r");
	num_seq = readintv(seq, fp, len_seq, seqname);
	fclose(fp);
	printf("# alignments input: %d.\n", num_seq);

/*	Remove duplicated alignment	*/

	label = (char *) ckalloc(num_seq * sizeof(char));
	posindex = (int *) ckalloc(num_seq * 6 * sizeof(int));
	for(m = 0; m < num_seq; m ++)	{
		transform_len(seqname[2 * m], chrname, num_chro, len_chro, &posindex[6 * m]);
		transform_len(seqname[2 * m + 1], chrname, num_chro, len_chro, &posindex[6 * m + 3]);
	}
	for(m = 0; m < num_seq; m ++)	{
		if(label[m] == 1)	{
			continue;
		}
		for(i = m + 1; i < num_seq; i ++)	{
			if(label[i] == 1)	{
				continue;
			}
			c = chk_aln(posindex, m, i);
			if(c)	{
				label[i] = 1;
				break;
			}
		}
	}
	free((void *) posindex);
	printf("Matrix computed\n");

/*	Write alignments	*/

	fp = ckopen(outfile, "w");
	for(i = 0; i < num_seq; i ++)	{
		if(label[i] == 0)	{
			outputseq(fp, seq[2 * i], len_seq[i], seqname[2 * i]);
			outputseq(fp, seq[2 * i + 1], len_seq[i], seqname[2 * i + 1]);
		}
	}
	fclose(fp);
	printf("Done...\n");

	free((void *) label);
	for(i = 0; i < 2 * num_seq; i ++)	{
		free((void **) seq[i]);
	}
	for(i = 0; i < 2 * MAX_NUM_ALN; i ++)	{
		free((void *) seqname[i]);
	}
	free((void **) seq);
	free((void **) seqname);
	free((void *) len_seq);
	chrname = free_name(chrname, 100);
}

char chk_aln(int *posindex, int index1, int index2)
{
	char	c1, c2;

	c1 = chk_ovp(&posindex[6 * index1], &posindex[6 * index2 + 3]);
	c2 = chk_ovp(&posindex[6 * index1 + 3], &posindex[6 * index2]);
	if(c1 == 1 && c2 == 1)	{
		return(1);
	}
	c1 = chk_ovp(&posindex[6 * index1], &posindex[6 * index2]);
	c2 = chk_ovp(&posindex[6 * index1 + 3], &posindex[6 * index2 + 3]);
	if(c1 == 1 && c2 == 1)	{
		return(1);
	}
	return(0);
}

char chk_ovp(int *pos1, int *pos2)
{
	int	i;

	for(i = 0; i < 3; i ++)	{
		if(pos1[i] != pos2[i])	break;
	}
	if(i == 3)	return(1);
	else		return(0);
}

void initenv(int argc, char **argv)
{
	int copt;
	int inpseq, outseq;
	extern char *optarg;

	min_leg = 500;
	min_id = 0.99;
	inpseq = qualinp = outseq = 0;

	while ((copt=getopt(argc,argv,"i:o:l:d:c:")) != EOF)	{
		switch(copt) {
			case 'i':
			  inpseq = 1;
			  sscanf(optarg,"%s", inpfile);
			  continue;
			case 'o':
			  outseq = 1;
			  sscanf(optarg,"%s", outfile);
			  continue;
			case 'c':
			  sscanf(optarg,"%s", lenfile);
			  continue;
			case 'l':
			  sscanf(optarg,"%d", &min_leg);
			  continue;
			case 'd':
			  sscanf(optarg,"%lf", &min_id);
			  continue;
			default:
			  printf("rem_dup -i InpFile -c LenFile -o outfile [-l min_leg -d min_id]\n");
			  printf("-i InpFile: The input file name of alignments\n");
			  printf("-c lenfile: inpput file for chromosome length.\n");
			  printf("-o OutFile: output alignment file\n");
			  printf("-l min_leg: minimum length of repeats.\n");
			  printf("-d min_id: minimum identity of repeats.\n");
			  exit(-1);
		}
		optind--;
	}

	if(inpseq == 0 || outseq == 0)	{
		printf("rem_dup -i InpFile -c LenFile -o outfile [-l min_leg -d min_id]\n");
		printf("-i InpFile: The input file name of alignments\n");
		printf("-c lenfile: inpput file for chromosome length.\n");
		printf("-o OutFile: output alignment file\n");
		printf("-l min_leg: minimum length of repeats.\n");
		printf("-d min_id: minimum identity of repeats.\n");
		exit(-1);
	}
}
