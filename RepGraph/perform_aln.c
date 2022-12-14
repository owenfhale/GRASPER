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

#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <perdef.h>
#include <param.h>
#include <extvab.h>
#include <extfunc.h>

#define MAX_NUM_ALN 40000

int min_length;
int tel_dist, centro_dist;

char inpfile[100], lenfile[100], alnfile[100], outfile[100];

void initenv(int argc, char **argv);
int segcompar(SEGMENT *a, SEGMENT *b);
int input_pair(FILE *fp, int *dir, int *coord1, int *coord2, char **chrname, int num_chro);
double cal_id(char **alnseq, int length);

main(int argc, char **argv)
{
	int	i, j, k, l, m, n, centro[200], len_chroseq[100];
	int	*coord1, *coord2, dir[2];
	char	**seq, **seqname;
	char	*seq1, *seq2;
	int	len_seq1, len_seq2;
	char	**alnseq;
	char	**chrname;
	int	num_pair;
	int	num_chro;
	int	length;
	SEGMENT	*segment;
	int	num_segment;
	int	*repeats, num_repeats;
	char	name[1000], temp[1000], c;
	double	id;
	FILE	*fp;

	readpar();
	initenv(argc, argv);

/*	Input chromsomal information	*/
	chrname = alloc_name(100, 100);
	fp = ckopen(lenfile, "r");
	num_chro = read_chro_centro(chrname, len_chroseq, centro, fp);
	fclose(fp);

/*	input alignment pairs	*/

	coord1 = (int *) ckalloc(3000 * sizeof(int));
	coord2 = (int *) ckalloc(3000 * sizeof(int));
	fp = ckopen(inpfile, "r");
	num_pair = input_pair(fp, dir, coord1, coord2, chrname, num_chro);
	fclose(fp);
	printf("Number of pair input: %d\n", num_pair);

	alnseq = (char **) ckalloc(2 * sizeof(char *));
	for(i = 0; i < num_pair; i ++)	{
		seq1 = (char *) ckalloc((abs(coord1[3 * i + 1] - coord1[3 * i + 2]) + 1) * sizeof(char));
		seq2 = (char *) ckalloc((abs(coord2[3 * i + 1] - coord2[3 * i + 2]) + 1) * sizeof(char));
		len_seq1 = getchseq(seq1, &coord1[3 * i], dir[0], chrname);
		len_seq2 = getchseq(seq2, &coord2[3 * i], dir[1], chrname);
printf("i %d len_seq1 %d len_seq2 %d\n", i, len_seq1, len_seq2);
		length = align2seq(seq1, len_seq1, seq2, len_seq2, alnseq);
printf("length %d\n", length);
		free((void *) seq1);
		free((void *) seq2);
		id = cal_id(alnseq, length);
		printf("%s %d %d %s %d %d: %f\n", chrname[coord1[3 * i]], coord1[3 * i + 1], coord1[3 * i + 2],
				chrname[coord2[3 * i]], coord2[3 * i + 1], coord2[3 * i + 2], id);
		sprintf(temp, "tmp/%s:%d-%d_%s:%d-%d.aln", chrname[coord1[3 * i]], coord1[3 * i + 1], coord1[3 * i + 2],
				chrname[coord2[3 * i]], coord2[3 * i + 1], coord2[3 * i + 2]);
		fp = ckopen(temp, "w");
		fprintf(fp, "%d %d\n", 2, length);
		sprintf(temp, "%s:%d-%d", chrname[coord1[3 * i]], coord1[3 * i + 1], coord1[3 * i + 2]);
		outputseq(fp, alnseq[0], length, temp);
		sprintf(temp, "%s:%d-%d", chrname[coord2[3 * i]], coord2[3 * i + 1], coord2[3 * i + 2]);
		outputseq(fp, alnseq[1], length, temp);
		fclose(fp);
		free((void *) alnseq[0]);
		free((void *) alnseq[1]);
	}

	free((void *) coord1);
	free((void *) coord2);
	free((void **) alnseq);
	chrname = free_name(chrname, 100);
}

int input_pair(FILE *fp, int *dir, int *coord1, int *coord2, char **chrname, int num_chro)
{
	int	i, j, k, l;
	char	str[500], tmp1[100], tmp2[100];

	k = 0;
	while(fgets(str, 490, fp))	{
		sscanf(str, "%s%d%d%s%d%d", tmp1, &coord1[3 * k + 1], &coord1[3 * k + 2],
			tmp2, &coord2[3 * k + 1], &coord2[3 * k + 2]);
		coord1[3 * k] = findgenname(tmp1, chrname, num_chro);
		coord2[3 * k] = findgenname(tmp2, chrname, num_chro);
		k ++;
	}
	return(k);
}

double cal_id(char **alnseq, int length)
{
	int	i, k;

	k = 0;
	for(i = 0; i < length; i ++)	{
		if(alnseq[0][i] == alnseq[1][i])	{
			k ++;
		}
	}
	return(((double) k) / length);
}

void initenv(int argc, char **argv)
{
	int copt;
	int inpseq, outseq;
	extern char *optarg;

	min_length = 2000;
	centro_dist = 5000000;
	inpseq = outseq = 0;

	while ((copt=getopt(argc,argv,"i:L:l:t:c:")) != EOF)	{
		switch(copt) {
			case 'i':
			  inpseq = 1;
			  sscanf(optarg,"%s", inpfile);
			  continue;
			case 'L':
			  sscanf(optarg,"%d", &min_length);
			  continue;
			case 'l':
			  sscanf(optarg,"%s", lenfile);
			  continue;
			default:
			  printf("perform_aln -i InpFile -l LenFile [-L min_leg ]\n");
			  printf("-i InpFile: The input file name of segments\n");
			  printf("-l lenfile: inpput file for chromosome length.\n");
			  printf("-L min_length: minimum interval length between repeats\n");
			  exit(-1);
		}
		optind--;
	}

	if(inpseq == 0)	{
		printf("perform_aln -i InpFile -l LenFile [-L min_leg ]\n");
		printf("-i InpFile: The input file name of segments\n");
		printf("-l lenfile: inpput file for chromosome length.\n");
		printf("-L min_length: minimum interval length between repeats\n");
		exit(-1);
	}
}
