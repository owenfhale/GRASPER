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
#define INTV -1000

char	inpfile[100], targfile[100], lenfile[100];
int	qualinp;
int	min_leg;
double	min_id;

void initenv(int argc, char **argv);
char chk_ovp(int *posindex1, int *posindex2);
void chk_list(char *label, int *posindex, int num_seq, int *intvlist);

main(int argc, char **argv)
{
	int	i, j, k, l, m, n;
	char	c;
	int	len_chro[100];
	int	num_chro, num_intv;
	char	temp[100];
	char	**chrname;
	char	**regname, c1, c2;
	char	*label;
	char	**seq, **seqname;
	int	*len_seq;
	int	num_seq;
	int	**intvlist;
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

/*	Input target intevals	*/

	fp = ckopen(targfile, "r");
	intvlist = (int **) ckalloc(100 * sizeof(int *));
	regname = (char **) ckalloc(100 * sizeof(char *));
	for(i = 0; i < 100; i ++)	{
		regname[i] = (char *) ckalloc(20 * sizeof(char));
		intvlist[i] = (int *) ckalloc(3 * sizeof(int));
	}
	num_intv = readtargintv(fp, intvlist, regname, len_chro, chrname, num_chro);
	fclose(fp);
	printf("number intevals: %d\n", num_intv);
	for(i = 0; i < num_intv; i ++)	{
		printf("%d %d %d\n", intvlist[i][0], intvlist[i][1], intvlist[i][2]);
	}

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
	printf("# alignments input: %d.\n", num_seq, len_chro, chrname, num_chro);

/*	Classify alignment	*/

	posindex = (int *) ckalloc(num_seq * 6 * sizeof(int));
	for(m = 0; m < num_seq; m ++)	{
		transform_len(seqname[2 * m], chrname, num_chro, len_chro, &posindex[6 * m]);
		transform_len(seqname[2 * m + 1], chrname, num_chro, len_chro, &posindex[6 * m + 3]);
	}

/*	Write alignments	*/

	label = (char *) ckalloc(num_seq * sizeof(char));
	for(i = 0; i < num_intv; i ++)	{
		for(j = 0; j < num_seq; j ++)	{
			label[j] = 0;
		}
		chk_list(label, posindex, num_seq, intvlist[i]);
		sprintf(temp, "%s.align", regname[i]);
		fp = ckopen(temp, "w");
		for(j = 0; j < num_seq; j ++)	{
			if(label[j] == 1)	{
				outputseq(fp, seq[2 * j], len_seq[j], seqname[2 * j]);
				outputseq(fp, seq[2 * j + 1], len_seq[j], seqname[2 * j + 1]);
			}
		}
		fclose(fp);
	}
	printf("Done...\n");

	free((void *) posindex);
	free((void *) label);

	for(i = 0; i < 100; i ++)	{
		free((void *) intvlist[i]);
		free((void *) regname[i]);
	}
	free((void **) intvlist);
	free((void **) regname);

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

void chk_list(char *label, int *posindex, int num_seq, int *intvlist)
{
	int	i, j, k;
	char	c1, c2;

	k = 0;
	for(j = 0; j < num_seq; j ++)	{
		if(label[j] == 1)	continue;
		c1 = chk_ovp(&posindex[6 * j], intvlist);
		c2 = chk_ovp(&posindex[6 * j + 3], intvlist);
		if(c1 || c2)	{
			label[j] = 1;
/*
			chk_list(label, posindex, num_seq, &posindex[6 * j]);
			chk_list(label, posindex, num_seq, &posindex[6 * j + 3]);
*/
			k ++;
		}
	}
	printf("%d alignments added\n", k);
}

char chk_ovp(int *posindex1, int *posindex2)
{
	int	i, j, k, l;
	int	range1[2], range2[2], nchro1, nchro2;

	nchro1 = posindex1[0];
	nchro2 = posindex2[0];
	if(nchro1 != nchro2)	{
		return(0);
	}
	range1[0] = posindex1[1];
	range1[1] = posindex1[2];
	range2[0] = posindex2[1];
	range2[1] = posindex2[2];
	if(range2[0] >= range1[0])	{
		if(range2[0] <= range1[1] - INTV)	{
			return(1);
		} else	{
			return(0);
		}
	} else 	{
		if(range1[0] <= range2[1] - INTV)	{
			return(1);
		} else	{
			return(0);
		}
	}
}

void initenv(int argc, char **argv)
{
	int copt;
	int inpseq, outseq;
	extern char *optarg;

	inpseq = qualinp = outseq = 0;

	while ((copt=getopt(argc,argv,"i:t:c:")) != EOF)	{
		switch(copt) {
			case 'i':
			  inpseq = 1;
			  sscanf(optarg,"%s", inpfile);
			  continue;
			case 't':
			  outseq = 1;
			  sscanf(optarg,"%s", targfile);
			  continue;
			case 'c':
			  sscanf(optarg,"%s", lenfile);
			  continue;
			default:
			  printf("grabaln -i InpFile -c LenFile -t targfile\n");
			  printf("-i InpFile: The input file name of alignments\n");
			  printf("-c lenfile: input file for chromosome length.\n");
			  printf("-t TargFile: input file of interesting regions\n");
			  exit(-1);
		}
		optind--;
	}

	if(inpseq == 0 || outseq == 0)	{
		printf("grabaln -i InpFile -c LenFile -t targfile\n");
		printf("-i InpFile: The input file name of alignments\n");
		printf("-c lenfile: input file for chromosome length.\n");
		printf("-t TargFile: input file of interesting regions\n");
		exit(-1);
	}
}
