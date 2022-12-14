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

#define MAX_NUM_ALN 40000
#define INTV 500

char	inpfile[100], outfile[100], qualfile[100], overlapfile[100], seqfile[100], lenfile[100], segfile[100],
	cshfile[100];
int	qualinp;
int	min_leg;
double	min_id;

void initenv(int argc, char **argv);
char chk_aln(int *posindex, int index1, int index2);
char chk_ovp(int *posindex1, int *posindex2);
int getclass(int *tmpclass, int size_tmpclass, char *label, int index, int **alnindex, int *num_aln);

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
	FILE	*fp, *fp1, *fp2;

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

/*	Classify alignment	*/

	posindex = (int *) ckalloc(num_seq * 6 * sizeof(int));
	for(m = 0; m < num_seq; m ++)	{
		transform_len(seqname[2 * m], chrname, num_chro, len_chro, &posindex[6 * m]);
		transform_len(seqname[2 * m + 1], chrname, num_chro, len_chro, &posindex[6 * m + 3]);
	}

	tmpaln = (int *) ckalloc(num_seq * sizeof(int));
	alnindex = (int **) ckalloc(num_seq * sizeof(int *));
	num_aln = (int *) ckalloc(num_seq * sizeof(int));
	n = 0;
	for(m = 0; m < num_seq; m ++)	{
		size_tmp = 0;
		for(i = 0; i < num_seq; i ++)	{
			if(i == m)	continue;
			c = chk_aln(posindex, m, i);
			if(c)	{
				tmpaln[size_tmp ++] = i;
			}
		}
		alnindex[m] = (int *) ckalloc(size_tmp * sizeof(int));
		num_aln[m] = size_tmp;
		for(i = 0; i < size_tmp; i ++)	{
			alnindex[m][i] = tmpaln[i];
		}
		n += size_tmp;
	}
	free((void *) tmpaln);
	printf("Matrix computed, total overlaps: %d\n", n);
	free((void *) posindex);

	alnclass = (int **) ckalloc(num_seq * sizeof(int *));
	tmpclass = (int *) ckalloc(num_seq * sizeof(int));
	size_class = (int *) ckalloc(num_seq * sizeof(int));
	label = (char *) ckalloc(num_seq * sizeof(char));
	num_class = 0;
	for(m = 0; m < num_seq; m ++)	{
		if(label[m] == 1)	continue;
		size_tmpclass = 0;
		size_tmpclass = getclass(tmpclass, size_tmpclass, label, m, alnindex, num_aln);
		alnclass[num_class] = (int *) ckalloc(size_tmpclass * sizeof(int));
		for(i = 0; i < size_tmpclass; i ++)	{
			alnclass[num_class][i] = tmpclass[i];
		}
		size_class[num_class] = size_tmpclass;
		num_class ++;
	}
	for(m = 0; m < num_seq; m ++)	{
		free((void *) alnindex[m]);
	}
	free((void **) alnindex);
	free((void *) num_aln);
	free((void *) label);
	free((void *) tmpclass);
	printf("%d class found\n", num_class);

/*	Write alignments	*/

	k = m = n = 0;
	for(i = 0; i < num_class; i ++)	{
		if(m < size_class[i])	{
			m = size_class[i];
		}
		if(size_class[i] == 1)	n ++;
		k += size_class[i];
	}
	printf("max %d --2 class %d, bigger class %d  total %d\n", m, n, num_class - n, k);
	fflush(stdout);

	fp2 = ckopen("watch", "w");
	for(i = 0; i < num_class; i ++)	{
		fprintf(fp2, "%d %d\n", i, size_class[i]);
	}
	fclose(fp2);

	sprintf(cshfile, "runall.csh");
	fp2 = ckopen(cshfile, "w");
	fprintf(fp2, "#! /usr/bin/csh\n\n");
	fp = ckopen(outfile, "w");
	for(i = 0; i < num_class; i ++)	{
		if(size_class[i] == 1)	{
			outputseq(fp, seq[2 * alnclass[i][0]], len_seq[alnclass[i][0]], seqname[2 * alnclass[i][0]]);
			outputseq(fp, seq[2 * alnclass[i][0] + 1], len_seq[alnclass[i][0]], seqname[2 * alnclass[i][0] + 1]);
		} else	{
			fprintf(fp2, "../run1by1.csh tmp-%d.aln\n", i + 1);
			sprintf(temp, "split/tmp-%d.aln", i + 1);
			fp1 = ckopen(temp, "w");
			for(j = 0; j < size_class[i]; j ++)	{
				outputseq(fp1, seq[2 * alnclass[i][j]], len_seq[alnclass[i][j]], seqname[2 * alnclass[i][j]]);
				outputseq(fp1, seq[2 * alnclass[i][j] + 1], len_seq[alnclass[i][j]], seqname[2 * alnclass[i][j] + 1]);
			}
			fclose(fp1);
		}
	}
	fclose(fp);
	fclose(fp2);
	printf("Done...\n");

	for(i = 0; i < num_class; i ++)	{
		free((void *) alnclass[i]);
	}
	free((void **) alnclass);
	free((void *) size_class);

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
	char	c;

	c = chk_ovp(&posindex[6 * index1], &posindex[6 * index2]);
	if(c == 1)	{
		return(c);
	}
	c = chk_ovp(&posindex[6 * index1], &posindex[6 * index2 + 3]);
	if(c == 1)	{
		return(c);
	}
	c = chk_ovp(&posindex[6 * index1 + 3], &posindex[6 * index2]);
	if(c == 1)	{
		return(c);
	}
	c = chk_ovp(&posindex[6 * index1 + 3], &posindex[6 * index2 + 3]);
	if(c == 1)	{
		return(c);
	}
	return(0);
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


int getclass(int *tmpclass, int size_tmpclass, char *label, int index, int **alnindex, int *num_aln)
{
	int	i, j, k, l, m, n;

	tmpclass[size_tmpclass ++] = index;
	label[index] = 1;
	for(i = 0; i < num_aln[index]; i ++)	{
		if(label[alnindex[index][i]] == 1)	{
			continue;
		}
		size_tmpclass = getclass(tmpclass, size_tmpclass, label, alnindex[index][i], alnindex, num_aln);
	}
	return(size_tmpclass);
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
			case 'l':
			  sscanf(optarg,"%d", &min_leg);
			  continue;
			case 'c':
			  sscanf(optarg,"%s", lenfile);
			  continue;
			case 'd':
			  sscanf(optarg,"%lf", &min_id);
			  continue;
			default:
			  printf("string2aln -i InpFile -c LenFile -o outfile [-l min_leg -d min_id]\n");
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
		printf("string2aln -i InpFile -c LenFile -o outfile [-l min_leg -d min_id]\n");
		printf("-i InpFile: The input file name of alignments\n");
		printf("-c lenfile: inpput file for chromosome length.\n");
		printf("-o OutFile: output alignment file\n");
		printf("-l min_leg: minimum length of repeats.\n");
		printf("-d min_id: minimum identity of repeats.\n");
		exit(-1);
	}
}
