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

char	inpfile[100], outfile[100], qualfile[100], overlapfile[100], seqfile[100], lenfile[100], segfile[100],
	cshfile[100], graphfile[100];
int	qualinp;
int	min_leg;
double	min_id;

void initenv(int argc, char **argv);
int regcompar(REGION *a, REGION *b);
int getclass_step(int *tmpclass, int size_tmpclass, int index, int *dupreg, REGION *region);
int getclass(int *tmpclass, int size_tmpclass, int *num_reg, int *label, int index, int **alnindex, int *num_aln, int *dupreg, REGION *region, int num_class);

main(int argc, char **argv)
{
	int	i, j, k, l, m, n;
	int	k1, k2, n1, n2, *last;
	char	c;
	int	len_chro[100];
	int	num_chro;
	int	maxi, maxk;
	char	temp[100], tmp1[500], tmp2[500];
	char	**chrname;
	int	*label;
	char	**seq, **seqname;
	int	*len_seq, *num_reg;
	int	num_seq;
	int	*tmpclass, **alnclass, *size_class, size_tmpclass, num_class, *dist;
	int	**tmpaln, *sizealn, **outaln, size_tmp;
	int	*dupreg, num_dupreg;
	REGION	*region;
	int	posindex[3];
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

	region = (REGION *) ckalloc(2 * num_seq * sizeof(REGION));
	n1 = n2 = 0;
	for(m = 0; m < num_seq; m ++)	{
		transform_len(seqname[2 * m], chrname, num_chro, len_chro, posindex);
		region[2 * m].chro = posindex[0] - 1;
		region[2 * m].range[0] = posindex[1];
		region[2 * m].range[1] = posindex[2];
		region[2 * m].index = 2 * m;
		transform_len(seqname[2 * m + 1], chrname, num_chro, len_chro, posindex);
		region[2 * m + 1].chro = posindex[0] - 1;
		region[2 * m + 1].range[0] = posindex[1];
		region[2 * m + 1].range[1] = posindex[2];
		region[2 * m + 1].index = 2 * m + 1;
		if(region[2 * m].chro == region[2 * m + 1].chro)	{
			n1 ++;
		} else	{
			n2 ++;
		}
	}
	printf("intrachromo: %d; interchromo: %d\n", n1, n2);
	qsort(region, 2 * num_seq, sizeof(REGION), (void *) regcompar);

/*	list duplicated regions	*/

	dupreg = (int *) ckalloc(2 * num_seq * sizeof(int));
	dist = (int *) ckalloc(num_chro * sizeof(int));
	num_dupreg = n = 0;
	m = 0;
	last = (int *) ckalloc(2 * num_seq * sizeof(int));
	for(i = 0; i < 2 * num_seq - 1; i ++)	{
		if(last[num_dupreg] < region[i].range[1])	last[num_dupreg] = region[i].range[1];
		if(region[i + 1].chro != region[i].chro || region[i + 1].range[0] - last[num_dupreg] > min_leg)	{
			k = last[num_dupreg] - region[m].range[0] + 1;
			if(k > n)	{
				n = k;
				posindex[0] = region[i].chro;
				posindex[1] = region[m].range[0];
				posindex[2] = last[num_dupreg];
			}
			dupreg[num_dupreg] = i + 1;
			dist[region[i].chro] ++;
			m = i + 1;
			num_dupreg ++;
		}
	}
	k = last[num_dupreg] - region[m].range[0] + 1;
	if(k > n)	{
		n = k;
		posindex[0] = region[i].chro;
		posindex[1] = region[m].range[0];
		posindex[2] = last[num_dupreg];
	}
	dupreg[num_dupreg ++] = 2 * num_seq - 1;
	printf("Total duplicated region: %d; Longest: %d (%s,%d-%d)\n", num_dupreg, n, chrname[posindex[0]], posindex[1], posindex[2]);
	for(i = 0; i < num_chro; i ++)	{
		printf("%s %d\n", chrname[i], dist[i]);
	}
	free((void *) last);

/*	Build the overlap graph	*/

	tmpaln = (int **) ckalloc(num_dupreg * sizeof(int *));
	outaln = (int **) ckalloc(num_dupreg * sizeof(int *));
	sizealn = (int *) ckalloc(num_dupreg * sizeof(int));
	k1 = 0;
	for(m = 0; m < num_dupreg; m ++)	{
		tmpaln[m] = (int *) ckalloc(num_dupreg * sizeof(int));
		outaln[m] = (int *) ckalloc(num_dupreg * sizeof(int));
		for(k = k1; k < dupreg[m]; k ++)	{
			n1 = region[k].index / 2;
			k2 = dupreg[m];
			for(i = m + 1; i < num_dupreg; i ++)	{
				for(j = k2; j < dupreg[i]; j ++)	{
					n2 = region[j].index / 2;
					if(n1 == n2)	{
						tmpaln[m][i] ++;
					}
				}
				k2 = dupreg[i];
			}
		}
		k1 = dupreg[m];
	}

	n = 0;
	for(m = 0; m < num_dupreg; m ++)	{
		for(i = m + 1; i < num_dupreg; i ++)	{
			tmpaln[i][m] = tmpaln[m][i];
			if(tmpaln[m][i] > 0)	{
				outaln[m][sizealn[m] ++] = i;
				outaln[i][sizealn[i] ++] = m;
				n ++;
			}
		}
	}
	printf("Matrix computed, total overlaps: %d\n", n);

	k = 0;
	for(m = 0; m < num_dupreg; m ++)	{
		if(sizealn[m] > k)	{
			k = sizealn[m];
			i = m;
		}
	}
	printf("Largest degree: %d, %s (%d) %d %d\n", k, chrname[region[dupreg[i]].chro], dupreg[i] - dupreg[i - 1] + 1,
			region[dupreg[i - 1]].range[0], region[dupreg[i]].range[1]);

/*	Derive overlap class	*/

	label = (int *) ckalloc(num_dupreg * sizeof(int));
	tmpclass = (int *) ckalloc(2 * num_seq * sizeof(int));
	alnclass = (int **) ckalloc(num_dupreg * sizeof(int *));
	size_class = (int *) ckalloc(2 * num_seq * sizeof(int));
	num_reg = (int *) ckalloc(num_dupreg * sizeof(int));
	num_class = 0;
	for(m = 0; m < num_dupreg; m ++)	{
		if(label[m] > 0)	continue;
		size_tmpclass = 0;
		size_tmpclass = getclass(tmpclass, size_tmpclass, num_reg, label, m, outaln, sizealn, dupreg, region, num_class);
		alnclass[num_class] = (int *) ckalloc(size_tmpclass * sizeof(int));
		for(i = 0; i < size_tmpclass; i ++)	{
			alnclass[num_class][i] = tmpclass[i];
		}
		size_class[num_class] = size_tmpclass;
		num_class ++;
	}

/*	Ouput overlap graph of the largest component	*/

	maxk = 0;
	for(i = 0; i < num_class; i ++)	{
		if(num_reg[i] > maxk)	{
			maxk = num_reg[i];
			maxi = i;
		}
	}

	l = 0;
	for(m = 0; m < num_dupreg - 1; m ++)	{
		posindex[1] = region[dupreg[m]].range[0];
		posindex[2] = region[dupreg[m + 1]].range[1];
		l += (posindex[2] - posindex[1] + 1);
	}
	printf("total length %d\n", l);

	n = l = 0;
	fp = ckopen(graphfile, "w");
	fprintf(fp, "digraph G {\n");
	fprintf(fp, "\tsize=\"8,8\";\n");
	for(m = 0; m < num_dupreg; m ++)	{
		if(label[m] != maxi + 1)	{
			continue;
		}
		posindex[0] = region[dupreg[m]].chro;
		posindex[1] = region[dupreg[m]].range[0];
		if(m < num_dupreg)	{
			posindex[2] = region[dupreg[m + 1]].range[1];
			l += (posindex[2] - posindex[1] + 1);
		}
		sprintf(tmp1, "%s,%dK", chrname[posindex[0]], posindex[1] / 1000);
		for(i = m + 1; i < num_dupreg; i ++)	{
			posindex[0] = region[dupreg[i]].chro;
			posindex[1] = region[dupreg[i]].range[0];
			sprintf(tmp2, "%s,%dK", chrname[posindex[0]], posindex[1] / 1000);
			tmpaln[i][m] = tmpaln[m][i];
			if(tmpaln[m][i] > 0)	{
				fprintf(fp, "\t%s -> %s [label=\"%d\"];\n", tmp1, tmp2, tmpaln[m][i]);
				n ++;
			}
		}
	}
	fprintf(fp, "}\n");
	fclose(fp);
	printf("Total length in the largest component: %d\n", l);
	printf("%d Largest component: %d alignments %d regions %d connections.\n", maxi, size_class[maxi], num_reg[maxi], n);
	fp2 = ckopen("watch", "w");
	for(i = 0; i < num_class; i ++)	{
		fprintf(fp2, "%d %d %d\n", i, size_class[i], num_reg[i]);
	}
	fclose(fp2);
printf("done\n");

	for(m = 0; m < num_dupreg; m ++)	{
		free((void *) tmpaln[m]);
	}
	free((void **) tmpaln);
	printf("Matrix computed, total overlaps: %d\n", n);

	free((void *) num_reg);
	for(m = 0; m < num_dupreg; m ++)	{
		free((void *) outaln[m]);
	}
	free((void **) outaln);
	free((void *) sizealn);

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

	free((void *) region);
	free((void *) dupreg);
	free((void *) dist);

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

int getclass(int *tmpclass, int size_tmpclass, int *num_reg, int *label, int index, int **alnindex, int *num_aln, int *dupreg, REGION *region, int num_class)
{
	int	i, j, k, l, m, n;

	num_reg[num_class] ++;
	size_tmpclass = getclass_step(tmpclass, size_tmpclass, index, dupreg, region);
	label[index] = num_class + 1;
	for(i = 0; i < num_aln[index]; i ++)	{
		if(label[alnindex[index][i]] > 0)	{
			continue;
		}
		size_tmpclass = getclass(tmpclass, size_tmpclass, num_reg, label, alnindex[index][i], alnindex, num_aln, dupreg, region, num_class);
	}
	return(size_tmpclass);
}

int getclass_step(int *tmpclass, int size_tmpclass, int index, int *dupreg, REGION *region)
{
	int	i, j, k, l;

	if(index == 0)	k = 0;
	else		k = dupreg[index - 1];
	for(i = k; i < dupreg[index]; i ++)	{
		if(region[i].index % 2 == 0)	{
			j = region[i].index / 2;
			tmpclass[size_tmpclass ++] = j;
		}
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

	while ((copt=getopt(argc,argv,"i:o:l:g:c:v:")) != EOF)	{
		switch(copt) {
			case 'i':
			  inpseq = 1;
			  sscanf(optarg,"%s", inpfile);
			  continue;
			case 'o':
			  outseq = 1;
			  sscanf(optarg,"%s", outfile);
			  continue;
			case 'g':
			  sscanf(optarg,"%s", graphfile);
			  continue;
			case 'l':
			  sscanf(optarg,"%d", &min_leg);
			  continue;
			case 'c':
			  sscanf(optarg,"%s", lenfile);
			  continue;
			default:
			  printf("classifyaln_seg -i InpFile -c LenFile -o outfile -g graphfile [-l min_leg]\n");
			  printf("-i InpFile: The input file name of alignments\n");
			  printf("-c lenfile: inpput file for chromosome length.\n");
			  printf("-o OutFile: output alignment file\n");
			  printf("-g GraphFile: output file of overlapgraph\n");
			  printf("-l min_leg: minimum length of allowed intervals.\n");
			  exit(-1);
		}
		optind--;
	}

	if(inpseq == 0 || outseq == 0)	{
		printf("classifyaln_seg -i InpFile -c LenFile -o outfile -g graphfile [-l min_leg]\n");
		printf("-i InpFile: The input file name of alignments\n");
		printf("-c lenfile: inpput file for chromosome length.\n");
		printf("-o OutFile: output alignment file\n");
		printf("-l min_leg: minimum length of allowed intervals.\n");
		printf("-g GraphFile: output file of overlapgraph\n");
		exit(-1);
	}
}

int regcompar(REGION *a, REGION *b)
{
	if(a -> chro > b -> chro || a -> chro == b -> chro &&
	   a -> range[0] > b -> range[0])	return(1);
	else				return(-1);
}
