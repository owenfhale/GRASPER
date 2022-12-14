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

#define par 30
#define MAX_NUM_ALN 40000
#define MIN_SEG_LEG 600
#define MAX_SUBREP 10000
#define MIN_SEG 2

int min_length;
int tel_dist, centro_dist;
int WIDTH;

char inpfile[100], lenfile[100], outfile[100];

void initenv(int argc, char **argv);
int segcompar(SEGMENT *a, SEGMENT *b);

main(int argc, char **argv)
{
	int	i, j, k, l, m, n, centro[200], len_chroseq[100];
	int	**coord;
	int	max_sub;
	int	num_seq;
	int	p1, p2;
	char	**subrepseq;
	int	*len_subrepseq;
	char	*treeseq;
	int	*stree;
	int	**intvlist, num_intv;
	int	MaxItem, nim;
	double  MaxEvoAve, *AveEvoDist, DiffEvo, EvoAve;
	double	*dist, *EvoDist;
	char	**chrname, **repnames;
	int	num_chro;
	int	**subrep, **rep, *index_subrep, *ancestor, num_subrep;
	char	**repseq;
	int	*len_repseq;
	SEGMENT	*segment;
	char	*singleton_mark;
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

/*	input segments	*/

	segment = (SEGMENT *) ckalloc(20000 * sizeof(SEGMENT));
	fp = ckopen(inpfile, "r");
	num_segment = input_segment(segment, chrname, num_chro, fp);
	fclose(fp);

/*	sort the segments	*/
	qsort((void *) segment, num_segment, sizeof(SEGMENT), (void *) segcompar);

/*	Define repeats from sub-repeats	*/

	repeats = (int *) ckalloc(num_segment * sizeof(int));
	num_repeats = 0;
	for(i = 0; i < num_segment; i ++)	{
		if(i == num_segment - 1 || segment[i + 1].pos[0] > segment[i].pos[1] + min_length ||
		   segment[i + 1].chro != segment[i].chro)	{
			repeats[num_repeats ++] = i;
		}
	}
	printf("%d repeats found.\n", num_repeats);
	printf("num_segment %d\n", num_segment);

/*	Mark non-paracentromic singletons	*/

	singleton_mark = (char *) ckalloc(num_segment * sizeof(char));
	k = 0;
	for(i = 0; i < num_repeats; i ++)	{
		p1 = segment[k].pos[0];
		p2 = segment[repeats[i]].pos[1];
		n = repeats[i] - k + 1;
		if((p1 > centro[2 * segment[k].chro + 1] + centro_dist &&
		    p2 < len_chroseq[segment[k].chro] - tel_dist ||
		    p2 < centro[2 * segment[k].chro] - centro_dist &&
		    p1 > tel_dist) && n <= MIN_SEG)	{
			for(j = k; j <= repeats[i]; j ++)	{
				singleton_mark[j] = 1;
			}
		}
		k = repeats[i];
	}

/*	get the sequences of each repeat	*/

	len_repseq = (int *) ckalloc(num_repeats * sizeof(int));
	repseq = (char **) ckalloc(num_repeats * sizeof(char *));
	coord = (int **) ckalloc(num_repeats * sizeof(int *));
	subrep = (int **) ckalloc(MAX_SUBREP * sizeof(int *));
	rep = (int **) ckalloc(MAX_SUBREP * sizeof(int *));
	for(i = 0; i < num_repeats; i ++)	{
		coord[i] = (int *) ckalloc(3 * sizeof(int));
	}
	for(i = 0; i < MAX_SUBREP; i ++)	{
		subrep[i] = (int *) ckalloc(num_segment * sizeof(int));
		rep[i] = (int *) ckalloc(num_segment * sizeof(int));
	}
	index_subrep = (int *) ckalloc(MAX_SUBREP * sizeof(int));
	k = -1;
	max_sub = 0;
	num_subrep = 0;
	for(i = 0; i < num_repeats; i ++)	{
		m = k + 1;
		coord[i][0] = segment[m].chro;
		k = repeats[i];
		for(j = m; j <= k; j ++)	{
			subrep[segment[j].eq_pos[0]][index_subrep[segment[j].eq_pos[0]]] = j;
			rep[segment[j].eq_pos[0]][index_subrep[segment[j].eq_pos[0]]] = i;
			index_subrep[segment[j].eq_pos[0]] ++;
			if(segment[j].eq_pos[0] > num_subrep)	{
				num_subrep = segment[j].eq_pos[0];
			}
		}
		coord[i][1] = segment[m].pos[0] + 1;
		coord[i][2] = segment[k].pos[1] + 2;
		repseq[i] = (char *) ckalloc((coord[i][2] - coord[i][1] + 1) * sizeof(char));
		len_repseq[i] = getchseq(repseq[i], coord[i], 0, chrname);
		if(max_sub < len_repseq[i])	{
			max_sub = len_repseq[i];
		}
	}
	num_subrep ++;
	printf("# subrepeats: %d\n", num_subrep);

/*	Evolutionary analysis of each sub-repeat	*/

	ancestor = (int *) ckalloc(num_subrep * sizeof(int));
	for(i = 0; i < num_subrep; i ++)	{
		ancestor[i] = -1;
	}
	fp = ckopen(outfile, "w");
	for(i = 0; i < num_subrep; i ++)	{
		if(index_subrep[i] <= 1)	continue;
		p1 = segment[subrep[i][0]].pos[0];
		p2 = segment[subrep[i][0]].pos[1];
		if(p2 - p1 < MIN_SEG_LEG)	{
			continue;
		}
		repnames = (char **) ckalloc(index_subrep[i] * sizeof(char *));
		len_subrepseq = (int *) ckalloc(index_subrep[i] * sizeof(int));
		subrepseq = (char **) ckalloc(index_subrep[i] * sizeof(char *));
		for(j = 0; j < index_subrep[i]; j ++)	{
			repnames[j] = (char *) ckalloc(100 * sizeof(char));
			subrepseq[j] = (char *) ckalloc(max_sub * sizeof(char));
		}
/*	Build subrepeat sequences	*/
		for(j = 0; j < index_subrep[i]; j ++)	{
			sprintf(repnames[j], "%s_%d-%d_R%d", chrname[segment[subrep[i][j]].chro],
				segment[subrep[i][j]].pos[0], segment[subrep[i][j]].pos[1],
				rep[i][j] + 1);
			p1 = segment[subrep[i][j]].pos[0] + 1;
			p2 = segment[subrep[i][j]].pos[1] + 1;
			len_subrepseq[j] = copyseq(subrepseq[j],
				 &repseq[rep[i][j]][p1 - coord[rep[i][j]][1]], p2 - p1 + 1);
			if(segment[subrep[i][j]].eq_pos[1] > 0)	{
				reverse_seq(subrepseq[j], len_subrepseq[j]);
			}
		}
/*	Compute K2d matrix	*/
		EvoDist = (double *) ckalloc(4 * pow1(index_subrep[i]) * sizeof(double));
		for(j = 0; j < index_subrep[i] - 1; j ++)	{
			for(k = j + 1; k < index_subrep[i]; k ++)	{
				m = numc(j, k);
				EvoDist[m] = Align2Dist(subrepseq[j], subrepseq[k],
						 len_subrepseq[j], len_subrepseq[k]);
			}
		}
		MaxEvoAve = EvoAve = nim = 0;
		AveEvoDist = (double *) ckalloc(index_subrep[i] * sizeof(double));
		for(j = 0; j < index_subrep[i]; j ++)	{
			for(k = 0; k < index_subrep[i]; k ++)	{
				if(j == k)	continue;
				m = numc(j, k);
				AveEvoDist[j] += EvoDist[m];
				nim ++;
			}
			AveEvoDist[j] /= (index_subrep[i] - 1);
			if(AveEvoDist[j] > MaxEvoAve)	{
				MaxEvoAve = AveEvoDist[j];
				MaxItem = j;
			}
			EvoAve += AveEvoDist[j];
		}
		EvoAve = (EvoAve - MaxEvoAve) / (index_subrep[i] - 1);
		DiffEvo = 0;
		for(j = 0; j < index_subrep[i]; j ++)	{
			if(j == MaxItem)	continue;
			DiffEvo += pow1(EvoAve - AveEvoDist[j]);
		}
		DiffEvo = sqrt(DiffEvo / (index_subrep[i] - 1));
		free((void *) AveEvoDist);

		fprintf(fp, "Sub-repeat s%d Total segments: %d\n", i, index_subrep[i]);
/*
		if(singleton_mark[subrep[i][MaxItem]] && MaxEvoAve > EvoAve + 3 * DiffEvo)	{
*/
		if(MaxEvoAve >= EvoAve + WIDTH * DiffEvo)	{
			fprintf(fp, "Ancester: %d %s\n", MaxItem + 1, repnames[MaxItem]);
			ancestor[i] = subrep[i][MaxItem];
		} else	{
			fprintf(fp, "Ancester: no significant hit\n");
		}
/*	Output K2d matrix	*/
		outputmatrix(fp, EvoDist, index_subrep[i], repnames);
/*	Neighbor-Join		*/
		stree = (int *) ckalloc(2 * index_subrep[i] * sizeof(int));
		dist = (double *) ckalloc(2 * index_subrep[i] * sizeof(double));
		NeighborJoin(index_subrep[i], EvoDist, stree, dist);
/*	Output tree		*/
		treeseq = (char *) ckalloc(index_subrep[i] * 1000 * sizeof(char));
		printNJtree(treeseq, stree, dist, index_subrep[i], repnames);
		fprintf(fp, "TREE\n%s\n", treeseq); 
		free((void *) EvoDist);
		free((void *) len_subrepseq);
		free((void *) dist);
		free((void *) stree);
		free((void *) treeseq);
		for(j = 0; j < index_subrep[i]; j ++)	{
			free((void *) subrepseq[j]);
			free((void *) repnames[j]);
		}
		free((void **) subrepseq);
		free((void **) repnames);
	}
	fclose(fp);

	sprintf(temp, "%s.anc", inpfile);
	fp = ckopen(temp, "w");
	fprintf(fp, "chromosome	chromStart	chromEnd	chromSize	Ancestors	AncestorChrom	AncestorStart	AncestorEnd	AncestorSize\n");
	for(i = 0; i < num_segment; i ++)	{
		k = ancestor[segment[i].eq_pos[0]];
		if(k >= 0)	{
		  if(segment[i].eq_pos[1] == segment[k].eq_pos[1])	{
			fprintf(fp, "%s	%d	%d	%d	s%d	%s	%d	%d	%d\n",
			  chrname[segment[i].chro], segment[i].pos[0] + 1, segment[i].pos[1] + 1, len_chroseq[segment[i].chro], 
			  segment[i].eq_pos[0], chrname[segment[k].chro], segment[k].pos[0] + 1,
			  segment[k].pos[1] + 1, segment[k].length);
		  } else	{
			fprintf(fp, "%s	%d	%d	%d	s%d	%s	%d	%d	%d\n",
			  chrname[segment[i].chro], segment[i].pos[0] + 1, segment[i].pos[1] + 1, len_chroseq[segment[i].chro], 
			  segment[i].eq_pos[0], chrname[segment[k].chro], segment[k].pos[1] + 1,
			  segment[k].pos[0] + 1, segment[k].length);
		  }
		} else	{
			fprintf(fp, "%s	%d	%d	%d	s%d	No hit\n",
			  chrname[segment[i].chro], segment[i].pos[0] + 1, segment[i].pos[1] + 1,
			  len_chroseq[segment[i].chro], segment[i].eq_pos[0]);
		}
	}
	fclose(fp);

	free((void *) singleton_mark);
	for(i = 0; i < MAX_SUBREP; i ++)	{
		free((void *) subrep[i]);
		free((void *) rep[i]);
	}
	free((void *) index_subrep);
	free((void *) ancestor);
	free((void **) subrep);
	free((void **) rep);
	for(i = 0; i < num_repeats; i ++)	{
		free((void *) repseq[i]);
		free((void *) coord[i]);
	}
	free((void **) repseq);
	free((void *) len_repseq);
	chrname = free_name(chrname, 100);
	free((void *) repeats);
	free((void *) segment);
}

int segcompar(SEGMENT *a, SEGMENT *b)
{
	if(a -> chro > b -> chro)	return(1);
	else if(a -> chro == b -> chro && a -> pos[0] > b -> pos[0])	return(1);
	else				return(-1);
}

void initenv(int argc, char **argv)
{
	int copt;
	int inpseq, outseq;
	extern char *optarg;

	WIDTH = 3;
	min_length = 5000;
	centro_dist = 5000000;
	inpseq = outseq = 0;

	while ((copt=getopt(argc,argv,"i:o:l:t:c:w:")) != EOF)	{
		switch(copt) {
			case 'i':
			  inpseq = 1;
			  sscanf(optarg,"%s", inpfile);
			  continue;
			case 'o':
			  outseq = 1;
			  sscanf(optarg,"%s", outfile);
			  continue;
			case 'L':
			  sscanf(optarg,"%d", &min_length);
			  continue;
			case 'l':
			  sscanf(optarg,"%s", lenfile);
			  continue;
			case 't':
			  sscanf(optarg,"%d", &tel_dist);
			  continue;
			case 'w':
			  sscanf(optarg,"%d", &WIDTH);
			  continue;
			case 'c':
			  sscanf(optarg,"%d", &centro_dist);
			  continue;
			default:
			  printf("prepare_aln -i InpFile -l LenFile -o outfile [-L min_leg -t Tel_dist -c Centro_dist -w WIDTH]\n");
			  printf("-i InpFile: The input file name of segments\n");
			  printf("-l lenfile: inpput file for chromosome length.\n");
			  printf("-o OutFile: output alignment files\n");
			  printf("-L min_length: minimum interval length between repeats\n");
			  printf("-t Tel_dist: distance from the telemere\n");
			  printf("-c Centro_dist: distance from the centromere\n");
			  printf("-w WIDTH: the width of the variance\n");
			  exit(-1);
		}
		optind--;
	}

	if(inpseq == 0 || outseq == 0)	{
		printf("prepare_aln -i InpFile -l LenFile -o outfile [-L min_leg -t Tel_dist -c Centro_dist -w WIDTH]\n");
		printf("-i InpFile: The input file name of segments\n");
		printf("-l lenfile: inpput file for chromosome length.\n");
		printf("-o OutFile: output alignment files\n");
		printf("-L min_length: minimum interval length between repeats\n");
		printf("-t Tel_dist: distance from the telemere\n");
		printf("-c Centro_dist: distance from the centromere\n");
		printf("-w WIDTH: the width of the variance\n");
		exit(-1);
	}
}
