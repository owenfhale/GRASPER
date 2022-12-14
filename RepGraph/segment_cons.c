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
#define MIN_SEG_LEG 1000
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
	char	*conseq;
	int	len_conseq;
	int	MaxItem, nim;
	char	**chrname;
	int	num_chro;
	int	**subrep, **rep, *index_subrep, num_subrep;
	char	**repseq;
	int	*len_repseq;
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

/*	input segments	*/

	segment = (SEGMENT *) ckalloc(70000 * sizeof(SEGMENT));
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
	max_sub += 1000;
	num_subrep ++;
	printf("# subrepeats: %d\n", num_subrep);

/*	make consensus of each sub-repeat	*/

	fp = ckopen(outfile, "w");
	for(i = 0; i < num_subrep; i ++)	{
		if(index_subrep[i] <= 1)	continue;
		p1 = segment[subrep[i][0]].pos[0];
		p2 = segment[subrep[i][0]].pos[1];
		if(p2 - p1 < MIN_SEG_LEG)	{
			continue;
		}
		len_subrepseq = (int *) ckalloc(index_subrep[i] * sizeof(int));
		subrepseq = (char **) ckalloc(index_subrep[i] * sizeof(char *));
		for(j = 0; j < index_subrep[i]; j ++)	{
			subrepseq[j] = (char *) ckalloc(max_sub * sizeof(char));
		}
/*	Build subrepeat sequences	*/
		for(j = 0; j < index_subrep[i]; j ++)	{
			p1 = segment[subrep[i][j]].pos[0] + 1;
			p2 = segment[subrep[i][j]].pos[1] + 1;
			len_subrepseq[j] = copyseq(subrepseq[j],
				 &repseq[rep[i][j]][p1 - coord[rep[i][j]][1]], p2 - p1 + 1);
			if(segment[subrep[i][j]].eq_pos[1] > 0)	{
				reverse_seq(subrepseq[j], len_subrepseq[j]);
			}
		}
/*	Make consensus of each subrep	*/
		conseq = (char *) ckalloc(max_sub * sizeof(char));
		len_conseq = buildcons(subrepseq, len_subrepseq, index_subrep[i], conseq);
		fprintf(fp, ">subrepeat s%d %d\n", i, len_conseq);
		for(j = 0; j < len_conseq; j ++)	{
			fprintf(fp, "%c", na_name[conseq[j]]);
			if(j % 60 == 59)	{
				fprintf(fp, "\n");
			}
		}
		if(j % 60 != 0)	{
			fprintf(fp, "\n");
		}
		free((void *) conseq);
		free((void *) len_subrepseq);
		for(j = 0; j < index_subrep[i]; j ++)	{
			free((void *) subrepseq[j]);
		}
		free((void **) subrepseq);
	}
	fclose(fp);

	for(i = 0; i < MAX_SUBREP; i ++)	{
		free((void *) subrep[i]);
		free((void *) rep[i]);
	}
	free((void *) index_subrep);
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

	while ((copt=getopt(argc,argv,"i:o:l:")) != EOF)	{
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
			  sscanf(optarg,"%s", lenfile);
			  continue;
			default:
			  printf("segment_cons -i InpFile -l LenFile -o outfile\n");
			  printf("-i InpFile: The input file name of segments\n");
			  printf("-l lenfile: inpput file for chromosome length.\n");
			  printf("-o OutFile: output alignment files\n");
			  exit(-1);
		}
		optind--;
	}

	if(inpseq == 0 || outseq == 0)	{
		printf("segment_cons -i InpFile -l LenFile -o outfile\n");
		printf("-i InpFile: The input file name of segments\n");
		printf("-l lenfile: inpput file for chromosome length.\n");
		printf("-o OutFile: output alignment files\n");
		exit(-1);
	}
}
