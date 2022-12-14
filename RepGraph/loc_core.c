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
#define MAX_SUBREP 10000
#define MIN_SEG 2

int min_length, MIN_COPY;
int tel_dist, centro_dist;
int WIDTH;

char inpfile[100], lenfile[100], outfile[100];

void initenv(int argc, char **argv);
int segcompar(SEGMENT *a, SEGMENT *b);

main(int argc, char **argv)
{
	int	i, j, k, l, m, n, centro[200], len_chroseq[100];
	int	max_sub;
	int	num_seq;
	int	p1, p2;
	char	**chrname;
	int	num_chro;
	SEGMENT	*segment;
	int	num_segment;
	int	*repcopy;
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

	segment = (SEGMENT *) ckalloc(80000 * sizeof(SEGMENT));
	fp = ckopen(inpfile, "r");
	num_segment = input_segment(segment, chrname, num_chro, fp);
	fclose(fp);
	printf("num_segment %d\n", num_segment);

/*	sort the segments	*/
	qsort((void *) segment, num_segment, sizeof(SEGMENT), (void *) segcompar);

/*	Compute the copy number of subrepeats	*/
	m = 0;
	for(i = 0; i < num_segment; i ++)	{
		if(m < segment[i].eq_pos[0])	{
			m = segment[i].eq_pos[0];
		}
	}
	printf("Max index %d\n", m);

	repcopy = (int *) ckalloc((m + 1) * sizeof(int));
	for(i = 0; i < num_segment; i ++)	{
		repcopy[segment[i].eq_pos[0]] ++;
	}

/*	Define repeats from sub-repeats	*/

	repeats = (int *) ckalloc(num_segment * sizeof(int));
	num_repeats = 0;
	m = 0;
	for(i = 0; i < num_segment; i ++)	{
		if(m < repcopy[segment[i].eq_pos[0]])	m = repcopy[segment[i].eq_pos[0]];
		if(i < num_segment - 1 && (segment[i + 1].pos[0] > segment[i].pos[1] + min_length ||
		   segment[i + 1].chro != segment[i].chro || (repcopy[segment[i + 1].eq_pos[0]] >= MIN_COPY
		   && repcopy[segment[i].eq_pos[0]] < MIN_COPY)))	{
			repeats[i] = 1;
			num_repeats ++;
		}
	}
	printf("Max copy %d\n", m);
	printf("%d repeats found.\n", num_repeats);

/*	Output subrepeats	*/
	fp = ckopen(outfile, "w");
	fprintf(fp, "chromosome	chromStart	chromEnd	chromSize	Repeats	RepeatStart	RepeatEnd	RepeatSize\n");
	for(i = 0; i < num_segment; i ++)	{
		if(repcopy[segment[i + 1].eq_pos[0]] < MIN_COPY)	continue;
		if(segment[i].eq_pos[1] == 0)	{
			fprintf(fp, "%s	%d	%d	%d	s%d	%d	%d	%d\n",
			  chrname[segment[i].chro], segment[i].pos[0] + 1, segment[i].pos[1] + 1, len_chroseq[segment[i].chro], 
			  segment[i].eq_pos[0], segment[i].src_pos[0], segment[i].src_pos[1], segment[i].length);
		} else	{
			fprintf(fp, "%s	%d	%d	%d	s%d	%d	%d	%d\n",
			  chrname[segment[i].chro], segment[i].pos[0] + 1, segment[i].pos[1] + 1, len_chroseq[segment[i].chro], 
			  segment[i].eq_pos[0], segment[i].src_pos[1], segment[i].src_pos[0], segment[i].length);
		}
		if(repeats[i] > 0)	{
			fprintf(fp, "------------------------------------------------------\n");
		}
	}
	fclose(fp);

	free((void *) repcopy);
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
	MIN_COPY = 10;
	centro_dist = 5000000;
	inpseq = outseq = 0;

	while ((copt=getopt(argc,argv,"i:o:l:c:")) != EOF)	{
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
			  sscanf(optarg,"%d", &MIN_COPY);
			  continue;
			case 'l':
			  sscanf(optarg,"%s", lenfile);
			  continue;
			default:
			  printf("loc_core -i InpFile -l LenFile -o outfile\n");
			  printf("-i InpFile: The input file name of segments\n");
			  printf("-l lenfile: inpput file for chromosome length.\n");
			  printf("-o OutFile: output alignment files\n");
			  exit(-1);
		}
		optind--;
	}

	if(inpseq == 0 || outseq == 0)	{
		printf("loc_core -i InpFile -l LenFile -o outfile\n");
		printf("-i InpFile: The input file name of segments\n");
		printf("-l lenfile: inpput file for chromosome length.\n");
		printf("-o OutFile: output alignment files\n");
		exit(-1);
	}
}
