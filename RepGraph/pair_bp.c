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

int min_length, min_seg_length;
int tel_dist, centro_dist;
int WIDTH;

char inpfile[100], lenfile[100], outfile[100];

void initenv(int argc, char **argv);
int segcompar(SEGMENT *a, SEGMENT *b);

main(int argc, char **argv)
{
	int	i, j, k, l, m, n, centro[200], len_chroseq[100];
	int	p1, p2;
	char	**chrname;
	int	num_chro;
	int	**rep;
	SEGMENT	*segment;
	int	num_segment;
	int	*repeats, num_repeats;
	int	*bp;
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
	printf("# segments: %d\n", num_segment);

	i = 0;
	while(i < num_segment)	{
		if(segment[i].src_pos[1] - segment[i].src_pos[0] + 1 < min_seg_length)	{
			segment[i] = segment[num_segment - 1];
			num_segment --;
		} else	 {
			i ++;
		}
	}
	printf("# long segments: %d\n", num_segment);

/*	sort the segments	*/
	qsort((void *) segment, num_segment, sizeof(SEGMENT), (void *) segcompar);

/*	Define repeats from sub-repeats	*/

	repeats = (int *) ckalloc(num_segment * sizeof(int));
	num_repeats = 0;
	for(i = 0; i < num_segment; i ++)	{
		if(i == num_segment - 1 || segment[i + 1].pos[0] > segment[i].pos[1] + min_length ||
		   segment[i + 1].chro != segment[i].chro)	{
			repeats[num_repeats ++] = i + 1;
		}
	}
	repeats[num_repeats ++] = num_segment;
	printf("%d repeats found.\n", num_repeats);

/*	get the symbolic representation of each repeat	*/
	rep = (int **) ckalloc(num_repeats * sizeof(int *));
	k = 0;
	for(i = 0; i < num_repeats; i ++)	{
		rep[i] = (int *) ckalloc((repeats[i] - k + 1) * sizeof(int));
		for(j = k; j < repeats[i]; j ++)	{
			if(segment[j].eq_pos[1])	{
				rep[i][j - k] = -(segment[j].eq_pos[0]);
			} else	{
				rep[i][j - k] = segment[j].eq_pos[0];
			}
		}
		k = repeats[i];
	}

/*	Compute pairwise # breakpoints	*/
	bp = (int *) ckalloc(num_repeats * (num_repeats - 1) / 2 * sizeof(int));
	k = 0;
	for(i = 0; i < num_repeats; i ++)	{
		for(j = i + 1; j < num_repeats; j ++)	{
			m = numc(i, j);
			bp[m] = comput_bp(rep[i], repeats[i] - k, rep[j], repeats[j] - repeats[j - 1]);
		}
		k = repeats[i];
	}

/*	Output ##breakspoint	*/
	fp = ckopen(outfile, "w");
	k = 0;
	for(i = 0; i < num_repeats; i ++)	{
		for(j = i + 1; j < num_repeats; j ++)	{
			m = numc(i, j);
			if(bp[m] > 0)	{
				fprintf(fp, "%s	%d	%d	%d	%s	%d	%d	%d	%d	%.4f\n",
					chrname[segment[k].chro], segment[k].pos[0], segment[repeats[i] - 1].pos[1], repeats[i] - k,
					chrname[segment[repeats[j] - 1].chro], segment[repeats[j - 1]].pos[0], segment[repeats[j] - 1].pos[1], repeats[j] - repeats[j - 1],
					bp[m], 1 - ((double) bp[m]) / (repeats[i] - k - 1));
			}
		}
		k = repeats[i];
	}
	fclose(fp);

	free((void *) bp);
	for(i = 0; i < num_repeats; i ++)	{
		free((void *) rep[i]);
	}
	free((void **) rep);
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
	min_length = 10000;
	centro_dist = 5000000;
	min_seg_length = 100;
	inpseq = outseq = 0;

	while ((copt=getopt(argc,argv,"i:o:l:g:s:")) != EOF)	{
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
			case 'g':
			  sscanf(optarg,"%d", &min_length);
			  continue;
			case 's':
			  sscanf(optarg,"%d", &min_seg_length);
			  continue;
			default:
			  printf("pair_bp -i InpFile -l LenFile -o outfile [-g min_gap_size -s min_seg_length\n");
			  printf("-i InpFile: The input file name of segments\n");
			  printf("-l lenfile: inpput file for chromosome length.\n");
			  printf("-o OutFile: output alignment files\n");
			  printf("-g min_gap_size: minimal gap size between segments (default: 10000)\n");
			  printf("-s min_seg_length: minimal length of segments (default 100)\n");
			  exit(-1);
		}
		optind--;
	}

	if(inpseq == 0 || outseq == 0)	{
		printf("pair_bp -i InpFile -l LenFile -o outfile [-g min_gap_size -s min_seg_length\n");
		printf("-i InpFile: The input file name of segments\n");
		printf("-l lenfile: inpput file for chromosome length.\n");
		printf("-o OutFile: output alignment files\n");
		printf("-g min_gap_size: minimal gap size between segments (default: 10000)\n");
		printf("-s min_seg_length: minimal length of segments (default 100)\n");
		exit(-1);
	}
}
