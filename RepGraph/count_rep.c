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

char lenfile[100], outfile[100];
char chriname[100];

int min_length, nchro;

void initenv(int argc, char **argv);
int segcompar(SEGMENT *a, SEGMENT *b);

main(int argc, char **argv)
{
	int	i, j, k, l, m, n, centro[200], len_chroseq[100];
	int	maxi;
	int	num_seq;
	char	**chrname, **repnames;
	int	num_chro;
	int	*counts, maxc;
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

	nchro = findgenname(chriname, chrname, num_chro);

/*	input segments	*/

	segment = (SEGMENT *) ckalloc(60000 * sizeof(SEGMENT));
	fp = stdin;
	num_segment = input_segment(segment, chrname, num_chro, fp);

/*	sort the segments	*/
/*
	qsort((void *) segment, num_segment, sizeof(SEGMENT), (void *) segcompar);
*/

/*	Define repeats from sub-repeats	*/

	repeats = (int *) ckalloc(num_segment * sizeof(int));
	num_repeats = 0;
	counts = (int *) ckalloc(10000 * num_segment * sizeof(int));
	n = m = 0;
	j = 0;

	for(i = 0; i < num_segment; i ++)	{
		if(i == num_segment - 1 || segment[i + 1].pos[0] > segment[i].pos[1] + min_length ||
		   segment[i + 1].chro != segment[i].chro)	{
			repeats[num_repeats ++] = i;
			if(segment[i].chro == nchro)	j ++;
		}
		if(segment[i].pos[1] - segment[i].pos[0] < min_length)	continue;
		counts[segment[i].eq_pos[0]] ++;
		if(segment[i].eq_pos[0] > n)	n = segment[i].eq_pos[0];
	}
	m += n;
	k = maxc = 0;
	for(i = 0; i < m; i ++)	{
		if(counts[i] > 1)	k ++;
		if(counts[i] > maxc)	{
			maxi = i;
			maxc = counts[i];
		}
	}
printf("m %d k %d\n", m, k);
	printf("%d repeats (%s %d) %d subrepeats %d max_multip subrepeats index %d.\n", num_repeats,
		chriname, j, k, maxc, maxi);

	free((void *) counts);
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
	extern char *optarg;

	strcpy(chriname, "chr1");

	min_length = 600;

	while ((copt=getopt(argc,argv,"i:l:s:")) != EOF)	{
		switch(copt) {
			case 'l':
			  sscanf(optarg,"%s", lenfile);
			  continue;
			case 's':
			  sscanf(optarg,"%s", chriname);
			  continue;
			default:
			  printf("count_rep -l LenFile\n");
			  printf("-l lenfile: inpput file for chromosome length.\n");
			  exit(-1);
		}
		optind--;
	}
}
