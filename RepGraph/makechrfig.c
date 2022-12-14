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
#include <math.h>
#include <perdef.h>
#include <param.h>
#include <extvab.h>
#include <extfunc.h>

#define MIN_NUM 6
#define BIGGAP 5000

int min_length;

int npages = 1;

double	scale[2];

int segcompar(SEGMENT *a, SEGMENT *b);

main(int argc, char **argv)
{
	int	i, j, k, l, m, n, len_seq[100], reallen[100];
	char	**chrname;
	int	*dir, *breaks, *colorindex, base[2];
	int	global_scale;
	SEGMENT	*segment, segment0;
	int	num_segment;
	int	range[2];
	int	num_chro;
	char	name[1000], temp[1000], c;
	char	pt;
	FILE	*fp, *fp1;

	if(argc < 5)	{
		printf("Usage: makechrfig intv_file output_file(ps format) chrolistfile min_leg\n");
		exit(-1);
	}
	min_length= atoi(argv[4]);

	chrname = alloc_name(100, 100);
	fp = ckopen(argv[3], "r");
	num_chro = readchrolist(len_seq, chrname, reallen, fp);
	fclose(fp);
	segment = (SEGMENT *) ckalloc(20000 * sizeof(SEGMENT));
	fp = ckopen(argv[1], "r");
	num_segment = input_segment(segment, chrname, num_chro, fp);
	fclose(fp);

/*	Sort the segments	*/
	qsort((void *) segment, num_segment, sizeof(SEGMENT), (void *) segcompar);
/*
for(i = 0; i < num_segment; i ++)	{
	printf("aftersort segment %d %d %d\n", segment[i].pos[0], segment[i].pos[1], segment[i].eq_pos[1]);
}
getchar();
*/

	global_scale = 25;
	breaks = (int *) ckalloc(2 * num_segment * sizeof(int));
	colorindex = (int *) ckalloc(2 * num_segment * sizeof(int));
	dir = (int *) ckalloc(2 * num_segment * sizeof(int));
	k = m = 0;
	fp = ckopen(argv[2], "w");
	prt_Header(fp, "Tang Hai-xu", "OUTPUT", "2-21-2004", 0);
	prt_Macro(fp);
	scale[0] = 0.7;
	scale[1] = 0.8;
	change_scale(fp, scale);
	base[0] = 50;
	base[1] = 850;
	range[0] = 1;
	range[1] = len_seq[segment[0].chro];
printf("range %d %d\n", range[0], range[1]);
	n = 0;
	breaks[n] = 0;
	for(i = 0; i < num_segment; i ++)	{
		if(segment[i].pos[1] - segment[i].pos[0] > 1000)	{
			if(segment[i].pos[0] - breaks[n - 1] > min_length)	{
				colorindex[n] = 0;
				dir[n] = 0;
				breaks[n ++] = segment[i].pos[0];
			}
			colorindex[n] = segment[i].eq_pos[0];
			dir[n] = segment[i].eq_pos[1];
			breaks[n ++] = segment[i].pos[1];
		}
	}
	dir[n] = 0;
	breaks[n ++] = len_seq[segment[i - 1].chro];
	l = output_region(fp, base, name, chrname[segment[0].chro], range,
		 global_scale, breaks, n, colorindex, dir);
	showpage(fp);
	page_trailer(fp);
	restore_scale(fp);
	prt_Trailer(fp, k + 1);
	fclose(fp);
	printf("%d regions output.\n", n - 1);

	chrname = free_name(chrname, 100);
	free((void *) breaks);
	free((void *) colorindex);
	free((void *) dir);
	free((void *) segment);
}

int segcompar(SEGMENT *a, SEGMENT *b)
{
	if(a -> chro > b -> chro)	return(1);
	else if(a -> chro == b -> chro && a -> pos[0] > b -> pos[0])	return(1);
	else				return(-1);
}
