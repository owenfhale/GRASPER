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

char allcolor[50][100] = {
	"red", "blue", "green", "yellow", "orange",
	"azure", "bisque", "brown", "cyan", "firebrick",
	"gold", "gray", "khaki", "magenta", "maroon",
	"orange", "orchid", "pink", "purple", "salmon",
	"sienna", "tan", "thistle", "tomato", "wheat",
	"olivedrab", "ivory", "indianred", "hotpink", "lightskyblue",
	"orangered", "plum", "seashell", "steelblue", "turquoise",
	"coral", "cornsilk", "aquamarine", "mediumorchid", "navajowhite",
	"rosybrown", "seagreen", "slateblue", "springgreen", "violetred",
	"darkgoldenrod", "darkseagreen", "deeppink", "skyblue", "goldenrod"};
	

#define MIN_NUM 6
#define BIGGAP 5000

int min_length;

int segcompar(SEGMENT *a, SEGMENT *b);

main(int argc, char **argv)
{
	int	i, j, k, l, m, n, len_seq[100], reallen[100];
	int	k1, k2;
	char	**chrname;
	int	*dir, *breaks, *colorname, *colorindex, base[2];
	int	global_scale;
	SEGMENT	*segment, segment0;
	int	num_segment;
	int	*eqpos, *mark;
	int	range[2];
	int	*multip, *first, *begin, *end, *length;
	int	num_chro;
	char	name[1000], temp[1000], c;
	char	pt;
	FILE	*fp, *fp1;

	if(argc < 5)	{
		printf("Usage: makechrrepg intv_file output_file(ps format) chrolistfile min_leg\n");
		exit(-1);
	}
	min_length= atoi(argv[4]);

	chrname = alloc_name(100, 100);
	fp = ckopen(argv[3], "r");
	num_chro = readchrolist(len_seq, chrname, reallen, fp);
	fclose(fp);
	first = (int *) ckalloc(20000 * sizeof(int));
	multip = (int *) ckalloc(20000 * sizeof(int));
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

	length = (int *) ckalloc(2 * num_segment * sizeof(int));
	breaks = (int *) ckalloc(2 * num_segment * sizeof(int));
	colorname = (int *) ckalloc(20000 * sizeof(int));
	colorindex = (int *) ckalloc(2 * num_segment * sizeof(int));
	dir = (int *) ckalloc(2 * num_segment * sizeof(int));
	k = m = 0;
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
			multip[segment[i].eq_pos[0]] ++;
			dir[n] = segment[i].eq_pos[1];
			length[n] = segment[i].length;
			breaks[n ++] = segment[i].pos[1];
		}
	}
	k = 0;
	for(i = 0; i < 3000; i ++)	{
		if(multip[i] > 0)	{
			colorname[i] = k;
			k ++;
		}
	}
printf("total color %d\n", k);

	dir[n] = 0;
	breaks[n] = len_seq[segment[i - 1].chro];
	begin = (int *) ckalloc((1 + n) * sizeof(int));
	end = (int *) ckalloc((1 + n) * sizeof(int));

	k1 = k2 = 1;
	for(i = 0; i <= n; i ++)	{
		if(i == 0 || i == n || colorindex[i] == 0)	{
			begin[i] = -k2;
			k1 ++;
			end[i] = k1;
			k2 ++;
		} else	{
			begin[i] = k1;
			k1 ++;
			end[i] = k1;
		}
		if(multip[colorindex[i]] > 1 && first[colorindex[i]] == 0) {
			first[colorindex[i]] = i + 1;
		}
	}

	eqpos = (int *) ckalloc((k1 + 1) * sizeof(int));
	for(i = 1; i <= k1; i ++)	{
		eqpos[i] = i;
	}
	for(i = 0; i <= n; i ++)	{
		if(multip[colorindex[i]] > 1) {
			eqpos[begin[i]] = begin[first[colorindex[i]] - 1];
			eqpos[end[i]] = end[first[colorindex[i]] - 1];
		}
	}
	mark = (int *) ckalloc((k1 + 1) * sizeof(int));
	for(i = 1; i <= k1; i ++)	{
		if(begin[i] < 0)	continue;
		mark[eqpos[begin[i]]] = 1;
		mark[eqpos[end[i]]] = 1;
	}

	fp = ckopen(argv[2], "w");
	fprintf(fp, "digraph G {\n");
	fprintf(fp, "size=\"8,8\";\n");
	fprintf(fp, "\tgraph [rankdir=LR]\n");
	fprintf(fp, "\t{ node [shape=circle label=\"\"]\n\t\t");
	k = 1;
	for(j = 1; j < k1; j ++)	{
		if(mark[j] == 0)	continue;
		fprintf(fp, "%d ", j);
		if(k % 10 == 0)	{
			fprintf(fp, "\n\t\t");
		}
		k ++;
	}
	fprintf(fp, "\n\t}\n");
	fprintf(fp, "	{ node [shape=circle ]\n\t\t");
	for(j = 1; j < k2; j ++)	{
		fprintf(fp, "s%d ", j);
		if(j % 10 == 0)	{
			fprintf(fp, "\n\t\t");
		}
	}
	fprintf(fp, "\n\t\t");
	for(j = 1; j < k2; j ++)	{
		fprintf(fp, "e%d ", j);
		if(j % 10 == 0)	{
			fprintf(fp, "\n\t\t");
		}
	}
	l = 0;
	fprintf(fp, "\n\t}\n");
	k1 = k2 = 1;
	for(i = 0; i <= n; i ++)	{
		if(begin[i] < 0)	{
			if(i == 0)	{
			   fprintf(fp, "\ts%d -> %d [style=bold,color=black];\n", -begin[i], eqpos[end[i]]);
			} else if(i == n)	{
			   fprintf(fp, "\t%d -> e%d [style=bold,color=black];\n", eqpos[end[i]], -begin[i]);
			   k2 ++;
			} else	{
			   fprintf(fp, "\t%d -> e%d [style=bold,color=black];\n",
				 eqpos[end[i - 1]], -begin[i]);
			   fprintf(fp, "\ts%d -> %d [style=bold,color=black];\n", -begin[i], eqpos[end[i]]);
			}
			l ++;
		} else if(multip[colorindex[i]] == 1 || first[colorindex[i]] == i + 1)	{
			k = colorname[colorindex[i]] / 50;
			m = colorname[colorindex[i]] % 50;
			if(k == 0)	{
				sprintf(temp, "%s", allcolor[m]);
			} else	{
				sprintf(temp, "%s%d", allcolor[m], k);
			}
			fprintf(fp, "\t%d -> %d [style=bold,color=%s,label=\"s%d(%d,%d)\"];\n",
				eqpos[begin[i]], eqpos[end[i]], temp,
				 colorindex[i], length[i], multip[colorindex[i]]);
			k1 ++;
		}
	}
	fprintf(fp, "}\n");
	fclose(fp);
	printf("%d regions output.\n", l - 1);

	chrname = free_name(chrname, 100);
	free((void *) breaks);
	free((void *) colorname);
	free((void *) colorindex);
	free((void *) dir);
	free((void *) segment);
	free((void *) first);
	free((void *) multip);
	free((void *) length);
	free((void *) mark);
	free((void *) eqpos);
	free((void *) begin);
	free((void *) end);
}

int segcompar(SEGMENT *a, SEGMENT *b)
{
	if(a -> chro > b -> chro)	return(1);
	else if(a -> chro == b -> chro && a -> pos[0] > b -> pos[0])	return(1);
	else				return(-1);
}
