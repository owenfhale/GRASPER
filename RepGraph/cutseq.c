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
#include <extfunc.h>

#define BIGGAP 5000
#define MAX_INTV 70000

int min_length = 100;

int segcompar(SEGMENT *a, SEGMENT *b);

main(int argc, char **argv)
{
	int	i, j, k, l, m, n, len_seq[100], reallen[100];
	char	**chrname;
	char	*seq;
	int	length;
	SEGMENT	*segment, *insertreg;
	int	num_segment, num_insertreg;
	int	num_chro;
	int	coord[3], dir;
	char	name[1000], temp[1000], c, str[500], newstr[500];
	FILE	*fp;

	if(argc < 5)	{
		printf("Usage: cutseq list hgs_directory outfile chrlist\n");
		exit(-1);
	}

	chrname = alloc_name(100, 100);
	fp = ckopen(argv[4], "r");
	num_chro = readchrolist(len_seq, chrname, reallen, fp);
	fclose(fp);
	segment = (SEGMENT *) ckalloc(MAX_INTV * sizeof(SEGMENT));
	fp = ckopen(argv[1], "r");
	num_segment = input_segment(segment, chrname, num_chro, fp);
	fclose(fp);
	printf("num_segment %d\n", num_segment);

	seq = (char *) ckalloc(30000 * sizeof(char));
	fp = ckopen(argv[3], "w");
	for(i = 0; i < num_segment; i ++)	{
		//sprintf(name, "%s/%s.fa", argv[2], chrname[segment[i].chro]);
		coord[0] = segment[i].chro;
		if(segment[i].eq_pos[1] == 0)	{
			coord[1] = segment[i].pos[0];
			coord[2] = segment[i].pos[1];
		} else	{
			coord[1] = segment[i].pos[1];
			coord[2] = segment[i].pos[0];
		}
		length = getchseq(seq, coord, segment[i].eq_pos[1], chrname);
printf("length %d\n", length);
		sprintf(str, "%s %d %d %d", chrname[segment[i].chro],
			coord[1], coord[2], segment[i].eq_pos[1]);
		writeseq(fp, seq, str, length);
	}
	fclose(fp);

	free((void *) seq);
	free((void *) segment);
}
