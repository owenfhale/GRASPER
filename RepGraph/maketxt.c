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

int min_length = 0;

int segcompar(SEGMENT *a, SEGMENT *b);

main(int argc, char **argv)
{
	int	i, j, k, l, m, n, len_seq[100], reallen[100];
	char	**chrname;
	int	global_scale;
	int	*insertlen, *insertpos, num_ins;
	SEGMENT	*segment, *insertreg;
	long int START_POS;
	int	num_segment, num_insertreg;
	int	**repeats;
	int	range[2];
	int	num_chro;
	char	name[1000], temp[1000], c, str[500], newstr[500];
	char	*label, pt;
	FILE	*fp;

	if(argc < 5)	{
		printf("Usage: maketxt intv_file uniq_reg_file txt_output chrolistfile\n");
		exit(-1);
	}

	chrname = alloc_name(100, 100);
	fp = ckopen(argv[4], "r");
	num_chro = readchrolist(len_seq, chrname, reallen, fp);
	fclose(fp);
	printf("num_chro %d\n", num_chro);
	segment = (SEGMENT *) ckalloc(MAX_INTV * sizeof(SEGMENT));
	fp = ckopen(argv[1], "r");
	num_segment = input_intevals(fp, segment, len_seq, num_chro);
	fclose(fp);
	printf("num_segment %d\n", num_segment);
	insertreg = (SEGMENT *) ckalloc(MAX_INTV * sizeof(SEGMENT));
	fp = ckopen(argv[2], "r");
	num_insertreg = input_reg(fp, insertreg);
	fclose(fp);
	printf("num_insertreg %d\n", num_insertreg);

	transform_segment(segment, num_segment, insertreg, num_insertreg);
	qsort((void *) segment, num_segment, sizeof(SEGMENT), (void *) segcompar);

/*	Statistics about contingency of sub-repeats	*/
	sprintf(temp, "%s.pair", argv[4]);
	fp = ckopen(temp, "w");
	fprintf(fp, "Subrepeat1	Subrepeat2	Chromosome	Subrepeat1 (From-to)	Subrepeat2 (From-to)	Gap\n");
	for(i = 0; i < num_segment - 1; i ++)	{
		if(segment[i + 1].chro != segment[i].chro || segment[i + 1].pos[0] > segment[i].pos[1] + BIGGAP)	{
			continue;
		}
		if(segment[i].eq_pos[1] == 0)	{
			fprintf(fp, "s%d	", segment[i].eq_pos[0]);
			fprintf(fp, "s%d	", segment[i + 1].eq_pos[0]);
			fprintf(fp, "%s	", chrname[segment[i + 1].chro]);
			fprintf(fp, "%d-%d	", segment[i].pos[0] + 1, segment[i].pos[1] + 1);
			if(segment[i + 1].eq_pos[1] == 0)	{
				fprintf(fp, "%d-%d	", segment[i + 1].pos[0] + 1, segment[i + 1].pos[1] + 1);
			} else	{
				fprintf(fp, "%d-%d	", segment[i + 1].pos[1] + 1, segment[i + 1].pos[0] + 1);
			}
			if(segment[i + 1].pos[0] <= segment[i].pos[1] + 1)	{
				fprintf(fp, "0\n");	
			} else if(segment[i + 1].pos[0] > segment[i].pos[1])	{
				fprintf(fp, "%d\n", segment[i + 1].pos[0] - segment[i].pos[1] - 1);	
			}
		} else	{
			fprintf(fp, "s%d	", segment[i + 1].eq_pos[0]);
			fprintf(fp, "s%d	", segment[i].eq_pos[0]);
			fprintf(fp, "%s	", chrname[segment[i + 1].chro]);
			if(segment[i + 1].eq_pos[1] == 0)	{
				fprintf(fp, "%d-%d	", segment[i + 1].pos[0] + 1, segment[i + 1].pos[1] + 1);
			} else	{
				fprintf(fp, "%d-%d	", segment[i + 1].pos[1] + 1, segment[i + 1].pos[0] + 1);
			}
			fprintf(fp, "%d-%d	", segment[i].pos[1] + 1, segment[i].pos[0] + 1);
			if(segment[i + 1].pos[0] <= segment[i].pos[1] + 1)	{
				fprintf(fp, "0\n");	
			} else if(segment[i + 1].pos[0] > segment[i].pos[1])	{
				fprintf(fp, "%d\n", segment[i + 1].pos[0] - segment[i].pos[1] - 1);	
			}
		}
	}
	fclose(fp);

/*	Output subrepeats	*/
	sprintf(temp, "%s", argv[3]);
	fp = ckopen(temp, "w");
	fprintf(fp, "chromosome	chromStart	chromEnd	chromSize	Repeats	RepeatStart	RepeatEnd	RepeatSize\n");
	for(i = 0; i < num_segment; i ++)	{
		if(segment[i].eq_pos[1] == 0)	{
			fprintf(fp, "%s	%d	%d	%d	s%d	%d	%d	%d\n",
			  chrname[segment[i].chro], segment[i].pos[0] + 1, segment[i].pos[1] + 1, reallen[segment[i].chro], 
			  segment[i].eq_pos[0], segment[i].src_pos[0], segment[i].src_pos[1], segment[i].length);
		} else	{
			fprintf(fp, "%s	%d	%d	%d	s%d	%d	%d	%d\n",
			  chrname[segment[i].chro], segment[i].pos[0] + 1, segment[i].pos[1] + 1, reallen[segment[i].chro], 
			  segment[i].eq_pos[0], segment[i].src_pos[1], segment[i].src_pos[0], segment[i].length);
		}
	}
	fclose(fp);

	chrname = free_name(chrname, 100);
	free((void *) segment);
	free((void *) insertreg);
}

int segcompar(SEGMENT *a, SEGMENT *b)
{
	if(a -> chro > b -> chro)	return(1);
	else if(a -> chro == b -> chro && a -> pos[0] > b -> pos[0])	return(1);
	else				return(-1);
}
