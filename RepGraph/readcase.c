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
#include <extvab.h>
#include <extfunc.h>

#define MAX_ITEM 100
#define MAX_ITEM_LEG 200

int readcase(ALIGN **align, char **src_seq, int *len_seq, int num_seq, FILE *fp, int MIN_LEG, double MIN_ID);
void extract_filename(char *item, char *filename);
void getaln(char *seq1, char *seq2, int len1, int len2, int *sapp, char *filename);

int readcase(ALIGN **align, char **src_seq, int *len_seq, int num_seq, FILE *fp, int MIN_LEG, double MIN_ID)
{
	int	i, j, k, l, m, n, s1, s2;
	int	len1, len2, pos1, pos2;
	int	nall;
	int	*sapp;
	char	**item;
	char	dir[10], filename[100];
	char	chrname[2][200];
	int	startpos[2], endpos[2];
	double	id;
	char	str[3000];
	char	**src_name;

	src_name = (char **) ckalloc(2 * sizeof(char *));
	for(i = 0; i < 2; i ++)	{
		src_name[i] = (char *) ckalloc(100 * sizeof(char));
	}
	strcpy(src_name[0], "Seq_forward");
	strcpy(src_name[1], "Seq_reverse");
	item = (char **) ckalloc(MAX_ITEM * sizeof(char *));
	for(i = 0; i < MAX_ITEM; i ++)	{
		item[i] = (char *) ckalloc(MAX_ITEM_LEG * sizeof(char));
	}
	sapp = (int *) ckalloc(2 * len_seq[0] * sizeof(int));

	nall = 0;
	while(fgets(str, 2990, fp))	{
		if(str[0] == '#')	continue;
		n = parsestr(str, item, "\t");
		strcpy(chrname[0], item[0]);
		startpos[0] = atoi(item[1]);
		endpos[0] = atoi(item[2]);
		strcpy(chrname[1], item[4]);
		startpos[1] = atoi(item[5]);
		endpos[1] = atoi(item[6]);
		extract_filename(item[17], filename);
/*
		s1 = getchr(chrname[0]);
		s2 = getchr(chrname[1]);
*/
		s1 = s2 = 0;
		if(startpos[0] > endpos[0] && startpos[1] > endpos[1])	{
			k = startpos[0];
			startpos[0] = endpos[0];
			endpos[0] = k;
			k = startpos[1];
			startpos[1] = endpos[1];
			endpos[1] = k;
		} else if(startpos[0] > endpos[0])	{
			k = startpos[0];
			startpos[0] = endpos[0];
			endpos[0] = k;
			k = startpos[1];
			startpos[1] = len_seq[s2] - endpos[1] - 1;
			endpos[1] = len_seq[s2] - startpos[1] - 1;
			s2 += num_seq;
		} else if(startpos[1] > endpos[1])	{
			startpos[1] = len_seq[s2] - startpos[1] - 1;
			endpos[1] = len_seq[s2] - endpos[1] - 1;
			s2 += num_seq;
		}
		len1 = endpos[0] - startpos[0] + 1;
		len2 = endpos[1] - startpos[1] + 1;
		if(len1 < MIN_LEG || len2 < MIN_LEG)	continue;
printf("filename %s\n", filename);
		getaln(&src_seq[s1][startpos[0]], &src_seq[s2][startpos[1]], len1, len2, sapp, filename);
		id = cal_identity(&src_seq[s1][startpos[0]], len1, &src_seq[s2][startpos[0]], len2, sapp);
printf("id1 %f\n", id);
		if(id >= MIN_ID)	{
			align[0] = new_align(sapp, align[s1], s1, s2, startpos[0], startpos[1], len1, len2);
/*
if(id < 0.60)	{
	printf("%d %d %c %d %d %.4f\n", startpos[0], endpos[1], 'F', startpos[1], endpos[1], id);
	output_align(align[0], src_name, src_seq, len_seq, num_seq);
}
*/
			nall ++;
		}
	}
	free((void *) sapp);
	for(i = 0; i < MAX_ITEM; i ++)	{
		free((void *) item[i]);
	}
	free((void *) item);
	for(i = 0; i < 2; i ++)	{
		free((void *) src_name[i]);
	}
	free((void **) src_name);
	return(nall);
}

void extract_filename(char *item, char *filename)
{
	int	i, j, k, l;

	l = strlen(item);
	for(i = l - 1; i >= 0; i --)	{
		if(item[i] == '/')	{
			break;
		}
	}
	sprintf(filename, "data/%s.indel", &item[i + 1]);
}

void getaln(char *seq1, char *seq2, int len1, int len2, int *sapp, char *filename)
{
	int	i, j, k, l;
	int	pos11, pos12, pos21, pos22;
	FILE	*fp;
	char	str[500];

	i = j = k = 0;
	fp = ckopen(filename, "r");
	while(fgets(str, 490, fp))	{
		if(str[0] == 'a')	continue;
		sscanf(str, "%*s%*s%*s%*s%d%d%d%d", &pos11, &pos12, &pos21, &pos22);
		while(i < pos11 && j < pos21)	{
			k ++;
			i ++;
			j ++; 
		}
		if(i == pos11)	{
			j = pos22;
			sapp[k ++] = pos22 - pos21 + 1;	
		} else if(j == pos21)	{
			i = pos12;
			sapp[k ++] = -(pos12 - pos11 + 1);
		}
	}
	fclose(fp);
	while(i < len1 && j < len2)	{
		k ++;
		i ++;
		j ++; 
	}
	if(i != len1 || j != len2)	{
		printf("File %s len1 %d len2 %d i %d j %d\n", filename, len1, len2, i, j);
		exit(0);
	}
}
