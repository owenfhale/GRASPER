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

int readreputer(ALIGN **align, char **src_seq, int *len_seq, FILE *fp, int MIN_LEG, double MIN_ID);
ALIGN *new_align(int *sapp, ALIGN *align, int r1, int r2, int pos1, int pos2, int len1, int len2);
double cal_identity(char *seq1, int len1, char *seq2, int len2, int *sapp);

int readreputer(ALIGN **align, char **src_seq, int *len_seq, FILE *fp, int MIN_LEG, double MIN_ID)
{
	int	i, j, k, l, m, n, s1, s2;
	int	len1, len2, pos1, pos2;
	int	nall;
	int	*sapp;
	char	dir[10];
	double	id;
	char	str[500];
	char	**src_name;

	src_name = (char **) ckalloc(2 * sizeof(char *));
	for(i = 0; i < 2; i ++)	{
		src_name[i] = (char *) ckalloc(100 * sizeof(char));
	}
	strcpy(src_name[0], "Seq_forward");
	strcpy(src_name[1], "Seq_reverse");
	sapp = (int *) ckalloc(len_seq[0] * sizeof(int));

	nall = 0;
	while(fgets(str, 400, fp))	{
		if(str[0] != '#')	{
			sscanf(str, "%d%d%s%d%d", &len1, &pos1, dir, &len2, &pos2);
			if(len1 < MIN_LEG || len2 < MIN_LEG)	continue;
			if(dir[0] == 'F')	{
				s1 = s2 = 0;
				band += abs(len1 - len2);
				k = ALIGN0(&src_seq[0][pos1 - 1], &src_seq[0][pos2 - 1], len1, len2,
					  -band, band, W, g, h, sapp, len1 + len2, len1 + len2);
				id = cal_identity(&src_seq[0][pos1], len1, &src_seq[0][pos2], len2, sapp);
				if(id >= MIN_ID)	{
					align[0] = new_align(sapp, align[0], s1, s2,
						    pos1, pos2, len1, len2);
/*
if(id < 0.60)	{
	printf("%d %d %c %d %d %.4f\n", pos1, pos1 + len1 - 1, dir[0], pos2, pos2 + len2 - 1, id);
	output_align(align[0], src_name, src_seq, len_seq, 2);
}
*/
					nall ++;
				}
			} else	{
				s1 = 0;
				s2 = 1;
				pos2 = len_seq[0] - pos2 - len2;
				k = ALIGN0(&src_seq[0][pos1 - 1], &src_seq[1][pos2 - 1], len1, len2,
					  -band, band, W, g, h, sapp, len1 + len2, len1 + len2);
				id = cal_identity(&src_seq[0][pos1], len1, &src_seq[1][pos2], len2, sapp);
				if(id >= MIN_ID)	{
					align[0] = new_align(sapp, align[0], s1, s2,
						    pos1, pos2, len1, len2);
/*
if(id < 0.60)	{
	printf("%d %d %c %d %d %.4f\n", pos1, pos1 + len1 - 1, dir[0], pos2, pos2 + len2 - 1, id);
	output_align(align[0], src_name, src_seq, len_seq, 2);
}
*/
					nall ++;
				}
			}
		}
	}
	free((void *) sapp);
	for(i = 0; i < 2; i ++)	{
		free((void *) src_name[i]);
	}
	free((void **) src_name);
	return(nall);
}

double cal_identity(char *seq1, int len1, char *seq2, int len2, int *sapp)
{
	int	i, j, k, l, m, n;
	int	op;
	double	p, s;

	s = 0;
	i = j = k = l = 0;
	while(i < len1 || j < len2)	{
		op = *sapp ++;
		if(op == 0)	{
			if(seq1[i] == seq2[j])	{
				k ++;
			} else	{
				l ++;
			}
			n ++;
			i ++;
			j ++;
		} else if(op > 0)	{
			l += op;
			j += op;
		} else 	{
			l -= op;
			i -= op;
		}
	}
	p = ((double) k) / (k + l);
	return(p);
}

ALIGN *new_align(int *sapp, ALIGN *align, int r1, int r2, int pos1, int pos2, int len1, int len2)
{
	int	i, j, k, l, m, n, len, mark;
	int	op;
	int	mis_match;
	int	*pos[2];
	ALIGN	*newalign, *aln, *aln_last;

	newalign = (ALIGN *) ckalloc(1 * sizeof(ALIGN));
	newalign -> reads[0] = r1;
	newalign -> reads[1] = r2;
	len = len1 + len2;
	for(i = 0; i < 2; i ++)	{
		pos[i] = (int *) ckalloc(len * sizeof(int));
	}
	i = j = k = mark = 0;
	mis_match = 0;
	while(i < len1 && j < len2)	{
		op = *sapp ++;
		if(op == 0)	{
			if(mark == 0)	{
				mark = 1;
				pos[0][k] = pos1 + i;
				pos[1][k] = pos2 + j;
				k ++;
			}
			i ++;
			j ++;
		} else if(op > 0)	{
			mark = 0;
			j += op;
		} else 	{
			mark = 0;
			i -= op;
		}
	}
	pos[0][k] = pos1 + i;
	pos[1][k] = pos2 + j;
	k ++;
	newalign -> length = k;
	newalign -> mis_match = mis_match;
	newalign -> pos[0] = (int *) ckalloc(k * sizeof(int));
	newalign -> pos[1] = (int *) ckalloc(k * sizeof(int));
	for(i = 0; i < k; i ++)	{
		newalign -> pos[0][i] = pos[0][i];
		newalign -> pos[1][i] = pos[1][i];
	}
	newalign -> next = align;
	align = newalign;

	free((void *) pos[0]);
	free((void *) pos[1]);
	return(align);
}
