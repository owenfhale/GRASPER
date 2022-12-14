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

int buildcons(char **subrepseq, int *len_subrepseq, int num_subrep, char *conseq);
void seq2ALIGN(char *seq1, int len_seq1, char *seq2, int len_seq2, int *alnpos);

int buildcons(char **subrepseq, int *len_subrepseq, int num_subrep, char *conseq)
{
	int	i, j, k, l, m, n;
	int	len_conseq;
	int	*len_subseq;
	int	**alnpos;
	char	**subseq, *subcons;
	int	len_subcons;
	int	counts[5];

	l = 0;
	for(i = 0; i < num_subrep; i ++)	{
		if(l < len_subrepseq[i])	{
			l = len_subrepseq[i];
		}
	}
	l = 2 * l + 100;
	subcons = (char *) ckalloc(l * sizeof(char));
	len_subseq = (int *) ckalloc(num_subrep * sizeof(int));
	subseq = (char **) ckalloc(num_subrep * sizeof(char *));
	alnpos = (int **) ckalloc(num_subrep * sizeof(int *));
	for(i = 0; i < num_subrep - 1; i ++)	{
		alnpos[i] = (int *) ckalloc(l * sizeof(int));
		subseq[i] = (char *) ckalloc(l * sizeof(char));
		seq2ALIGN(subrepseq[num_subrep - 1], len_subrepseq[num_subrep - 1], subrepseq[i], len_subrepseq[i], alnpos[i]);
	}
	len_conseq = 0;
	for(i = 0; i < len_subrepseq[num_subrep - 1]; i ++)	{
		for(j = 0; j < 5; j ++)	{
			counts[j] = 0;
		}
		counts[subrepseq[num_subrep - 1][i]] ++;
		for(j = 0; j < num_subrep - 1; j ++)	{
			if(alnpos[j][i] == 0)	{
				counts[4] ++;
			} else	{
				counts[subrepseq[j][alnpos[j][i] - 1]] ++;
			}
		}
		m = counts[4];
		for(j = 0; j < 4; j ++)	{
			if(counts[j] > m)	{
				m = counts[j];
				conseq[len_conseq] = j;
			}
		}
		if(m != counts[4])	{
			len_conseq ++;
		}
		if(i == len_subrepseq[num_subrep - 1] - 1)	continue;

		m = 1;
		n = 0;
		for(j = 0; j < num_subrep - 1; j ++)	{
			if(alnpos[j][i] > 0 && alnpos[j][i + 1] - alnpos[j][i] > 1)	{
				for(k = alnpos[j][i]; k < alnpos[j][i + 1] - 1; k ++)	{
					subseq[n][len_subseq[n] ++] = subrepseq[j][k];
				}
				n ++;
			} else	{
				m ++;
			}
		}
		if(m < num_subrep / 2 && n >= 2)	{
			len_subcons = buildcons(subseq, len_subseq, n, subcons);
			for(j = 0; j < len_subcons; j ++)	{
				conseq[len_conseq ++] = subcons[j];
			}
		}
		for(j = 0; j < num_subrep; j ++)	{
			len_subseq[j] = 0;
		}
	}
	for(i = 0; i < num_subrep - 1; i ++)	{
		free((void *) alnpos[i]);
		free((void *) subseq[i]);
	}
	free((void **) alnpos);
	free((void **) subseq);
	free((void *) len_subseq);
	free((void *) subcons);
	return(len_conseq);
}

void seq2ALIGN(char *seq1, int len_seq1, char *seq2, int len_seq2, int *alnpos)
{
	int	i, j, k, l, m, n;
	int	*sapp;
	int	windows;
	int	gscore, ALIGN0();

	sapp = (int *) ckalloc((len_seq1 + len_seq2 + 1) * sizeof(int));
	windows = abs(len_seq1 - len_seq2) + 1000;
	gscore = ALIGN0(&seq1[-1], &seq2[-1], len_seq1, len_seq2, -windows, windows, W, 16, 4, sapp, len_seq1 + len_seq2 + 1,
			 len_seq1 + len_seq2 + 1);
	i = j = k = 0;
	l = 0;
	while(i < len_seq1 && j < len_seq2)	{
		if(sapp[k] == 0)	{
			alnpos[i] = j + 1;
			i ++;
			j ++;
			l ++;
		} else if(sapp[k] > 0)	{
			j += sapp[k];
		} else	{
			i -= sapp[k];
		}
		k ++;
	}
	free((void *) sapp);
}
