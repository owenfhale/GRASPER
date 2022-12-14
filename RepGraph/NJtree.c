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

#define MAX_DIST 20000

void NeighborJoin(int num_seq, double *EvoDist, int *tree, double *dist);
void printNJtree(char *treeseq, int *tree, double *dist, int num_seq, char **segname);

void NeighborJoin(int num_seq, double *EvoDist, int *tree, double *dist)
{
	int	i, j, k, l, m, n, k1;
	int	mini, minj;
	int	num_step;
	char	*label;
	double	mins, d, *ut;

	if(num_seq <= 1)	{
		tree[0] = tree[1] = 0;
		dist[0] = dist[1] = EvoDist[0];
		return;
	}
	label = (char *) ckalloc(2 * num_seq * sizeof(char));
	ut = (double *) ckalloc(2 * num_seq * sizeof(double));
	num_step = 0;
	n = num_seq;
	while(num_step < n - 1)	{
		mini = minj = -1;
		for(i = 0; i < num_seq; i ++)	{
			if(label[i] == 1)	continue;
			if(mini >= 0)	minj = i;
			else		mini = i;
		}
		if(n - num_step == 2)	{
			k = numc(mini, minj);
			tree[2 * num_step] = mini;
			tree[2 * num_step + 1] = minj;
			dist[2 * num_step] = EvoDist[k];
			dist[2 * num_step + 1] = EvoDist[k];
			num_step ++;
			continue;
		}
		for(i = 0; i < num_seq; i ++)	{
			ut[i] = 0;
			if(label[i] == 1)	continue;
			for(j = 0; j < num_seq; j ++)	{
				if(label[j] == 1)	continue;
				if(i == j)	continue;
				k = numc(i, j);
				ut[i] += EvoDist[k];
			}
			ut[i] /= (n - num_step - 2);
		}
		mins = MAX_DIST;
		for(i = 0; i < num_seq; i ++)	{
			if(label[i] == 1)	continue;
			for(j = i + 1; j < num_seq; j ++)	{
				if(label[j] == 1)	continue;
				k = numc(i, j);
				d = EvoDist[k] - ut[i] - ut[j];
				if(d < mins)	{
					mins = d;
					mini = i;
					minj = j;
				}
			}
		}
		k = numc(mini, minj);
		tree[2 * num_step] = mini;
		tree[2 * num_step + 1] = minj;
		dist[2 * num_step] = 0.5 * EvoDist[k] + 0.5 * (ut[mini] - ut[minj]);
		dist[2 * num_step + 1] = 0.5 * EvoDist[k] + 0.5 * (ut[minj] - ut[mini]);
		label[mini] = label[minj] = 1;
		num_step ++;
		for(i = 0; i < num_seq; i ++)	{
			if(label[i] == 1)	continue;
			l = numc(mini, i);
			j = numc(minj, i);
			k = numc(mini, minj);
			k1 = numc(i, num_seq);
			EvoDist[k1] = (EvoDist[l] + EvoDist[j] - EvoDist[k]) / 2;
		}
		num_seq ++;
	}
	free((void *) ut);
	free((void *) label);
}

void printNJtree(char *treeseq, int *tree, double *dist, int num_seq, char **segname)
{
	int	i, j, k, l;
	char	**tmpstr, *tmpstr1;
	int	len;

	if(num_seq <= 1)	{
		printf("Can't build the tree: sequence less than 2.\n");
		exit(0);
	}

	tmpstr = (char **) ckalloc(num_seq * sizeof(char *));
	tmpstr1 = (char *) ckalloc(num_seq * 1000 * sizeof(char));
	len = 0;
	for(i = 0; i < num_seq - 1; i ++)	{
		tmpstr[i] = (char *) ckalloc((2 * len + 1000) * sizeof(char));
		if(tree[2 * i] < num_seq)	{
			sprintf(tmpstr[i], "(%s:%.4f, ", segname[tree[2 * i]], dist[2 * i]);
		} else	{
			sprintf(tmpstr[i], "(%s:%.4f, ", tmpstr[tree[2 * i] - num_seq], dist[2 * i]);
		}
		if(tree[2 * i + 1] < num_seq)	{
			sprintf(tmpstr1, "%s:%.4f)", segname[tree[2 * i + 1]], dist[2 * i + 1]);
		} else	{
			sprintf(tmpstr1, "%s:%.4f)", tmpstr[tree[2 * i + 1] - num_seq], dist[2 * i + 1]);
		}
		strcat(tmpstr[i], tmpstr1);
		l = strlen(tmpstr[i]);
		if(l > len)	{
			len = l;
		}
	}
	strcpy(treeseq, tmpstr[num_seq - 2]);
	for(i = 0; i < num_seq - 1; i ++)	{
		free((void *) tmpstr[i]);
	}
	free((void **) tmpstr);
	free((void *) tmpstr1);
}
