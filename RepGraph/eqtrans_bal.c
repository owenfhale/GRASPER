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

int nsuper, *numtangle, *maxmultip, *maxlength, *avelength, *mlength;

int eqtrans_bal(NODES **vertex, int num_vertex, PATH *path, int num_path);

int eqtrans_bal(NODES **vertex, int num_vertex, PATH *path, int num_path)
{
	int	i, j, k, l, n, m;
	int	tot_edge, tot_edge_old, tot_path, tot_path_old;
	int	**num_pa;
	int	votethresh;
	EDGE	*edge, **all_edge, *edge1, *edge2, *edge0;

	votethresh = LOW_COV;

	num_pa = (int **) ckalloc(MAX_BRA * sizeof(int *));
	for(i = 0; i < MAX_BRA; i ++)	{
		num_pa[i] = (int *) ckalloc(MAX_BRA * sizeof(int));
	}

	tot_edge = count_edge_simp(vertex, num_vertex, num_pa);
	printf("In the beginning, %d vertices, %d edges, %d sources/sinks.\n",
		num_vertex, tot_edge, num_pa[1][0] * 2);
	do	{
		tot_edge_old = tot_edge;
		tot_path_old = tot_path;
		num_vertex = vertex_run(vertex, num_vertex, path, num_path, votethresh);
		num_vertex = merge_graph_path(vertex, num_vertex, path, num_path, votethresh);
		tot_edge = count_edge_simp(vertex, num_vertex, num_pa);
		tot_path = count_path(path, num_path);
		printf("Transforming EULER graph: %d vertices, %d edges, %d sources and %d sinks.\n",
			num_vertex, tot_edge, num_pa[0][1], num_pa[1][0]);
	} while(tot_edge_old > tot_edge || tot_path_old > tot_path);

	num_vertex = splitbeg(vertex, num_vertex, path, num_path);
	num_vertex = merge_graph_path(vertex, num_vertex, path, num_path, votethresh);
	tot_edge = count_edge_simp(vertex, num_vertex, num_pa);
	printf("Finally %d vertices, %d edges, %d sources/sinks.\n", num_vertex, tot_edge, num_pa[0][1] * 2);

	for(i = 0; i < MAX_BRA; i ++)	{
		free((void *) num_pa[i]);
	}
	free((void **) num_pa);
	return(num_vertex);
}
