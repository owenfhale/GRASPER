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

EDGE *detach_bal(EDGE *edge1, EDGE *edge2, PATH *path, int num_path, int *f_edge);

EDGE *detach_bal(EDGE *edge1, EDGE *edge2, PATH *path, int num_path, int *f_edge)
{
	int	i, j, k, l, m, n, q, n1, n2, c, min_c, min_l, min_i, num_r, l1, l2;
	double	r, min_r, min_min_r;
	int	num_endpath, num_startpath, num_midforpath, num_midaftpath;
	int	num_startmatch, num_endmatch;
	int	*num_match1, *num_match2, *num_cont;
	int	*edgematch1, *edgematch2;
	PATH	*endpath, *startpath;
	PATH	*midforpath, *midaftpath;
	NODES	*begin, *end, *vertex;
	EDGE	*edge, *newedge;

	if(edge1 -> end != edge2 -> begin)	{
		printf("edge1 %d %d edge2 %d %d\n", edge1, edge1 -> end, edge2, edge2 -> begin);
		printf("Edge not connected\n");
		return;
	}

	f_edge[0] = f_edge[1] = 0;

	vertex = edge1 -> end;

/*	num_startmatch--number of branches consistent with edge1
	num_endmatch--number of branches consistent with edge2	*/

	num_startmatch = countstartmatch(edge1, vertex, path, num_path);
	num_endmatch = countendmatch(edge2, vertex, path, num_path);
	if(num_startmatch == 0 || num_endmatch == 0)	{
		printf("Path not found %d(%d,%d-%d) %d %d(%d,%d-%d) %d \n", edge1, edge1 -> length,
			edge1 -> begin, edge1 -> end,
			num_startmatch, edge2, edge2 -> length, edge2 -> begin, edge2 -> end, num_endmatch);
		printf("vertex %d %d %d num_path %d\n", vertex, vertex -> num_lastedge, vertex -> num_nextedge,
			vertex -> num_path);
		exit(-1);
	}

	endpath = (PATH *) ckalloc(vertex -> num_path * sizeof(PATH));
	startpath = (PATH *) ckalloc(vertex -> num_path * sizeof(PATH));
	num_endpath = collectendpaths(vertex, edge1, endpath, path, num_path); /* P->x */
	num_startpath = collectstartpaths(vertex, edge2, startpath, path, num_path); /* Py-> */

	midforpath = (PATH *) ckalloc(vertex -> num_path * sizeof(PATH));
	midaftpath = (PATH *) ckalloc(vertex -> num_path * sizeof(PATH));
	num_cont = (int *) ckalloc((num_startpath + num_endpath) * sizeof(int));
	for(j = 0; j < num_endpath; j ++)	{
		num_cont[j] = 0;
	}

/*	num_match1--number of branches consistend with P->x	*/

	num_match1 = (int *) ckalloc(num_endpath * sizeof(int));
	edgematch1 = (int *) ckalloc(num_endpath * sizeof(int));
	for(i = 0; i < vertex -> num_nextedge; i ++)	{
		edge = vertex -> nextedge[i];
/*	if x and y' are not connected, num_midforpath <- 0		*/
		num_midforpath = collect2forpaths(vertex, edge1, edge, midforpath, path, num_path);  /* Px,y' */
		if(num_midforpath == 0)	{
			continue;
		}
		for(j = 0; j < num_endpath; j ++)	{
/*	check consistency of P->x and Px,y'	*/
			c = chk_consist(&endpath[j], midforpath, num_midforpath, &q);
/*	contained <--> c = 0;	consistent <--> c-1 = # of consistent;	inconsistent <--> c = -1	*/
			if(c >= 0)	{
/*	Px,y' and P->x are consistent	*/
				if(c == 0 && (num_cont[j] == 0 || edge == edge2))	{
					num_cont[j] = 1;
					if(edge == edge2)	{
						edgematch1[j] = 1;
					} else	{
						edgematch1[j] = 0;
					}
					num_match1[j] = 1;
				} else if(c > 0 && num_cont[j] == 0)	{
					num_match1[j] ++;
					if(edge == edge2)	{
						edgematch1[j] = 1;
					}
				}
			}
		}
	}

/*	n1--number of P->x matching more than one Px,y'	*/
	n1 = 0;
	for(j = 0; j < num_endpath; j ++)	{
		if(endpath[j].len_path > 0 && (num_match1[j] == 0 || (num_match1[j] > 1 && edgematch1[j] == 1)))	{
			n1 ++;
		}
	}


	for(j = 0; j < num_startpath; j ++)	{
		num_cont[j] = 0;
	}

/*	num_match2--number of branches consistend with Py->	*/

	num_match2 = (int *) ckalloc(num_startpath * sizeof(int));
	edgematch2 = (int *) ckalloc(num_startpath * sizeof(int));
	for(i = 0; i < vertex -> num_lastedge; i ++)	{
		edge = vertex -> lastedge[i];
/*	if x' and y are not connected, num_midaftpath <- 0		*/
		num_midaftpath = collect2aftpaths(vertex, edge, edge2, midaftpath, path, num_path);  /* Px',y */
		if(num_midaftpath == 0)	{
			continue;
		}
		l = 1;
		for(j = 0; j < num_startpath; j ++)	{
/*	check consistency of Py-> and Px',y	*/
			c = chk_consist(&startpath[j], midaftpath, num_midaftpath, &q);
/*	contained <--> c = 0;	consistent <--> c-1 = # of consistent;	inconsistent <--> c = -1	*/
/*	x',y and y-> are consistent	*/
			if(c >= 0)	{
				if(c == 0 && (num_cont[j] == 0 || edge == edge1))	{
					num_cont[j] = 1;
					if(edge == edge1)	{
						edgematch2[j] = 1;
					} else	{
						edgematch2[j] = 0;
					}
					num_match2[j] = 1;
				} else if(c > 0 && num_cont[j] == 0)	{
					num_match2[j] ++;
					if(edge == edge1)	{
						edgematch2[j] = 1;
					}
				}
			}
		}
	}

/*	n2--number of Py-> matching more than one Px',y	*/

	n2 = 0;
	for(j = 0; j < num_startpath; j ++)	{
		if(startpath[j].len_path > 0 && (num_match2[j] == 0 || (num_match2[j] > 1 && edgematch2[j] == 1)))	{
			n2 ++;
		}
	}

	free((void *) num_cont);

	newedge = NULL;
	if(n1 == 0 && n2 == 0)	{
		begin = edge1 -> begin;
		end = edge2 -> end;
		newedge = new_edge(vertex, begin, end, edge1, edge2);
		reducepath(path, num_path, vertex, edge1, edge2, newedge, edgematch1, edgematch2);
		if(num_startmatch == 1 && vertex -> num_lastedge > 1)	{
			n = searcherase(begin -> nextedge, edge1, begin -> num_nextedge);
			erasenext(begin, n);
			n = searcherase(vertex -> lastedge, edge1, vertex -> num_lastedge);
			eraselast(vertex, n);
			free((void *) edge1);
			f_edge[0] = 1;
		}
		if(num_endmatch == 1 && edge2 != edge1 && vertex -> num_nextedge > 1)	{
			n = searcherase(vertex -> nextedge, edge2, vertex -> num_nextedge);
			erasenext(vertex, n);
			n = searcherase(end -> lastedge, edge2, end -> num_lastedge);
			eraselast(end, n);
			free((void *) edge2);
			f_edge[1] = 1;
		}
	}

	free((void *) num_match1);
	free((void *) edgematch1);
	free((void *) num_match2);
	free((void *) edgematch2);
	for(i = 0; i < num_endpath; i ++)	{
		free((void **) endpath[i].edge);
	}
	free((void *) endpath);
	for(i = 0; i < num_startpath; i ++)	{
		free((void **) startpath[i].edge);
	}
	free((void *) startpath);
	for(i = 0; i < num_midforpath; i ++)	{
		free((void **) midforpath[i].edge);
	}
	for(i = 0; i < num_midaftpath; i ++)	{
		free((void **) midaftpath[i].edge);
	}
	free((void *) midforpath);
	free((void *) midaftpath);

	return(newedge);
}
