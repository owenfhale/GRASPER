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

int collect2forpaths(NODES *vertex, EDGE *edge1, EDGE *edge2, PATH *midforpath,
		     PATH *path, int num);
int collect2aftpaths(NODES *vertex, EDGE *edge1, EDGE *edge2, PATH *midaftpath,
		     PATH *path, int num);
int collectstartpaths(NODES *vertex, EDGE *edge, PATH *startpath, PATH *path, int num);
int collectendpaths(NODES *vertex, EDGE *edge, PATH *endpath, PATH *path, int num);

int collect2forpaths(NODES *vertex, EDGE *edge1, EDGE *edge2, PATH *midforpath,
		     PATH *path, int num)
{
	int	i, j, k, l, m, n, c, q;
	int	num_path;
	PATH	tmppath;

	tmppath.edge = (EDGE **) ckalloc(MAX_BRA * sizeof(EDGE *));

	num_path = 0;
	for(i = 0; i < vertex -> num_path; i ++)	{
		k = vertex -> path_index[i];
		j = vertex -> path_pos[i];
		if(j != 0 && j != path[k].len_path && path[k].edge[j - 1] == edge1 && path[k].edge[j] == edge2)	{
			tmppath.len_path = j - 1;
			for(m = 0; m < j - 1; m ++)	{
				tmppath.edge[m] = path[k].edge[j - 2 - m];
			}
			c = chk_consist(&tmppath, midforpath, num_path, &q);
			if(c > 0)	{
				midforpath[c - 1].len_path = tmppath.len_path;
				free((void **) midforpath[c - 1].edge);
				midforpath[c - 1].edge = (EDGE **) ckalloc((midforpath[c - 1].len_path + 1) * sizeof(EDGE *));
				for(m = 0; m < tmppath.len_path; m ++)	{
					midforpath[c - 1].edge[m] = tmppath.edge[m];
				}
			} else if(c < 0)	{
				midforpath[num_path].len_path = tmppath.len_path;
				midforpath[num_path].edge = (EDGE **) ckalloc((midforpath[num_path].len_path + 1) * sizeof(EDGE *));
				for(m = 0; m < tmppath.len_path; m ++)	{
					midforpath[num_path].edge[m] = tmppath.edge[m];
				}
				num_path ++;
			}
		}
	}

	free((void **) tmppath.edge);
	return(num_path);
}

int collect2aftpaths(NODES *vertex, EDGE *edge1, EDGE *edge2, PATH *midaftpath,
		     PATH *path, int num)
{
	int	i, j, k, l, m, n, c, q;
	int	num_path;
	PATH	tmppath;

	tmppath.edge = (EDGE **) ckalloc(MAX_BRA * sizeof(EDGE *));

	num_path = 0;
	for(i = 0; i < vertex -> num_path; i ++)	{
		k = vertex -> path_index[i];
		j = vertex -> path_pos[i];
		if(j != 0 && j != path[k].len_path && path[k].edge[j - 1] == edge1 && path[k].edge[j] == edge2)	{
			tmppath.len_path = path[k].len_path - j - 1;
			for(m = j + 1; m < path[k].len_path; m ++)	{
				tmppath.edge[m - j - 1] = path[k].edge[m];
			}
			c = chk_consist(&tmppath, midaftpath, num_path, &q);
			if(c > 0)	{
				midaftpath[c - 1].len_path = tmppath.len_path;
				free((void **) midaftpath[c - 1].edge);
				midaftpath[c - 1].edge = (EDGE **) ckalloc((midaftpath[c - 1].len_path + 1) * sizeof(EDGE *));
				for(m = 0; m < tmppath.len_path; m ++)	{
					midaftpath[c - 1].edge[m] = tmppath.edge[m];
				}
			} else if(c < 0)	{
				midaftpath[num_path].len_path = tmppath.len_path;
				midaftpath[num_path].edge = (EDGE **) ckalloc((midaftpath[num_path].len_path + 1) * sizeof(EDGE *));
				for(m = 0; m < tmppath.len_path; m ++)	{
					midaftpath[num_path].edge[m] = tmppath.edge[m];
				}
				num_path ++;
			}
		}
	}

	free((void **) tmppath.edge);
	return(num_path);
}

int collectendpaths(NODES *vertex, EDGE *edge, PATH *endpath, PATH *path, int num)
{
	int	i, j, k, l, m, n, c;
	int	num_path;

	num_path = 0;
	for(i = 0; i < vertex -> num_path; i ++)	{
		k = vertex -> path_index[i];
		j = vertex -> path_pos[i];
		if(path[k].len_path > 0 && j == path[k].len_path && path[k].edge[j - 1] == edge)	{
			endpath[num_path].len_path = j - 1;
			endpath[num_path].edge = (EDGE **) ckalloc((endpath[num_path].len_path + 1) * sizeof(EDGE *));
			for(m = 0; m < j - 1; m ++)	{
				endpath[num_path].edge[m] = path[k].edge[j - 2 - m];
			}
			num_path ++;
		}
	}

	return(num_path);
}

int collectstartpaths(NODES *vertex, EDGE *edge, PATH *startpath, PATH *path, int num)
{
	int	i, j, k, l, m, n, c;
	int	num_path;

	num_path = 0;
	for(i = 0; i < vertex -> num_path; i ++)	{
		k = vertex -> path_index[i];
		j = vertex -> path_pos[i];
		if(path[k].len_path > 0 && j == 0 && path[k].edge[j] == edge)	{
			startpath[num_path].len_path = path[k].len_path - 1;
			startpath[num_path].edge = (EDGE **) ckalloc((startpath[num_path].len_path + 1) * sizeof(EDGE *));
			for(m = 0; m < startpath[num_path].len_path; m ++)	{
				startpath[num_path].edge[m] = path[k].edge[m + 1];
			}
			num_path ++;
		}
	}

	return(num_path);
}
