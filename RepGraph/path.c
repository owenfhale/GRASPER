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

void resetpath(PATH *path, int num_path, NODES *vertex, EDGE *lastedge, EDGE *nextedge, EDGE *newedge, int begtag, int endtag);
void replacepath(PATH *path, int num_path, NODES *vertex, EDGE *edge, EDGE *newedge);
void searchpath(NODES *vertex, PATH *path, int num_path, int **match, int *match1, int *match2, int *num_match);
int count_path(PATH *path, int num_path);
int cutbegpath(PATH *path, int num_path, NODES *vertex, EDGE *edge, EDGE *nextedge);
int cutendpath(PATH *path, int num_path, NODES *vertex, EDGE *edge, EDGE *lastedge);
void remove_edge(PATH *path, int path_index, int path_pos);
int chk_path(int p, int *vp, int nvp);

int cutbegpath(PATH *path, int num_path, NODES *vertex, EDGE *edge, EDGE *nextedge)
{
	int	i, j, k, n, m, l, label, *pos;

	pos = (int *) ckalloc(num_path * sizeof(int));
	for(i = 0; i < vertex -> num_path; i ++)	{
		pos[i] = MAX_BRA;
	}

	m = 0;
	for(i = vertex -> num_path - 1; i >= 0; i --)	{
		k = vertex -> path_index[i];
		j = vertex -> path_pos[i];
		if(path[k].len_path <= 0)	{
			for(n = i; n < vertex -> num_path - 1; n ++)	{
				vertex -> path_index[n] = vertex -> path_index[n + 1];
				vertex -> path_pos[n] = vertex -> path_pos[n + 1];
			}
			vertex -> num_path --;
			continue;
		}
		if(j < path[k].len_path && path[k].edge[j] == nextedge)	{
			if(j > m)	{
				m = j;
			}
			for(l = 0; l < m; l ++)	{
				remove_edge(&path[k], k, 0);
			}
			pos[k] = j;
		}
	}
	for(i = vertex -> num_path - 1; i >= 0; i --)	{
		k = vertex -> path_index[i];
		j = vertex -> path_pos[i];
		if(j < pos[k])	{
			for(n = i; n < vertex -> num_path - 1; n ++)	{
				vertex -> path_index[n] = vertex -> path_index[n + 1];
				vertex -> path_pos[n] = vertex -> path_pos[n + 1];
			}
			vertex -> num_path --;
		}
	}

	free((void *) pos);

	return(m);
}

int cutendpath(PATH *path, int num_path, NODES *vertex, EDGE *edge, EDGE *lastedge)
{
	int	i, j, k, n, m, label, *pos;

	pos = (int *) ckalloc(num_path * sizeof(int));
	for(i = 0; i < vertex -> num_path; i ++)	{
		pos[i] = MAX_BRA;
	}

	m = 0;
	for(i = vertex -> num_path - 1; i >= 0; i --)	{
		k = vertex -> path_index[i];
		j = vertex -> path_pos[i];
		if(path[k].len_path <= 0)	{
			for(n = i; n < vertex -> num_path - 1; n ++)	{
				vertex -> path_index[n] = vertex -> path_index[n + 1];
				vertex -> path_pos[n] = vertex -> path_pos[n + 1];
			}
			vertex -> num_path --;
			continue;
		}
		if(j > 0 && path[k].edge[j - 1] == lastedge)	{
			if(path[k].len_path - j > m)	{
				m = path[k].len_path - j;
			}
			path[k].len_path = j;
			pos[k] = j;
		}
	}
	for(i = vertex -> num_path - 1; i >= 0; i --)	{
		k = vertex -> path_index[i];
		j = vertex -> path_pos[i];
		if(j > pos[k])	{
			for(n = i; n < vertex -> num_path - 1; n ++)	{
				vertex -> path_index[n] = vertex -> path_index[n + 1];
				vertex -> path_pos[n] = vertex -> path_pos[n + 1];
			}
			vertex -> num_path --;
		}
	}

	free((void *) pos);

	return(m);
}

void resetpath(PATH *path, int num_path, NODES *vertex, EDGE *lastedge, EDGE *nextedge, EDGE *newedge, int begtag, int endtag)
{
	int	i, j, k, l, n, m, c;
	int	label;
	int	*beg_path_visited, *end_path_visited, num_beg_path_visited, num_end_path_visited;
	NODES	*begin;

	beg_path_visited = (int *) ckalloc(vertex -> num_path * sizeof(int));
	end_path_visited = (int *) ckalloc(vertex -> num_path * sizeof(int));

	num_beg_path_visited = num_end_path_visited = 0;
	for(i = vertex -> num_path - 1; i >= 0; i --)	{
		k = vertex -> path_index[i];
		j = vertex -> path_pos[i];
		label = 0;
		if(path[k].len_path > 0)	{
			if(j > 0 && j < path[k].len_path && path[k].edge[j - 1] == lastedge && path[k].edge[j] == nextedge)	{
				path[k].edge[j - 1] = newedge;
				remove_edge(&path[k], k, j);
				label = 1;
			} else if(j == 0 && path[k].edge[j] == nextedge)	{
				if(begtag == 1)	{
					path[k].edge[j] = newedge;
					add_path(newedge -> begin, k, j);
					if(vertex != newedge -> begin)	{
						label = 1;
					}
				} else if(begtag == 0)	{
					c = chk_path(k, beg_path_visited, num_beg_path_visited);
					if(!c)	{
						beg_path_visited[num_beg_path_visited ++] = k;
						remove_edge(&path[k], k, j);
						label = 1;
					}
				}
			} else if(j == path[k].len_path && path[k].edge[j - 1] == lastedge)	{
				if(endtag == 1)	{
					path[k].edge[j - 1] = newedge;
					add_path(newedge -> end, k, j);
					if(vertex != newedge -> end)	{
						label = 1;
					}
				} else if(endtag == 0)	{
					c = chk_path(k, end_path_visited, num_end_path_visited);
					if(!c)	{
						end_path_visited[num_end_path_visited ++] = k;
						remove_edge(&path[k], k, j - 1);
						label = 1;
					}
				}
			}
		} else	{
			label = 1;
		}
		if(label == 1)	{
			for(n = i; n < vertex -> num_path - 1; n ++)	{
				vertex -> path_index[n] = vertex -> path_index[n + 1];
				vertex -> path_pos[n] = vertex -> path_pos[n + 1];
			}
			vertex -> num_path --;
		}
	}

	free((void *) beg_path_visited);
	free((void *) end_path_visited);
}

int chk_path(int p, int *vp, int nvp)
{
	int	i;

	for(i = 0; i < nvp; i ++)	{
		if(vp[i] == p)	{
			break;
		}
	}

	if(i == nvp)	{
		return(0);
	} else	{
		return(1);
	}
}

void replacepath(PATH *path, int num_path, NODES *vertex, EDGE *edge, EDGE *newedge)
{
	int	i, j, k, l, n, m;
	NODES	*begin;

	for(i = vertex -> num_path - 1; i >= 0; i --)	{
		k = vertex -> path_index[i];
		j = vertex -> path_pos[i];
		if(path[k].len_path <= 0)	{
			for(n = i; n < vertex -> num_path - 1; n ++)	{
				vertex -> path_index[n] = vertex -> path_index[n + 1];
				vertex -> path_pos[n] = vertex -> path_pos[n + 1];
			}
			vertex -> num_path --;
			continue;
		}
		if(j < path[k].len_path)	{
			if(path[k].edge[j] == edge)	{
				path[k].edge[j] = newedge;
			}
		}
	}
}

void remove_edge(PATH *path, int path_index, int path_pos)
{
	int	i, j, k, l, n;
	int	**pos;
	NODES	*vertex, **vertex_temp;

	vertex_temp = (NODES **) ckalloc(path -> len_path * sizeof(NODES *));
	pos = (int **) ckalloc((path -> len_path - 1) * sizeof(int *));
	for(i = 0; i < path -> len_path - 1; i ++)	{
		pos[i] = (int *) ckalloc(2 * max_seq * sizeof(int));
	}

	for(i = path_pos; i < path -> len_path - 1; i ++)	{
		path -> edge[i] = path -> edge[i + 1];
		vertex = vertex_temp[i - path_pos] = path -> edge[i] -> begin;
		for(j = 0; j < vertex -> num_path; j ++)	{
			if(vertex -> path_index[j] == path_index && vertex -> path_pos[j] == i + 1)	{
				pos[i - path_pos][j] = i;
			} else	{
				pos[i - path_pos][j] = -1;
			}
		}
	}

	vertex = path -> edge[path -> len_path - 1] -> end;
	for(j = vertex -> num_path - 1; j >= 0; j --)	{
		if(vertex -> path_index[j] == path_index && vertex -> path_pos[j] == path -> len_path)	{
			vertex -> path_pos[j] = path -> len_path - 1;
			break;
		}
	}

	path -> len_path -= 1;

	for(i = path_pos; i < path -> len_path; i ++)	{
		vertex = vertex_temp[i - path_pos];
		for(j = 0; j < vertex -> num_path; j ++)	{
			if(pos[i - path_pos][j] >= 0)	{
				vertex -> path_pos[j] = pos[i - path_pos][j];
			}
		}
	}

	for(i = 0; i < path -> len_path; i ++)	{
		free((void *) pos[i]);
	}
	free((void **) pos);
	free((void **) vertex_temp);
}


void searchpath(NODES *vertex, PATH *path, int num_path, int **match, int *match1, int *match2, int *num_match)
{
	int	i, j, k, n1, n2, m, label;
	NODES	*begin, *end;
	EDGE	*nextedge, *lastedge;

	for(i = 0; i < vertex -> num_path; i ++)	{
		k = vertex -> path_index[i];
		j = vertex -> path_pos[i];
		if(path[k].len_path <= 0)	{
			continue;
		}
		label = 0;
		for(n1 = 0; n1 < vertex -> num_lastedge; n1 ++)	{
			for(n2 = 0; n2 < vertex -> num_nextedge; n2 ++)	{
				lastedge = vertex -> lastedge[n1];
				nextedge = vertex -> nextedge[n2];
				if(j > 0 && j < path[k].len_path && path[k].edge[j - 1] == lastedge && path[k].edge[j] == nextedge)	{
					match[n1][n2] = 1;
					label = 1;
				} else if(j == path[k].len_path && path[k].edge[j - 1] == lastedge)	{
					match1[n1] = 1;
					num_match[0] ++;
					label = 1;
				} else if(j == 0 && path[k].edge[j] == nextedge)	{
					match2[n2] = 1;
					num_match[1] ++;
					label = 1;
				}
			}
		}
	}
}

int count_path(PATH *path, int num_path)
{
	int	i, n, m, q, j, k;
	NODES	*vertex;

	n = 0;
	m = 0;
	q = 0;
	for(i = 0; i < num_path; i ++)	{
		if(path[i].len_path > 1)	{
			q += path[i].len_path;
			n ++;
			if(path[i].len_path == 1)	{
				m ++;
			}
		} else if(path[i].len_path < 0)	{
			printf("Path length error!\n");
			exit(-1);
		}
	}

	return(q);
}
