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

NODES *chk_merge_node(LIST **list, int read1, int read2, int pos1, int pos2);
void insert_position(NODES *node, int read, int pos);
void merge(int num_seq, int *len_seq, ALIGN **eq_class, int num_class, LIST **list);
NODES *combine_nodes(NODES *node1, NODES *node2);
NODES *free_nodes(NODES *node);
void update_link(NODES *node, int r1, int s1, int r2, int s2);

void merge(int num_seq, int *len_seq, ALIGN **eq_class, int num_class, LIST **list)
{
	int 	i, j, k, l, m, n;
	char	c1, c2;
	ALIGN	*align;
	int	read1, read2, read_rev1, read_rev2, pos1, pos2, pos_rev1, pos_rev2;
	NODES	*node, *node_rev, *node_next, *node1, *node2;

	n = 0;
	for(i = 0; i < 2 * num_seq; i ++)	{
		align = eq_class[i];
		while(align)	{
			read1 = align -> reads[0];
			read2 = align -> reads[1];
			read_rev1 = reverse_read(read1, num_seq);
			read_rev2 = reverse_read(read2, num_seq);
			for(j = 0; j < align -> length - 1; j ++)	{
				pos1 = align -> pos[0][j];
				pos2 = align -> pos[1][j];
				pos_rev1 = len_seq[read1] - pos1 - 1;
				pos_rev2 = len_seq[read2] - pos2 - 1;
				while(pos1 < align -> pos[0][j + 1] && pos2 < align -> pos[1][j + 1])	{
					node = chk_merge_node(list, read1, read2, pos1, pos2);
					node_rev = chk_merge_node(list, read_rev1, read_rev2, pos_rev1, pos_rev2);
					node -> bal_node = node_rev;
					if(node != node_rev)	{
						node_rev -> bal_node = node;
					}
					pos1 ++;
					pos2 ++;
					pos_rev1 --;
					pos_rev2 --;
				}
			}
			n ++;
			align = align -> next;
		}
	}
	printf("# merged overlaps: %d\n", n);
}

NODES *chk_merge_node(LIST **list, int read1, int read2, int pos1, int pos2)
{
	int	i, j, k, l;
	char	c;
	POSITION *position;
	NODES	*node1, *node2, *node3, *new_node;

	if(pos1 >= 0 && pos2 >= 0)	{
		node1 = list[read1][pos1].node;
		node2 = list[read2][pos2].node;
		if(node1 == node2)	{
			update_link(node1, read1, pos1, read2, pos2);
			return(node1);
		} else	{
			new_node = combine_nodes(node1, node2);
			update_link(new_node, read1, pos1, read2, pos2);
			position = new_node -> position;
			while(position)	{
				read2 = position -> readindex;
				pos2 = position -> position;
/*	Link lists to new node	*/
				if(pos2 >= 0)	{
					list[read2][pos2].node = new_node;
				}
				position = position -> next;
			}
			free((void *) node1);
			free((void *) node2);
			return(new_node);
		}
	} else	{
		printf("Negative positions: %d %d %d %d\n", read1, pos1, read2, pos2);
		exit(-1);
	}
}

void insert_position(NODES *node, int read, int pos)
{
	POSITION *position, *position1;

	position1 = node -> position;
	while(position1)	{
		if(position1 -> readindex == read && position1 -> position == pos)	{
			return;
		}
		position1 = position1 -> next;
	}
	position = (POSITION *) ckalloc(1 * sizeof(POSITION));
	position -> readindex = read;
	position -> position = pos;
	position -> next = node -> position;
	node -> position = position;
	node -> npos ++;
}

void update_link(NODES *node, int r1, int s1, int r2, int s2)
{
	int	i, j, k, l, t1, t2;
	POSITION	*position;

/*
	t1 = t2 = -1;
	position = node -> position;
	while(position)	{
		if(position -> readindex == r1 && position -> position == s1)	{
			t1 = i;
		}
		if(position -> readindex == r2 && position -> position == s2)	{
			t2 = i;
		}
		position = position -> next;
	}
	if(t1 < 0 || t2 < 0)	{
		printf("node %d Not found position %d %d %d %d %d %d\n", node, r1, s1, r2, s2, t1, t2);
		exit(-1);
	}
*/
}

NODES *free_nodes(NODES *node)
{
	POSITION	*position;

	while(node -> position)	{
		position = node -> position -> next;
		free((void *) node -> position);
		node -> position = position;
	}
	if(node -> lastedge)	free((void **) node -> lastedge);
	if(node -> nextedge)	free((void **) node -> nextedge);
	free((void *) node);
	return((NODES *) NULL);
}

NODES *combine_nodes(NODES *node1, NODES *node2)
{
	int	i, k, j, l, n;
	POSITION *pos, *pos1, *pos2, *pos3;
	NODES	*new_node, *node01, *node02;

	new_node = (NODES *) ckalloc(1 * sizeof(NODES));
	new_node -> npos = node1 -> npos + node2 -> npos;
	if(node1 -> npos > node2 -> npos)	{
		node01 = node2;
		node02 = node1;
	} else	{
		node01 = node1;
		node02 = node2;
	}
	pos = new_node -> position = node01 -> position;
	while(pos -> next)	{
		pos = pos -> next;
	}
	pos -> next = node02 -> position;
	pos = new_node -> position;
	while(pos)	{
		if(pos -> position < 0)	{
			pos2 = pos;
			pos1 = pos -> next;
			while(pos1)	{
				if(pos1 -> readindex == pos -> readindex &&
				   pos1 -> position == pos -> position)	{
					pos2 -> next = pos1 -> next;
					free((void *) pos1);
					new_node -> npos --;
				} else	{
					pos2 = pos1;
				}
				pos1 = pos2 -> next;
			}
		}
		pos = pos -> next;
	}
	return(new_node);
}
