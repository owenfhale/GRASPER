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

#define pow1(x) ((x) * (x))
#define rev(x) ((2 + (x)) % 4)
#define reverse_read(x, y) ((x) >= (y) ? (x - y) : (x + y))
#define max(x,y) ((x) > (y) ? (x) : (y))
#define min(x,y) ((x) < (y) ? (x) : (y))
#define numc(i, j) ((i) < (j)? (j)*((j)-1)/2+(i): (i)*((i)-1)/2+(j))
#define numl(i, j) ((i) < (j)? (j)*((j)+1)/2+(i): (i)*((i)+1)/2+(j))

#define NUM_CHRO 24

typedef struct index {
	int index;
	struct index *next;
} INDEX;

typedef struct position {
	char readindex;
	int  position;
	struct position *next;
} POSITION;

typedef struct POSPOINT {
	POSITION *position;
} POSPOINT;

typedef struct nodes {
	POSITION *position;
	int	npos;
	struct nodes	*bal_node;
	int	visit;
	int	num_path;
	int	num_nextedge;
	int	num_lastedge;
	struct	edge **nextedge;
	struct	edge **lastedge;
} NODES;

typedef struct readinterval {
	int	begin;
	int	eq_read;
	int	length;
	int	offset;
	char	cov;
	struct readinterval *next;
} READINTERVAL;

typedef struct edge {
	NODES	*begin, *end;
	char	*seq;
	char	visit;
	int	start_cover;
	int	subrepeat;
	int	length;
	int	multip;
	struct edge	*bal_edge;
	READINTERVAL	*readinterval;
} EDGE;

typedef struct eqclass {
	READINTERVAL	*readinterval1;
	READINTERVAL	*readinterval2;
	struct readinterval *next;
} EQCLASS;

typedef struct align {
	int	reads[2];
	int	*pos[2];
	int	length;
	int	mis_match;
	struct align *next;
	struct align *last;
} ALIGN;

typedef struct linkpos {
	int readindex, position;
	struct linkpos *next;
} LINKPOS;

typedef struct hash {
	LINKPOS    *linkpos;
	struct hash *next;
} HASH;

typedef struct table {
	HASH    *prev;
} TABLE;

typedef struct block {
	struct block *bal_block;
	struct block *last;
	struct block *next;
	int	multiplicity;
	READINTERVAL	*readinterval;
} BLOCK;

typedef struct bintree {
	struct bintree *left, *right;
	int	readindex;
	int	minpos, maxpos;
} BINTREE;

typedef struct list {
	NODES	*node;
} LIST;

typedef struct bindex {
	int	index;
	int	index_mul;
	int	begin;
} BINDEX;

typedef struct segment {
	int	chro;
	long	length;
	long 	pos[2];
	long 	src_pos[2];
	long	eq_pos[2];
} SEGMENT;

typedef struct path {
	EDGE   **edge;
	int len_path;
	int begin_length;
	int end_length;
} PATH;

typedef struct line {
        int     bgnpoint[2];
        int     endpoint[2];
} LINE;

typedef struct reg {
	int	index;
	int	chro;
	int	range[2];
} REGION;
