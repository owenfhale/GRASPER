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
#include <color.h>
#include <param.h>
#include <extfunc.h>

char	inpfile[100], outfile[100], seqfile[100];

void initenv(int argc, char **argv);

main(int argc, char **argv)
{
	int	i, j, k, l, m, n;
	int	dist[20];
	int	reads;
	int	num_vertex, num_edge;
	int	*len_seq, num_seq, num_remain;
	int	**num_pa;
	char	**src_seq, **src_name;
	char	temp[100];
	int	num_position;
	int	*pos1, *pos2, *dir1, *dir2;
	int	pos_rev1, pos_rev2, read_rev1, read_rev2;
	EDGE	**edge, *edge1, *edge2, *bal_edge1, *bal_edge2;
	PATH	*path;
	int	num_path;
	char	c;
	NODES	**vertex, *begin, *node, *node_rev, *node_next, **start_node;
	LIST	**list;
	READINTERVAL	*readinterval;
	POSITION	*position;
	FILE	*fp, *fp1;


	readpar();
	random1(&idum);
	initenv(argc, argv);

/*	Input the length of the reads (required) */

	len_seq = (int *) ckalloc(20 * sizeof(int));
	src_seq = (char **) ckalloc(20 * sizeof(char *));
	src_name = (char **) ckalloc(10 * sizeof(char *));
	for(i = 0; i < 10; i ++)	{
		src_name[i] = (char *) ckalloc(100 * sizeof(char));
	}

/*	for each component, build-repeat graph	*/

	fp = ckopen(seqfile, "r");
	num_seq = readseq1by1(src_seq, src_name, len_seq, fp);
	fclose(fp);
	printf("num_seq %d\n", num_seq);
	for(i = 0; i < num_seq; i ++)	{
		printf("Genome length: %d %d\n", i, len_seq[i]);
	}

/*	Make reverse complements of input sequences rev(i) --> i + num_seq	*/

	for(i = 0; i < num_seq; i ++)	{
		len_seq[i + num_seq] = len_seq[i];
		src_seq[num_seq + i] = (char *) ckalloc((len_seq[i] + 1) * sizeof(char));
		for(j = 0; j < len_seq[i]; j ++)	{
			src_seq[num_seq + i][j] = rev(src_seq[i][len_seq[i] - j - 1]);
		}
	}

/*	Input equivalent positions from read mapping --	*/

	fp = ckopen(inpfile, "r");
	fscanf(fp, "%d", &n);
	printf("Total equivalent positions: %d\n", n);
	n += 10;
	pos1 = (int *) ckalloc(n * sizeof(int));
	pos2 = (int *) ckalloc(n * sizeof(int));
	dir1 = (int *) ckalloc(n * sizeof(int));
	dir2 = (int *) ckalloc(n * sizeof(int));
	printf("Read equivalent positions...\n");
	num_position = readposition(pos1, pos2, dir1, dir2, len_seq, num_seq, fp);
	fclose(fp);
	printf("# equivalent positions input: %d\n", num_position);

/*	Initialize the nodes: each position in each read is assigned
	as a new node. An array of "list" is set up for each read	*/

	list = (LIST **) ckalloc(2 * num_seq * sizeof(LIST *));
	for(i = 0; i < 2 * num_seq; i ++)	{
		list[i] = (LIST *) ckalloc(len_seq[i] * sizeof(LIST));
	}
	printf("intitialize nodes...\n");
	initialize(list, len_seq, num_seq);
	printf("done.\n");
	n = countnode(list, len_seq, 2 * num_seq);
	printf("# of nodes before merge: %d\n", n);

/*	Glue together two nodes if their corresponding positions are defined
	as equivalent in a pairwise alignment		*/

/*
	for(i = 0; i < num_position; i ++)	{
		printf("read %d %d pos %d %d\n", dir1[i], dir2[i], pos1[i], pos2[i]);
	}
	getchar();
*/

	printf("Merge...\n");
	n = 0;
	for(i = 0; i < num_position; i ++)	{
		if(list[dir1[i]][pos1[i]].node == list[dir2[i]][pos2[i]].node)		{
			continue;
		}

// check if two nodes share two positions from the same sequence and within a distance
		c = check_width(list[dir1[i]][pos1[i]].node, list[dir2[i]][pos2[i]].node);
		if(c)	continue;

		node = chk_merge_node(list, dir1[i], dir2[i], pos1[i], pos2[i]);
		read_rev1 = reverse_read(dir1[i], num_seq);
		read_rev2 = reverse_read(dir2[i], num_seq);
		pos_rev1 = len_seq[dir1[i]] - 1 - pos1[i];
		pos_rev2 = len_seq[dir2[i]] - 1 - pos2[i];
		node_rev = chk_merge_node(list, read_rev1, read_rev2, pos_rev1, pos_rev2);
		node -> bal_node = node_rev;
		if(node != node_rev)	{
			node_rev -> bal_node = node;
		}
		n ++;
	}
	printf("done. %d pairs merged.\n", n);
	free((void *) pos1);
	free((void *) pos2);
	free((void *) dir1);
	free((void *) dir2);

/*      Compute the width of each node  */

        for(i = 0; i < 2 * num_seq; i ++)       {
                for(j = 0; j < len_seq[i]; j ++)        {
                        if(!list[i][j].node -> visit)   {
                                list[i][j].node -> num_path = countthickness(list[i][j].node);
                                list[i][j].node -> visit = 1;
                        }
                }
        }
	cleannode(list, len_seq, 2 * num_seq);
	n = countnode(list, len_seq, 2 * num_seq);
	printf("# of nodes after merge: %d\n", n);

/*	Add edges to the graph		*/
	edge = (EDGE **) ckalloc(n * sizeof(EDGE *));
	num_edge = graph(num_seq, len_seq, list, edge);
	printf("# edges: %d\n", num_edge);
	start_node = (NODES **) ckalloc(2 * num_seq * sizeof(NODES *));
	for(i = 0; i < 2 * num_seq; i ++)	{
		start_node[i] = list[i][0].node;
		free((void *) list[i]);
	}
	free((void **) list);

	vertex = (NODES **) ckalloc(2 * num_edge * sizeof(NODES *));
	num_vertex = count_vertex(edge, num_edge, vertex);
	free((void **) edge);

/*	Assign the complementary edges of each edge	*/
	for(i = 0; i < num_vertex; i ++)	{
		for(j = 0; j < vertex[i] -> num_nextedge; j ++)	{
			edge1 = vertex[i] -> nextedge[j];
			edge1 -> bal_edge = find_bal_edge(edge1, len_seq, num_seq, i);
		}
	}

/*	Remove bulges in the graph	*/
	printf("shave..\n");
	num_vertex = shave_graph(vertex, num_vertex);
	printf("done.\n");

/*	remove short edges	*/
/*
	printf("Remove shortedges...\n");
	num_vertex = rem_short_edge(vertex, num_vertex, len_seq);
	printf("done.\n%d vertices remained.\n", num_vertex);
*/

	num_pa = (int **) ckalloc(MAX_BRA * sizeof(int *));
	for(i = 0; i < MAX_BRA; i ++)	{
		num_pa[i] = (int *) ckalloc(MAX_BRA * sizeof(int));
	}
	num_edge = count_edge_simp(vertex, num_vertex, num_pa);
	printf("%d vertices %d edges (%d source %d sinks) remained.\n", num_vertex, num_edge,
		num_pa[0][1], num_pa[1][0]);

/*	Allocate the spaces for paths	*/
	printf("Allocating paths...\n");
	for(i = 0; i < num_vertex; i ++)	{
		vertex[i] -> num_path = 0;
	}

/*	Build sequence paths	*/
	printf("Define paths...\n");
	m = 0;
	for(i = 0; i < num_vertex; i ++)	{
		for(j = 0; j < vertex[i] -> num_nextedge; j ++)	{
			m += vertex[i] -> nextedge[j] -> multip;
		}
	}
	path = (PATH *) ckalloc(2 * num_seq * sizeof(PATH));
	for(i = 0; i < 2 * num_seq; i ++)	{
		path[i].edge = (EDGE **) ckalloc(m * sizeof(EDGE *));
	}
	printf("num_seq %d\n", num_seq);
	num_path = readpath(start_node, path, num_seq);
	printf("num_path %d\n", num_path);
	free((void **) start_node);
	num_edge = count_edge_simp(vertex, num_vertex, num_pa);
	m = l = 0;
	for(i = 0; i < num_vertex; i ++)	{
		for(j = 0; j < vertex[i] -> num_nextedge; j ++)	{
/*
if(vertex[i] -> nextedge[j] -> length == 1076)	{
	for(k = 0; k < vertex[i] -> nextedge[j] -> multip; k ++)	{
		printf("range %d %d\n", vertex[i] -> nextedge[j] -> readinterval[k].eq_read,
			vertex[i] -> nextedge[j] -> readinterval[k].begin);
	}
	getchar();
}
*/
			l += vertex[i] -> nextedge[j] -> length;
			if(vertex[i] -> nextedge[j] -> length > m)	{
				m = vertex[i] -> nextedge[j] -> length;
			}
		}
	}
	printf("%d vertics %d edges (%d source %d sinks) remained: total length %d (maximal %d).\n", num_vertex, num_edge,
	 	num_pa[0][1], num_pa[1][0], l, m);

/*	Make consensus of edges	*/
	initial_edge(vertex, num_vertex, src_seq, num_seq);
	printf("edge initialed\n");

/*	Output sequence path	*/

	n = 0;
	for(i = 0; i < num_vertex; i ++)	{
		vertex[i] -> visit = i;
		for(j = 0; j < vertex[i] -> num_nextedge; j ++)	{
			vertex[i] -> nextedge[j] -> start_cover = n;
			n ++;
		}
	}
	for(m = 0; m < num_seq; m ++)	{
		printf("SEQ%d ", m + 1);
		for(i = 0; i < path[m].len_path; i ++)	{
			printf("%d -- %d(%d,%d) --> ", path[m].edge[i] -> begin -> visit,
				path[m].edge[i] -> start_cover, path[m].edge[i] -> multip,
				path[m].edge[i] -> length);
			if(i % 5 == 4)	{
				printf("\n");
			}
		}
		printf("%d\n", path[m].edge[i - 1] -> end -> visit);
	}

/*	Output graph & contigs	*/
	sprintf(temp, "%s.edge", seqfile);
	fp = ckopen(temp, "w");
	sprintf(temp, "%s.graph", seqfile);
	fp1 = ckopen(temp, "w");
	write_graph(vertex, num_vertex, fp, fp1);
	fclose(fp);
	fclose(fp1);

/*	Output read intervals in each edge	*/
	sprintf(temp, "%s.intv", seqfile);
	fp = ckopen(temp, "w");
	write_interval(vertex, num_vertex, fp);
	fclose(fp);

/*	Output graphviz format graph	*/

	sprintf(temp, "%s", outfile);
	fp = ckopen(temp, "w");
	output_graph(vertex, num_vertex, fp);
	fclose(fp);

	for(i = 0; i < MAX_BRA; i ++)	{
		free((void *) num_pa[i]);
	}
	free((void **) num_pa);
	for(i = 0; i < 2 * num_seq; i ++)	{
		free((void **) path[i].edge);
	}
	free((void *) path);
	free_graph(vertex, num_vertex);
	for(i = 0; i < 2 * num_seq; i ++)	{
		free((void *) src_seq[i]);
	}
	for(i = 0; i < num_seq; i ++)	{
		free((void *) src_name[i]);
	}
	free((void **) src_seq);
	free((void **) src_name);
	free((void *) len_seq);
}

void initenv(int argc, char **argv)
{
	int copt;
	int inpseq, outseq;	
	extern char *optarg;

	inpseq = outseq = 0;

	while ((copt=getopt(argc,argv,"i:o:s:l:")) != EOF)	{
		switch(copt) {
			case 'i':
			  inpseq = 1;
			  sscanf(optarg,"%s", inpfile);
			  continue;
			case 's':
			  sscanf(optarg,"%s", seqfile);
			  continue;
			case 'o':
			  outseq = 1;
			  sscanf(optarg,"%s", outfile);
			  continue;
			case 'w':
			  sscanf(optarg, "%d", &MIN_WHIRL_SIZE);
			  continue;
			case 'l':
			  sscanf(optarg, "%d", &SHORTLEG);
			  continue;
			default:
			  printf("repeat_sin_read -i InpFile -s SeqFile -o outfile [-l shortleg -w MIN_WHIRL_SIZE]\n");
			  printf("-i InpFile: The input file name of alignments\n");
			  printf("-s SeqFile: The input file name of reads\n");
			  printf("-o Outfile: The output file name of contigs\n");
			  printf("-l shortleg (optional): default 50\n");
			  exit(-1);
		}
		optind--;
	}

	if(inpseq == 0 || outseq == 0)	{
		printf("repeat_sin_read -i InpFile -s SeqFile -o outfile [-l shortleg -w MIN_WHIRL_SIZE]\n");
		printf("-i InpFile: The input file name of alignments\n");
		printf("-s SeqFile: The input file name of reads\n");
		printf("-o Outfile: The output file name of contigs\n");
		printf("-l shortleg (optional): default 50\n");
		exit(-1);
	}
}
