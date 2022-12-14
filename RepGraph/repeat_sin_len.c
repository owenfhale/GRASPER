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

#define MAX_NUM 100

char	inpfile[100], outfile[100], lenfile[100], seqfile[100];

void initenv(int argc, char **argv);

main(int argc, char **argv)
{
	int	i, j, k, l, m, n;
	int	dist[20];
	int	reads;
	int	num_vertex, num_class, num_edge;
	int	*len_seq, num_seq, num_remain;
	int	**num_pa;
	char	**src_seq, **src_name;
	char	temp[100];
	ALIGN	**eq_class, *align;
	EDGE	**edge, *edge1, *edge2, *bal_edge1, *bal_edge2;
	PATH	*path;
	int	num_path;
	NODES	**vertex, *begin, *node, *node_next, **start_node;
	LIST	**list;
	READINTERVAL	*readinterval;
	POSITION	*position;
	FILE	*fp, *fp1;

	readpar();
	random1(&idum);
	initenv(argc, argv);
	printf("%d %d %d\n", sizeof(POSITION), sizeof(NODES), sizeof(LIST));

/*	Input the length of the genome (required) */

	len_seq = (int *) ckalloc(2 * MAX_NUM * sizeof(int));
	src_name = alloc_name(MAX_NUM, 100);
	fp = ckopen(lenfile, "r");
	num_seq = readlen(fp, len_seq, src_name);
	fclose(fp);

	src_seq = (char **) ckalloc(2 * num_seq * sizeof(char *));
	l = 0;
	printf("Genome length: ");
	for(i = 0; i < num_seq; i ++)	{
		l += len_seq[i];
		printf("%d ", len_seq[i]);
	}
	printf("\n");
	printf("Total length: %d\n", l);

/*	Make reverse complements of input sequences rev(i) --> i + num_seq	*/

	for(i = 0; i < num_seq; i ++)	{
		len_seq[i + num_seq] = len_seq[i];
		src_seq[i] = (char *) ckalloc(len_seq[i] * sizeof(char));
		src_seq[i + num_seq] = (char *) ckalloc(len_seq[i] * sizeof(char));
		for(j = 0; j < len_seq[i]; j ++)	{
			src_seq[num_seq + i][j] = rev(src_seq[i][len_seq[i] - j - 1]);
		}
	}

/*	Input equivalent readintervales between reads --
	see the format of the equivalent readinterval files	*/

	printf("Read equivalent readintervales...\n");
	eq_class = (ALIGN **) ckalloc(2 * num_seq * sizeof(ALIGN *));
	fp = ckopen(inpfile, "r");
	num_class = readclass(eq_class, num_seq, fp);
	fclose(fp);
	printf("# equivalent readintervales input: %d\n", num_class);

/*
	for(i = 0; i < 2 * num_seq; i ++)	{
		align = eq_class[i];
		while(align)	{
			printf("See: \n");
			output_align(align, src_name, src_seq, len_seq, num_seq);
			getchar();
			align = align -> next;
		}
	}
*/

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

	printf("Merge...\n");
	merge(num_seq, len_seq, eq_class, num_class, list);
	printf("done.\n");
	for(i = 0; i < num_seq; i ++)	{
		while(eq_class[i])	{
			eq_class[i] = free_align(eq_class[i]);
		}
	}
	free((void **) eq_class);

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
	start_node = (NODES **) ckalloc(num_seq * sizeof(NODES *));
	for(i = 0; i < num_seq; i ++)	{
		if(len_seq[i] > 0)	{
			start_node[i] = list[i][0].node;
		} else	{
			start_node[i] = (NODES *) NULL;
		}
	}
	for(i = 0; i < 2 * num_seq; i ++)	{
		free((void *) list[i]);
	}
	free((void **) list);

	vertex = (NODES **) ckalloc(2 * num_edge * sizeof(NODES *));
	num_vertex = count_vertex(edge, num_edge, vertex);
	free((void **) edge);

	num_pa = (int **) ckalloc(MAX_BRA * sizeof(int *));
	for(i = 0; i < MAX_BRA; i ++)	{
		num_pa[i] = (int *) ckalloc(MAX_BRA * sizeof(int));
	}
	num_edge = count_edge_simp(vertex, num_vertex, num_pa);
	printf("%d vertices %d edges (%d source %d sinks) remained.\n", num_vertex, num_edge,
		num_pa[0][1], num_pa[1][0]);

/*	Assign the complementary edges of each edge	*/
	for(i = 0; i < num_vertex; i ++)	{
		for(j = 0; j < vertex[i] -> num_nextedge; j ++)	{
			edge1 = vertex[i] -> nextedge[j];
			edge1 -> bal_edge = find_bal_edge(edge1, len_seq, num_seq, i);
		}
	}

/*	Remove bulges in the graph	*/
	printf("Shave...\n");
	num_vertex = shave_graph(vertex, num_vertex);
	printf("done.\n");

/*      Remove cycles shorter than some threshold in the graph  */
/*
        printf("Shaving graph...\n");
        num_vertex = rem_cycle(vertex, num_vertex);
        printf("done.\n%d vertices remained.\n", num_vertex);
*/

/*	remove short edges	*/
/*
	printf("Remove shortedges...\n");
	num_vertex = rem_short_edge(vertex, num_vertex, len_seq);
	printf("done.\n%d vertices remained.\n", num_vertex);
	fflush(stdout);
*/

	num_edge = count_edge_simp(vertex, num_vertex, num_pa);
	printf("%d vertices %d edges (%d source %d sinks) remained.\n", num_vertex, num_edge,
		num_pa[0][1], num_pa[1][0]);
	fflush(stdout);

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
	num_path = readpath(start_node, path, num_seq);
	free((void **) start_node);
	num_edge = count_edge_simp(vertex, num_vertex, num_pa);
	m = l = 0;
	for(i = 0; i < num_vertex; i ++)	{
		for(j = 0; j < vertex[i] -> num_nextedge; j ++)	{
			l += vertex[i] -> nextedge[j] -> length;
			if(vertex[i] -> nextedge[j] -> length > m)	{
				m = vertex[i] -> nextedge[j] -> length;
			}
		}
	}
	printf("%d vertics %d edges (%d source %d sinks) remained: total length %d (maximal %d).\n", num_vertex, num_edge,
	 	num_pa[0][1], num_pa[1][0], l, m);
	fflush(stdout);

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
		printf("len_path %d\n", path[m].len_path);
		printf("Sequence%d: ", m + 1);
		for(i = 0; i < path[m].len_path; i ++)	{
			printf("%d -- %d(%d,%d) --> ", path[m].edge[i] -> begin -> visit,
				path[m].edge[i] -> start_cover, path[m].edge[i] -> multip,
				path[m].edge[i] -> length);
			if(i % 5 == 4)	{
				printf("\n");
			}
		}
		if(path[m].len_path > 0)	{
			printf("%d\n", path[m].edge[i - 1] -> end -> visit);
		} else	{
			printf("\n");
		}
		fflush(stdout);
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
		if(path[i].len_path > 0)	{
			free((void **) path[i].edge);
		}
	}
	free((void *) path);
	free_graph(vertex, num_vertex);
	for(i = 0; i < 2 * num_seq; i ++)	{
		free((void *) src_seq[i]);
	}
	free((void **) src_seq);
	free_name(src_name, MAX_NUM);
	free((void *) len_seq);
}

void initenv(int argc, char **argv)
{
	int copt;
	int inpseq, outseq;	
	extern char *optarg;

	inpseq = outseq = 0;

	while ((copt=getopt(argc,argv,"i:o:s:l:w:L:")) != EOF)	{
		switch(copt) {
			case 'i':
			  inpseq = 1;
			  sscanf(optarg,"%s", inpfile);
			  continue;
			case 's':
			  sscanf(optarg,"%s", seqfile);
			  continue;
			case 'L':
			  sscanf(optarg,"%s", lenfile);
			  continue;
			case 'o':
			  outseq = 1;
			  sscanf(optarg,"%s", outfile);
			  continue;
			case 'w':
			  sscanf(optarg,"%d", &MIN_WHIRL_SIZE);
			  continue;
			case 'l':
			  sscanf(optarg, "%d", &SHORTLEG);
			  continue;
			default:
			  printf("repeat_sin_len -i InpFile -L LenFile -s SeqFile -o outfile [-l shortleg -w MIN_WHIRL_SIZE]\n");
			  printf("-i InpFile: The input file name of alignments\n");
			  printf("-L InpFile: The input file name of chromosome length\n");
			  printf("-o Outfile: The output file name of contigs\n");
			  printf("-l shortleg (optional): default 50\n");
			  exit(-1);
		}
		optind--;
	}

	if(inpseq == 0 || outseq == 0)	{
		printf("repeat_sin -i InpFile -L LenFile -s SeqFile -o outfile [-l shortleg -w MIN_WHIRL_SIZE]\n");
		printf("-i InpFile: The input file name of alignments\n");
		printf("-L InpFile: The input file name of chromosome length\n");
		printf("-o Outfile: The output file name of contigs\n");
		printf("-l shortleg (optional): default 50\n");
		exit(-1);
	}
}
