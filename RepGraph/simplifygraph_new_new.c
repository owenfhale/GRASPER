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
#define TANDEM_LEG 1000

char	inpfile[100], graphfile[100], edgefile[100], outfile[100], seqfile[100], intvfile[100], lenfile[100];
int	min_multip, min_length, min_multip2;
extern int thresh1;
extern char partmark;

void initenv(int argc, char **argv);

main(int argc, char **argv)
{
	int	i, j, k, l, m, n;
	int	dist[20];
	int	*len_seq;
	int	**num_pa;
	char	**src_name;
	BINDEX	*index;
	int	num_vertex, num_edge, num_seq, num_imp, num_vertex_new, num_path;
	char	temp[100];
	EDGE	**edge, *edge1, *edge2, *bal_edge1, *bal_edge2, **impedges;
	NODES	**vertex, **vertex_new, **start_node;
	READINTERVAL	*readinterval;
	PATH	*path;
	POSITION	*position;
	FILE	*fp, *fp1;

	readpar();
	random1(&idum);
	initenv(argc, argv);

/*	Input the length of the genome (required) */

	len_seq = (int *) ckalloc(2 * MAX_NUM * sizeof(int));
	src_name = alloc_name(MAX_NUM, 100);
	fp = ckopen(lenfile, "r");
	num_seq = readlen(fp, len_seq, src_name);
	fclose(fp);

/*	input graph and edge files	*/
	read_graph_file(edgefile, graphfile,
		  &num_vertex, &vertex,
		  &num_edge, &edge);
	printf("num_vertex %d\n", num_vertex);

/*	input interval files	*/
	fp = ckopen(intvfile, "r");
	read_interval(vertex, num_vertex, fp);
	fclose(fp);
	printf("num_vertex %d\n", num_vertex);

	num_pa = (int **) ckalloc(MAX_BRA * sizeof(int *));
	for(i = 0; i < MAX_BRA; i ++)	{
		num_pa[i] = (int *) ckalloc(MAX_BRA * sizeof(int));
	}
	num_edge = count_edge_simp(vertex, num_vertex, num_pa);
	printf("%d vertices %d edges (%d source %d sinks) in the beginning\n", num_vertex, num_edge,
		num_pa[0][1], num_pa[1][0]);

/*	Find important edges, i.e. multiplicity >= min_multip && Length >= min_length	*/

	num_vertex_new = 0;
	impedges = (EDGE **) ckalloc(num_edge * sizeof(EDGE *));
	vertex_new = (NODES **) ckalloc(num_vertex * sizeof(NODES *));
	num_imp = 0;
	for(i = 0; i < num_vertex; i ++)	{
		for(j = 0; j < vertex[i] -> num_nextedge; j ++)	{
			edge1 = vertex[i] -> nextedge[j];
			if(edge1 -> visit == 0 && edge1 -> multip >= min_multip)	{
				num_vertex_new = create_impedges(impedges, edge1, vertex_new, num_vertex_new, vertex, &num_imp);
			}
		}
	}
	free_graph(vertex, num_vertex);
	printf("num_vertex_new %d\n", num_vertex_new);
	num_edge = count_edge_simp(vertex_new, num_vertex_new, num_pa);
	printf("Subgraph of important edges only -- ");
	printf("%d vertices %d edges (%d source %d sinks) remained.\n", num_vertex_new, num_edge,
		num_pa[0][1], num_pa[1][0]);
	n = 0;
	for(i = 0; i < num_imp; i ++)	{
		n += impedges[i] -> multip;
	}
	printf("Important edges: %d. Total multiplicity: %d\n", num_imp, n);

/*	Output graph & contigs	*/
	sprintf(temp, "%s.new", edgefile);
	fp = ckopen(temp, "w");
	sprintf(temp, "%s.new", graphfile);
	fp1 = ckopen(temp, "w");
	write_graph(vertex_new, num_vertex_new, fp, fp1);
	fclose(fp);
	fclose(fp1);

/*	Output read intervals in each edge	*/
	sprintf(temp, "%s.new", intvfile);
	fp = ckopen(temp, "w");
	write_interval_single(vertex_new, num_vertex_new, fp);
	fclose(fp);

/*	Output graphviz format graph	*/

	sprintf(temp, "%s", outfile);
	fp = ckopen(temp, "w");
	output_graph(vertex_new, num_vertex_new, fp);
	fclose(fp);

	for(i = 0; i < MAX_BRA; i ++)	{
		free((void *) num_pa[i]);
	}
	free((void **) num_pa);
	free((void **) impedges);
	free_graph(vertex_new, num_vertex_new);
	for(i = 0; i < num_seq; i ++)	{
		free((void *) src_name[i]);
	}
	free((void **) src_name);
	free((void *) len_seq);
}

void initenv(int argc, char **argv)
{
	int copt;
	int inpseq, outseq;	
	extern char *optarg;

	inpseq = outseq = 0;
	min_multip = 1;
	min_multip2 = 2;
	min_length = 100;
	partmark = 0;

	while ((copt=getopt(argc,argv,"i:o:e:g:v:m:l:L:M:t:p")) != EOF)	{
		switch(copt) {
			case 'i':
			  inpseq = 1;
			  sscanf(optarg,"%s", inpfile);
			  sprintf(edgefile, "%s.edge", inpfile);
			  sprintf(intvfile, "%s.intv", inpfile);
			  sprintf(graphfile, "%s.graph", inpfile);
			  continue;
			case 'p':
			  partmark = 1;
			  continue;
			case 'L':
			  sscanf(optarg,"%s", lenfile);
			  continue;
			case 'e':
			  sscanf(optarg,"%s", edgefile);
			  continue;
			case 'g':
			  sscanf(optarg,"%s", graphfile);
			  continue;
			case 'v':
			  sscanf(optarg,"%s", intvfile);
			  continue;
			case 't':
			  sscanf(optarg,"%d", &thresh1);
			  continue;
			case 'o':
			  outseq = 1;
			  sscanf(optarg,"%s", outfile);
			  continue;
			case 'l':
			  sscanf(optarg, "%d", &min_length);
			  continue;
			case 'm':
			  sscanf(optarg, "%d", &min_multip);
			  continue;
			case 'M':
			  sscanf(optarg, "%d", &min_multip2);
			  continue;
			default:
			  printf("simplifygraph_new_len -i InpFile -o outfile -L lenfile [-e edgefile -g graphfile -v intvfile -M min_multiplicity -l minimal_edge_length -p]\n");
			  printf("-i InpFile: The input file name of sequences\n");
			  printf("-L lenfile: The input file name of chromosome length\n");
			  printf("-o Outfile: The output file name of graphs\n");
			  printf("-e edgefile (optional): default using inpfile.edge\n");
			  printf("-g graphfile (optional): default using inpfile.graph\n");
			  printf("-l min_multiplicity (optional): default 2\n");
			  printf("-M min_edge_length (optional): default 100\n");
			  exit(-1);
		}
		optind--;
	}

	if(inpseq == 0 || outseq == 0)	{
		  printf("simplifygraph_new_len -i InpFile -o outfile -L lenfile [-e edgefile -g graphfile -v intvfile -mM min_multiplicity -l minimal_edge_length -p]\n");
		  printf("-i InpFile: The input file name of sequences\n");
		  printf("-L lenfile: The input file name of chromosome length\n");
		  printf("-o Outfile: The output file name of graphs\n");
		  printf("-e edgefile (optional): default using inpfile.edge\n");
		  printf("-g graphfile (optional): default using inpfile.graph\n");
		  printf("-l min_edge_length (optional): default 100\n");
		  printf("-M min_multiplicity (optional): default 2\n");
		  exit(-1);
	}
}
