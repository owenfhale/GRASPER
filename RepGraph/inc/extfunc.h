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

extern ALIGN *AlignReads(char **src_seq, char **score, int *len_seq, int num_seq, READINTERVAL *readinterval, int index);
extern int overalign(char *seq1, char *seq2, char *score1, char *score2, 
	      int len_seq1, int len_seq2, int offset, int *sapp, int *cutp, int *seql,
	      int index1, int index2, int num_seq);
extern double perc_qual(char *score, int s1, int s2);
extern int readreputer(ALIGN **align, char **src_seq, int *len_seq, FILE *fp, int MIN_LEG, double MIN_ID);
extern ALIGN *new_align(int *sapp, ALIGN *align, int r1, int r2, int pos1, int pos2, int len1, int len2);
extern double cal_identity(char *seq1, int len1, char *seq2, int len2, int *sapp);
extern void buildblock(int num_seq, int *len_seq, EQCLASS **eq_readinterval, LIST **blocklist, BLOCK *block);
extern char char2intgen(char c);
extern char char2int(char c);
extern void *ckalloc(int amount);
extern FILE *ckopen(char *name, char *mode);
extern int alignlength(ALIGN *align);
extern ALIGN *remove_trans(ALIGN *align1, ALIGN *align2, ALIGN *align3);
extern ALIGN *chkconn(int i1, int i2, ALIGN **eq_readinterval, int num_seq);
extern char chkedgelink(EDGE *edge);
extern char chksingle(NODES *node, EDGE *edge, EDGE *sedge, int offset, int cumoff);
extern int count_edge(NODES **vertex, int num_vertex, int **num_pa);
extern int count_edge_simp(NODES **vertex, int num_vertex, int **num_pa);
extern int forwardtanglelen(EDGE *edge, int *ave, int *multip);
extern int backtanglelen(EDGE *edge, int *ave, int *multip);
extern int count_tangle(NODES **vertex, int num_vertex, int disttangle[][7]);
extern int count_bal(NODES *vertex);
extern void getmaxedge(EDGE *edge, EDGE **maxedge);
extern int count_vertex(EDGE **edge, int num_edge, NODES **vertex);
extern int collect_vertex(NODES *v, NODES **vertex, int num_vertex);
extern EDGE *find_bal_edge(EDGE *edge, int *len_seq, int num_seq, int index);
extern char chk_readinterval(READINTERVAL *readinterval1, READINTERVAL *readinterval2, int n1, int n2, int *len_seq, int num_seq);
extern int singstatnode(NODES *node, int **nstat);
extern void statnode(LIST **list, int *len_seq, char **src_seq, int num_seq);
extern int countnode(LIST **list, int *len_seq, int num_seq);
extern int cleannode(LIST **list, int *len_seq, int num_seq);
extern char del_trans(EDGE *edge1, EDGE *edge2);
extern long trans_seq(char *seq, int len);
extern LINKPOS *ins_linkpos(LINKPOS *linkpos, int i1, int i2);
extern HASH *free_hash(HASH *hash);
extern HASH *ins_hash(HASH *hash, int i1, int i2, char **src_seq, int num_seq);
extern READINTERVAL *ins_readinterval(READINTERVAL *readinterval, int readindex, int offset);
extern int detectoverlap(char **src_seq, char **score, int *len_seq, int num_seq, READINTERVAL **readinterval);
extern int comp_word(char *w1, char *w2, int len);
extern READINTERVAL *new_readinterval(HASH *hash, int i1, int i2, char **src_seq, READINTERVAL *readinterval, int num_seq);
extern void erasedge(EDGE *edge);
extern void erasenext(NODES *vertex, int n);
extern void eraselast(NODES *vertex, int n);
extern int searcherase(EDGE **edge, EDGE *e, int num);
extern void errcorrt(LIST **list, int *len_seq, char **src_seq, char **score, char **src_name,
	      int num_seq, FILE *fp, FILE *fp1);
extern READINTERVAL *free_readinterval(READINTERVAL *readinterval);
extern int size_readinterval(READINTERVAL *readinterval);
extern void free_graph(NODES **vertex, int num_vertex);
extern int GALIGN(char *A, char *B, char *SA, char *SB, int M, int N, int low, int up, int W[][15],
	  int G, int H, int *S, int MW, int MX);
extern INDEX *insert_index(INDEX *index, int p);
extern INDEX *free_index(INDEX *index);
extern void initialize(LIST **list, int *len_seq, int num_seq);
extern EDGE *insert_edge(NODES *node, NODES *node_next, int read, int pos1, int pos2);
extern void  insert_readinterval(EDGE *edge, READINTERVAL readinterval_new);
extern void add_nextedge(NODES *node, EDGE *edge);
extern void add_lastedge(NODES *node, EDGE *edge);
extern int LOCAL_ALIGN(char *A, char *B, char *SA, char *SB, int M, int N, int low, int up, int W[][15],
		int G, int H, int *psi, int *psj, int *pei, int *pej, int MW);
extern int makedge(NODES *node, EDGE **edge, int num_edge, LIST **list);
extern EDGE *newedge(EDGE **midedge, int num_midedge, LIST **list);
extern int readintervalcompar(READINTERVAL *a, READINTERVAL *b);
extern int copyreadinterval(READINTERVAL *readinterval, READINTERVAL *readcov);
extern void sortreadinterval(READINTERVAL *readinterval, int multip);
extern void insert_readcov(READINTERVAL **readcov, POSITION *position, int offset);
extern void sort_nodepos(NODES *node);
extern int poscompar(POSPOINT *a, POSPOINT *b);
extern void ins_node(LIST **list, int read1, int read2, int pos1, int pos2, int endpos);
extern NODES *chk_merge_node(LIST **list, int read1, int read2, int pos1, int pos2);
extern void insert_position(NODES *node, int read, int pos);
extern void merge(int num_seq, int *len_seq, ALIGN **eq_readinterval, int num_readinterval, LIST **list);
extern NODES *combine_nodes(NODES *node1, NODES *node2);
extern NODES *free_nodes(NODES *node);
extern void update_link(NODES *node, int r1, int s1, int r2, int s2);
extern void output_graph(NODES **vertex, int num_vertex, FILE *fp);
extern int  ran_number(int n, int *idum);
extern double random1(int *idum);
extern int readreadinterval(ALIGN **eq_readinterval, int num_seq, FILE *fp);
extern int readoverlap(READINTERVAL **readinterval, FILE *fp);
extern void allocpath(PATH *path, int len);
extern void findpath(PATH *path, EDGE *edge, int *len_seq, int num_seq, char *label);
extern void backpath(PATH *path, NODES *vertex, int reads, int end, int *len_seq);
extern void forpath(PATH *path, NODES *vertex, int reads, int begin, int *len_seq);
extern int readqual(char **val, int *len, FILE *fp);
extern int readscore(char **val, int *len, FILE *fp);
extern int readseq1by1(char **src_seq, char **src_name, int *len_seq, FILE *fp);
extern int readseq1by1gen(char **src_seq, char **src_name, int *len_seq, FILE *fp);
extern READINTERVAL *insert_unitig(READINTERVAL *unitig, NODES *node, EDGE *edge, int *len_seq);
extern int search_unitig(NODES **nodeall, NODES *node, READINTERVAL **unitig, int nunitig, int num_seq, int *len_seq);
extern int merge_graph(NODES **vertex, int num_vertex);
extern EDGE *merge_vertex(NODES *vertex);
extern int shave_graph(NODES **vertex, int num_vertex);
void combine_readinterval(EDGE *edge1, EDGE *edge2, EDGE *newedge);
extern EDGE *new_edge(NODES *vertex, NODES *begin, NODES *end, EDGE *edge1, EDGE *edge2);
extern void combine_readinterval(EDGE *edge1, EDGE *edge2, EDGE *newedge);
extern int findposition(READINTERVAL *readinterval, int num, int readindex, int position);
extern void visit_comp(NODES **node, NODES *begin_node, int *comp, int num_comp, int num_seq);
extern int countwidth(NODES *node);
extern int countthickness(NODES *node);
extern void countallnode(LIST **list, int *len_seq, int num_seq);
extern ALIGN *free_align(ALIGN *align);
extern int size_align(ALIGN *align);
extern char chkedge(NODES *node1, NODES *node2);
extern EDGE *merge_vertex_path(NODES *vertex, PATH *path, int num_path);
extern int merge_graph_path(NODES **vertex, int num_vertex, PATH *path, int num_path, int small_cover);
extern int eqtrans_bal(NODES **vertex, int num_vertex, PATH *path, int num_path);
extern int rm_edge(NODES *vertex, PATH *path, int num_path);
extern void replace1edge(PATH *path, int num_path, EDGE *edge1, EDGE *edge2);
extern void reducepath(PATH *path, int num_path, NODES *vertex, EDGE *edge1, EDGE *edge2, EDGE *newedge, 
		int *edgematch1, int *edgematch2);
extern int numedge(EDGE *edge, PATH *path);
extern EDGE *detach_bal(EDGE *edge1, EDGE *edge2, PATH *path, int num_path, int *f_edge);
extern int searchlast(NODES *vertex, EDGE *edge, PATH *path, int num_path, int *match);
extern int searchnext(NODES *vertex, EDGE *edge, PATH *path, int num_path, int *match);
extern int countmatch(EDGE *edge1, EDGE *edge2, PATH *path, int num_path);
extern int countstartmatch(EDGE *edge, NODES *vertex, PATH *path, int num_path);
extern int countendmatch(EDGE *edge, NODES *vertex, PATH *path, int num_path);
extern int collect2forpaths(NODES *vertex, EDGE *edge1, EDGE *edge2, PATH *midforpath,
		     PATH *path, int num);
extern int collect2aftpaths(NODES *vertex, EDGE *edge1, EDGE *edge2, PATH *midaftpath,
		     PATH *path, int num);
extern int collectstartpaths(NODES *vertex, EDGE *edge, PATH *startpath, PATH *path, int num);
extern int collectendpaths(NODES *vertex, EDGE *edge, PATH *endpath, PATH *path, int num);
extern void resetpath(PATH *path, int num_path, NODES *vertex, EDGE *lastedge, EDGE *nextedge, EDGE *newedge, int begtag, int endtag);
extern void replacepath(PATH *path, int num_path, NODES *vertex, EDGE *edge, EDGE *newedge);
extern void searchpath(NODES *vertex, PATH *path, int num_path, int **match, int *match1, int *match2, int *num_match);
extern int count_path(PATH *path, int num_path);
extern int cutbegpath(PATH *path, int num_path, NODES *vertex, EDGE *edge, EDGE *nextedge);
extern int cutendpath(PATH *path, int num_path, NODES *vertex, EDGE *edge, EDGE *lastedge);
extern void remove_edge(PATH *path, int path_index, int path_pos);
extern int chk_path(int p, int *vp, int nvp);
extern int splitbeg(NODES **vertex, int num_vertex, PATH *path, int num_path);
extern NODES *new_vertex(NODES *vertex);
extern void move_path_next(EDGE *nextedge, NODES *vertex, NODES *vertex0, PATH *path);
extern void move_path_last(EDGE *lastedge, NODES *vertex, NODES *vertex0, PATH *path);
extern int chk_consist(PATH *startpath, PATH *midpath, int num_midpath, int *n);
extern void add_path(NODES *vertex, int path_index, int path_pos);
extern void rem_path(NODES *vertex, int path_index, int path_pos);
extern void set_path(NODES **vertex, int num_vertex, PATH *path, int num_path);
extern void merge_gap(char **src_seq, int num_seq, int *len_seq, ALIGN **eq_readinterval, int num_readinterval, LIST **list, int gap_k);
extern void errcorrt_pair(ALIGN **eq_readinterval, int **multip, int *len_seq, char **src_seq, char **score, char **src_name,
		   int num_seq, int LOW_COV, int intv, FILE *fp, FILE *fp1);
extern void trim_align(ALIGN *align, char **src_seq, int *len_seq);
extern int graph(int num_seq, int *len_seq, LIST **list, EDGE **edge);
extern void movereadinterval(EDGE *edge1, EDGE *edge);
extern int realmultip(EDGE *edge);
extern int readblast(ALIGN **align, FILE *fp, int MIN_LEG, double MIN_ID);
extern void count_multip(NODES **vertex, int num_vertex);
extern void output_contig(NODES **vertex, int num_vertex, char **src_name, char **src_seq,
		 int num_seq, FILE *fp, FILE *fp1, FILE *fp2, FILE *fp3);
extern void writeseq(FILE *fp, char *seq, char *name, int length);

extern void initial_edge(NODES **vertex, int num_vertex, char **src_seq, int num_seq);
extern void initedge(EDGE *edge, char **src_seq);
extern int readpath(NODES **start_node, PATH *path, int num_seq);
extern int singlepath(NODES *start_node, PATH *path, int reads, int begin);
extern int trans_loc(int loc, int *range, int num);
extern int rem_short_edge(NODES **vertex, int num_vertex, int *len);
extern int readposition(int *pos1, int *pos2, int *dir1, int *dir2, int *len_seq, int num_seq, FILE *fp);
extern void extedgeseq(EDGE *edge1, EDGE *edge2);
extern int readins(int *insertpos, int *insertlen, FILE *fp);
extern int readcase(ALIGN **align, char **src_seq, int *len_seq, int num_seq, FILE *fp, int MIN_LEG, double MIN_ID);
extern void extract_filename(char *item, char *filename);
extern void getaln(char *seq1, char *seq2, int len1, int len2, int *sapp, char *filename);
extern int add_int_edge(EDGE **impedges, BINDEX *index, int num, NODES **vertex, int num_vertex, int pos, int num_seq, int len_seq, NODES **start_node);
extern void buildedge(NODES *node1, NODES *node2, int index, int sp, int ep, int len_seq, int num_seq);
extern EDGE *single_edge(NODES *node1, NODES *node2, int index, int sp, int ep);
extern int create_impedges(EDGE **impedges, EDGE *edge, NODES **vertex_new, int num_vertex_new, NODES **vertex, int *num_imp);
extern NODES *initiate_node(NODES *node1, NODES *node);
extern EDGE *initiate_edge(EDGE *edge1, EDGE *edge);
extern int getstr(char *str, int *range, int *len_seq, char **chrname, int num_chro);
extern int getseq(char *str, char *seq, int len);
extern int readph(ALIGN **align, char **src_seq, int *len_seq, FILE *fp, int MIN_LEG, double MIN_ID);
extern double seq2sapp(char *seq1, char *seq2, int len1, int len2, int *sapp, int *len);
extern int readstring(ALIGN **align, FILE *fp, int MIN_LEG, double MIN_ID, int *len_seq, char **chrname, int num_chro);
extern int readlen(FILE *fp, int *len_seq, char **chrname);
extern char **alloc_name(int n1, int n2);
extern char **free_name(char **chrname, int n1);
extern double Kimura2D(double TransitionRate, double TranversionRate);
extern double Align2Dist(char *seq1, char *seq2, int len1, int len2);
extern int copyseq(char *seq1, char *seq2, int len_seq);
extern int reverse_seq(char *seq, int len);
extern double simulate_pair(int length1, int length2, int copynum1, int copynum2, int lengththresh,
	double *distchr, int *num_chrseg, int num_chr, int num_seg, int seglength);
extern double simulate_pair_all(int length1, int length2, int copynum1, int copynum2, int lengththresh,
	double *distchr, int *num_chrseg, int num_chr, int num_seg, int seglength, int simtimes);
