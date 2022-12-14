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

int copyseq(char *seq1, char *seq2, int len_seq);
int reverse_seq(char *seq, int len);

int copyseq(char *seq1, char *seq2, int len_seq)
{
	int	i, len;

	len = 0;
	for(i = 0; i < len_seq; i ++)	{
		seq1[i] = seq2[i];
	}
	return(len_seq);
}

int reverse_seq(char *seq, int len)
{
	int	i;
	char	c;

	for(i = 0; i < len / 2; i ++)	{
		c = rev(seq[i]);
		seq[i] = rev(seq[len - 1 - i]);
		seq[len - 1 - i] = c;
	}
	if(len % 2 == 1)	{
		seq[len / 2] = rev(seq[len / 2]);
	}
}
