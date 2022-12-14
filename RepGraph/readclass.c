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

int readclass(ALIGN **eq_class, int num_seq, FILE *fp);

int readclass(ALIGN **eq_class, int num_seq, FILE *fp)
{
	int	i, j, k, l, n;
	int	num_class;
	ALIGN	*class;

	l = 0;
	for(i = 0; i < 2 * num_seq; i ++)	{
		fread(&n, sizeof(int), 1, fp);
		l += n;
		for(j = 0; j < n; j ++)	{
			class = (ALIGN *) ckalloc(1 * sizeof(ALIGN));
			class -> reads[0] = i;
			fread(&(class -> reads[1]), sizeof(int), 1, fp);
			fread(&(class -> mis_match), sizeof(int), 1, fp);
			fread(&(class -> length), sizeof(int), 1, fp);
			for(k = 0; k < 2; k ++)	{
				class -> pos[k] = (int *) ckalloc(class -> length * sizeof(int));
			}
			fread(class -> pos[0], sizeof(int), class -> length, fp);
			fread(class -> pos[1], sizeof(int), class -> length, fp);
			class -> next = eq_class[i];
			class -> last = NULL;
			if(eq_class[i])	{
				eq_class[i] -> last = class;
			}
			eq_class[i] = class;
		}
	}
	return(l);
}
