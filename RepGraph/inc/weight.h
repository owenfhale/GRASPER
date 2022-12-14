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

/*                g  a  c  t  n  r  y  w  s  m  k  h  b  v  d  */
int W[15][15] = { 4,-3,-3,-3, 1, 1,-3,-3, 1,-3, 1,-3, 1, 1, 1,   /* g */
		 -3, 4,-3,-3, 1, 1,-3, 1,-3, 1,-3, 1,-3, 1, 1,   /* a */
		 -3,-3, 4,-3, 1,-3, 1,-3, 1, 1,-3, 1, 1, 1,-3,   /* c */
		 -3,-3,-3, 4, 1,-3, 1, 1,-3,-3, 1, 1, 1,-3, 1,   /* t */
		  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,   /* n */
		  1, 1,-3,-3, 1, 1,-3, 1, 1, 1, 1, 1, 1, 1, 1,   /* r */
		 -3,-3, 1, 1, 1,-3, 1, 1, 1, 1, 1, 1, 1, 1, 1,   /* y */
		 -3, 1,-3, 1, 1, 1, 1, 1,-3, 1, 1, 1, 1, 1, 1,   /* w */
		  1,-3, 1,-3, 1, 1, 1,-3, 1, 1, 1, 1, 1, 1, 1,   /* s */
		 -3, 1, 1,-3, 1, 1, 1, 1, 1, 1,-3, 1, 1, 1, 1,   /* m */
		  1,-3,-3, 1, 1, 1, 1, 1, 1,-3, 1, 1, 1, 1, 1,   /* k */
		 -3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,   /* h */
		  1,-3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,   /* b */
		  1, 1, 1,-3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,   /* v */
		  1, 1,-3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};  /* d */
