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

int total_color = 21;

double color_box[21][3] = {1, 0, 0, 0, 1, 0, 0, 0, 1,
			   1, 1, 0, 1, 0, 1, 0, 1, 1,
			   0.3, 1, 0.3, 0.3, 1, 0.6, 0.3, 1, 1,
			   0.6, 1, 0.3, 0.6, 1, 0.6, 0.6, 1, 1,
			   1, 1, 0.3, 1, 1, 0.6, 0, 1, 0.3,
			   0, 1, 0.6, 1, 0.5, 0, 0, 0.5, 1,
			   1, 0.5, 1, 0, 0.5, 0, 0.5, 0.5, 0.5};
