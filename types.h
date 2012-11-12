/**
    GapMis: a tool for pairwise sequence alignment with a single gap.
    Copyright (C) 2011 Solon P. Pissis, Tomas Flouri

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
**/

#ifndef TYPES_H
#define TYPES_H

#define max(a,b) ((a) > (b)) ? (a) : (b)
#define min(a,b) ((a) < (b)) ? (a) : (b)
#define SWAP(x,y) { ( x ) ^= ( y ); ( y ) ^= ( x ); ( x ) ^= ( y ); }

struct TSeq
 {
   char               * header;
   char               * data;
 };

struct TSwitch
 {
   char               * seq_a;
   char               * seq_b;
   char               * out_file;
   char               * matrix;
   double               gap_open_pen;
   double               gap_extend_pen;
   int                  max_gap;
 };

#endif
