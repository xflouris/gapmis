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

#ifndef OUTPUT_H
#define OUTPUT_H

#include "types.h"

#define LINE_LNG 50

void usage ( void );

unsigned int results ( const char * filename, struct TSeq * t, unsigned int n, struct TSeq * p, unsigned int m, double MAXscore, unsigned int MINgap, unsigned int where, unsigned int gap_pos, unsigned int swap, unsigned int matrix, double gap_penalty, double ext_penalty );

unsigned int print_alignment( const char * seqa, unsigned int seqa_len, const char * seqb, unsigned int seqb_len, unsigned int MINgap, unsigned int gap_pos, char * seq_gap, char * mark_mis, unsigned int * MINmis );

void print_header ( FILE * out, const char * filename, unsigned int matrix, double gap_penalty, double ext_penalty );

void wrap ( struct TSeq * s1, struct TSeq * s2, const char * diff, int len, FILE * output );

void wrap2 ( struct TSeq * s1, struct TSeq * s2, const char * diff, int len, FILE * output );

void print_line ( const char * s, int start, int stop, int * nr_gaps, int end, int diff, FILE * output, const char * header );

#endif
