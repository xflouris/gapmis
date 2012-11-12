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

#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <stdio.h>
#include <getopt.h>
#include "types.h"

#define MAX_SIZE 1024
#define BUFFER_SIZE 128
#define NUC_SCORING_MATRIX_SIZE 15		
#define PRO_SCORING_MATRIX_SIZE 24		
#define ERR 24					//error number returned if char_to_index returns an invalid index

struct TSeq * read_fasta_file ( const char * szReadsFile );

int decode_switches ( int argc, char ** argv, struct TSwitch * );

unsigned int dp_algorithm ( double ** G, unsigned int ** H, char * t, unsigned int n, char * p, unsigned int m, unsigned int scoring_matrix, unsigned int MAXgap );

int nuc_delta ( char a, char b );

int pro_delta ( char a, char b );

unsigned int nuc_char_to_index ( char a );

unsigned int pro_char_to_index ( char a );

unsigned int i_limits ( unsigned int n, unsigned int m, unsigned int * up, unsigned int * down, unsigned int MAXgap );

unsigned int j_limits ( unsigned int i, unsigned int m, unsigned int * left, unsigned int * right, unsigned int MAXgap );

unsigned int opt_solution ( double ** G, unsigned int n, unsigned int m, unsigned int MAXgap, double gap_open_penalty, double gap_extend_penalty, double * MAXscore, unsigned int * MINgap, unsigned int * where, unsigned int * start );

double total_scoring ( unsigned int gap, double current_score, double gap_open_penalty, double gap_extend_penalty );

unsigned int backtracing ( unsigned int ** H, unsigned int m, unsigned int n, unsigned int start, unsigned int where, unsigned int * gap_pos );

unsigned int swap_txt_pat ( struct TSeq ** seqa, unsigned int * n, struct TSeq ** seqb, unsigned int * m );

#endif
