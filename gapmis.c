/**
    GapMis: a tool for pairwise sequence alignment with a single gap.
    Copyright (C) 2011 Solon P. Pissis, Tomas Flouri

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
**/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "functions.h"
#include "output.h"

int main ( int argc, char ** argv)
 {
   struct TSwitch  sw;
   char   * out_file;

   unsigned int MINgap;		//to be computed
   double MAXscore;
   
   unsigned int MAXgap;		//input arguments		
   double gap_open_pen;
   double gap_extend_pen;
   unsigned int scoring_matrix;
   
   unsigned int start;		//where to start backtracing
   unsigned int gap_pos;	//position of the gap
   unsigned int where;		//where is the gap: text or pattern
   
   struct TSeq * t;		//text t
   unsigned int n; 		//length of text
   struct TSeq * p;		//pattern p
   unsigned int m; 		//length of pattern
   
   double **       G; 		//dynamic programming matrix
   unsigned int ** H; 		//backtracing matrix
   	
   unsigned int swap;           //swap the text and the pattern in case m < n        
   unsigned int i;

   /* checks the arguments */
   i = decode_switches ( argc, argv, &sw );

   if ( i < 5 || ! sw . seq_a || ! sw . seq_b ) 
    {
      usage ();
      return ( 1 );
    }
   else 
    {
      gap_open_pen   = - sw . gap_open_pen;	//the penalties should have a negative value
      gap_extend_pen = - sw . gap_extend_pen;
      out_file       =   sw . out_file;

      if ( ! strcmp ( "EDNAFULL", sw . matrix ) )       scoring_matrix = 0;
      else if ( ! strcmp ( "EBLOSUM62", sw . matrix ) ) scoring_matrix = 1;
      else
       {
         fprintf ( stderr, "Error: scoring matrix argument should be `EDNAFULL' for nucleotide sequences or `EBLOSUM62' for protein sequences!!!\n" );
         return ( 1 );
       }
    }
   
   /* reads the text */
   t = read_fasta_file ( sw . seq_a );
   if ( ! t )
    {
      fprintf (stderr, "Error: cannot read file %s!!!\n", sw . seq_a );
      return ( 1 );
    }
   if ( ! t -> header ) t -> header = strdup ( "Seq A" );
   
   /* reads the pattern */
   p = read_fasta_file ( sw . seq_b );
   if ( ! p )
    {
      fprintf( stderr, "Error: cannot read file %s!!!\n", sw . seq_b );
      return ( 1 );
    }
   if ( ! p -> header ) p -> header = strdup ( "Seq B" );
   
   /* calculate text's and pattern's length */
   n = strlen ( t -> data );
   m = strlen ( p -> data );
   
   /* checks the lengths of text and pattern and swaps if needed */
   if ( m > n )
    {
      swap_txt_pat ( &t, &n, &p, &m );
      swap = 1;
    }
   else
      swap = 0;
   
   /* checks the max gap length allowed: MAXgap < n */
   MAXgap =  ( sw . max_gap <= -1 ) ?  n - 1 : sw . max_gap;
   
   if( MAXgap >= n )
    {
      fprintf ( stderr, "Error: the max gap length should be less than the length of the text!!!\n" );
      return ( 1 );
    }
   
   /* 2d dynamic memory allocation for matrices G and H*/
   if ( ! ( G = ( double ** ) malloc ( ( n + 1 ) * sizeof ( double * ) ) ) )
    {
      fprintf ( stderr, "Error: DP matrix could not be allocated!!!\n" );
      return ( 1 );
    } 
   
   if ( ! ( G[0] = ( double * ) calloc ( ( n + 1 ) * ( m + 1 ), sizeof ( double ) ) ) )
    {
      fprintf ( stderr, "Error: DP matrix could not be allocated!!!\n" );
      return ( 1 );
    } 
   
   for ( i = 1; i < n + 1; ++ i )
     G[i] = ( void * ) G[0] + i * ( m + 1 ) * sizeof ( double );
   
   if ( ! ( H = ( unsigned int ** ) malloc ( ( n + 1 ) * sizeof ( unsigned int * ) ) ) )
    {
      fprintf ( stderr, "Error: DP matrix could not be allocated!!!\n" );
      return ( 1 );
    } 
   
   if( ! ( H[0] = ( unsigned int * ) calloc ( ( n + 1 ) * ( m + 1 ) , sizeof ( unsigned int ) ) ) )
    {
      fprintf( stderr, "Error: DP matrix could not be allocated!!!\n" );
      return ( 1 );
    }
   
   for ( i = 1 ; i < n + 1 ; ++ i )
     H[i] = ( void * ) H[0] + i * ( m + 1 ) * sizeof ( unsigned int );
   
   /* dynamic programming algorithm */
   if ( ! ( dp_algorithm( G, H, t -> data, n, p -> data, m, scoring_matrix, MAXgap ) ) )
    {
      fprintf ( stderr, "Error: dp_algorithm() failed!!!\n" );
      return ( 1 );	
    }
   
   /* computes the optimal alignment based on the matrix score and the gap function */
   opt_solution ( G, n, m, MAXgap, gap_open_pen, gap_extend_pen, &MAXscore, &MINgap, &where, &start );
     
   /* computes the position of the gap */
   if ( MINgap > 0 ) backtracing ( H, m, n, start, where, &gap_pos );
   else gap_pos = 0;
   
   /* outputs the results */
   if ( ! ( results( out_file, t, n, p, m, MAXscore, MINgap, where, gap_pos, swap, scoring_matrix, gap_open_pen, gap_extend_pen ) ) )
    {
      fprintf(stderr, "Error: results() failed!!!\n");
      return ( 1 );	
    }
   
   free ( G[0] );
   free ( H[0] );
   free ( G );
   free ( H );
   
   free ( t -> data );
   free ( t -> header );
   free ( p -> data );
   free ( p -> header );
   free ( t );
   free ( p );
   free ( sw . out_file );
   free ( sw . matrix );
   return ( 0 );
 }



