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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "types.h"
#include "output.h"

/* Usage of the tool */
void usage ( void )
 {
   fprintf ( stdout, "  ____             __  __ _\n" );
   fprintf ( stdout, " / ___| __ _ _ __ |  \\/  (_)___\n" );
   fprintf ( stdout, "| |  _ / _` | '_ \\| |\\/| | / __|\n" );
   fprintf ( stdout, "| |_| | (_| | |_) | |  | | \\__ \\\n" );
   fprintf ( stdout, " \\____|\\__,_| .__/|_|  |_|_|___/\n" );
   fprintf ( stdout, "            |_|\n\n" );

   fprintf ( stdout, "Usage: gapmis <options>\n" );
   fprintf ( stdout, "Standard (Mandatory):\n" );
   fprintf ( stdout, "  -a, --sequence-a          <str>     Sequence A filename.\n" );
   fprintf ( stdout, "  -b, --sequence-b          <str>     Sequence B filename.\n" );
   fprintf ( stdout, "Optional:\n" );
   fprintf ( stdout, "  -g, --gap-open-penalty    <float>   The gap open penalty is the score taken\n"
                     "                                      away when a gap is created. The best\n"
                     "                                      value depends on the choice of comparison\n"
                     "                                      matrix.   The   default   value   assumes\n"
                     "                                      you  are  using  the  EBLOSUM62  matrix\n"
                     "                                      protein  sequences,  and  the  EDNAFULL\n"
                     "                                      matrix for nucleotide sequences. Floating\n"
                     "                                      point number from 1.0 to 100.0. (default:\n"
                     "                                      10.0)\n" );


   fprintf ( stdout, "  -e, --gap-extend-penalty  <float>   The gap extension penalty is added to\n"
                     "                                      the standard gap penalty for each base or\n"
                     "                                      residue in the gap. This is how long gaps\n"
                     "                                      are penalized. Floating point number from\n"
                     "                                      0.0  to  10.0.  (default:  0.5)\n" );

   fprintf ( stdout, "  -o, --output-file         <str>     Output   alignement   filename   (default:\n"
                     "                                      gapmis.out)\n" );

   fprintf ( stdout, "  -d, --data-file           <str>     This  is  the  scoring  matrix  used  when\n"
                     "                                      comparing  sequences.  It  can  be  either\n"
                     "                                      `EBLOSUM62'  (for  protein  sequences)  or\n" 
                     "                                      `EDNAFULL'   (for  nucleotide  sequences).\n"  
                     "                                      (default: EDNAFULL)\n" );
   fprintf ( stdout, "  -m, --max-gap             <int>     Limit the maximum gap size to this value\n"
                     "                                      (default: length of the longest sequence\n"
                     "                                      minus  1)\n\n" );

 }

void print_header ( FILE * out, const char * filename, unsigned int matrix, double gap_penalty, double ext_penalty )
 {
   time_t               t;

   time ( &t );

   fprintf ( out, "####################################\n" );
   fprintf ( out, "# Program: GapMis\n" );
   fprintf ( out, "# Rundate: %s", ctime ( &t ) );
   fprintf ( out, "# Report file: %s\n", filename );
   fprintf ( out, "# Matrix: %s\n", ( matrix ? "BLOSUM62" : "EDNAFULL" ) );
   fprintf ( out, "# Gap penalty: %.3f\n", - gap_penalty );
   fprintf ( out, "# Extend penalty: %.3f\n", - ext_penalty );
   fprintf ( out, "####################################\n\n" );
 }

/* Creates the output file with the alignment */
unsigned int results ( const char * filename,    /* output filename */
                       struct TSeq * t,          /* text */          
                       unsigned int n,           /* length of t */
                       struct TSeq * p,          /* pattern */
                       unsigned int m,           /* length of p */
                       double max_score, 
                       unsigned int min_gap, 
                       unsigned int where, 
                       unsigned int gap_pos, 
                       unsigned int swap, 
                       unsigned int matrix,
                       double gap_penalty,
                       double ext_penalty )
 {

   FILE          * output;
   unsigned int    min_mis = 0;        //the number of mismatches in the alignment
   char          * seq_gap;            //the sequence with the inserted gap 
   char          * mark_mis;           //a string with the mismatches marked as '|' (and the matches as ' ')
   
   
   /* dynamic memory allocation for seq_gap */
   if ( ! ( seq_gap = ( char * ) calloc ( n + 1 + min_gap + 1, sizeof ( char ) ) ) )
    {
      fprintf ( stderr, "Error: seq_gap could not be allocated!!!\n" );
      return ( 0 );
    } 
   
   /* dynamic memory allocation for mark_mis */
   if ( ! ( mark_mis = ( char* ) calloc ( n + min_gap + 1, sizeof( char ) ) ) )
    {
      fprintf ( stderr, "Error: mark_mis could not be allocated!!!\n" );
      return ( 0 );
    } 
   
   if ( ! ( output = fopen ( filename, "w" ) ) )
    {
      fprintf ( stderr, "Error: cannot open file %s!!!\n", filename );
      return ( 0 );
    }
   
   print_header ( output, filename, matrix, gap_penalty, ext_penalty );
   
   if ( where == 1 )                //gap is in the text
    {
      print_alignment ( t -> data, n, p -> data, m, min_gap, gap_pos, seq_gap, mark_mis, &min_mis );
      free ( t -> data );
      t -> data = seq_gap;
      if ( ! swap )
       {       
         wrap ( t, p, mark_mis, LINE_LNG, output ); 
       }
      else
       {        
         wrap ( p, t, mark_mis, LINE_LNG, output ); 
       }
    }
   else                            //gap is in the pattern
    {
      print_alignment ( p -> data, m, t -> data, n, min_gap, gap_pos, seq_gap, mark_mis, &min_mis );
      free ( p -> data );
      p -> data = seq_gap;
      if ( ! swap )
       {        
         wrap ( t, p, mark_mis, LINE_LNG, output ); 
       }
      else        
       {        
         wrap ( p, t, mark_mis, LINE_LNG, output ); 
       }
    }
   
   fprintf ( output, "\n" );
   fprintf ( output, "Alignment score: %lf\n", max_score );
   fprintf ( output, "Number of mismatches: %d\n", min_mis );
   fprintf ( output, "Length of gap: %d\n", min_gap );
   
   if( min_gap > 0 )
    {
      fprintf ( output, "The gap was inserted in %.13s after position %d\n", ( where == 1 ) ? t -> header : p -> header, gap_pos );
    }
   
   if ( fclose ( output ) ) 
           fprintf ( stderr, "Error: cannot close file %s!!!\n", filename );
   
   free ( mark_mis );
   
   return ( 1 );	
 }

/*
Creates seq_gap and mark_mis, and computes min_mis
*/
unsigned int print_alignment ( const char * seqa, unsigned int seqa_len, const char * seqb, unsigned int seqb_len, unsigned int min_gap, unsigned int gap_pos, char * seq_gap, char * mark_mis, unsigned int* min_mis )
 {
   unsigned int i, j;

   if ( min_gap > 0 )
    {

      for ( i = 0; i < gap_pos; ++ i )
       {
         seq_gap[i] = seqa[i];
         if ( seqa[i] != seqb[i] )	
          {
            ++ ( *min_mis );
            mark_mis[i] = '.';
          }
         else				
           mark_mis[i] = '|';
       }

      for ( j = 0; j < min_gap; ++ j )
       {
         seq_gap[ j + i ] = '-'; 
         mark_mis[ j + i ] = ' ';
       }

      for ( ; i < seqb_len - min_gap && i < seqa_len ; ++ i )
       {
         seq_gap[j + i] = seqa[i];
         if ( seqa[i] != seqb[i + min_gap] )	
          {
            ++ ( *min_mis );
            mark_mis[ j + i ] = '.';
          }
         else
           mark_mis[j + i] = '|';
       }
      
      for ( ; i < seqa_len; ++ i )
       {
         seq_gap[j + i] = seqa[i];
         mark_mis[ j + i ] = '|';
       }
    }
   else
    {
      for ( i = 0; i < seqa_len; ++ i )
       {
         seq_gap[i] = seqa[i];
         if ( seqa[i] != seqb[i] )
          {
            ++ ( *min_mis );
            mark_mis[i] = '.';
          }
         else			
           mark_mis[i] = '|';
       }
    }	
   
   return ( 1 );
 }

void print_line ( const char * s, int start, int stop, int * nr_gaps, int end, int diff, FILE * output, const char * header )
 {
   int                  k;

   if ( start == stop ) return;

   if ( diff )
    {
      fprintf ( output, "%25s", "" );
    }
   else
    {
     if ( header )
       fprintf ( output, "%-13.13s %10d ", header, start + 1 - *nr_gaps );
     else
       fprintf ( output, "%-13.13s %10d ", "", start + 1 - *nr_gaps );
    }

   for ( ; start < stop; ++ start )
    {
      fputc ( s[start], output );
      if ( s[start] == '-' && ! diff ) ++ ( *nr_gaps );
    }

   if ( stop != end )
    {
      for ( k = stop; k < end; ++ k )
       {
         fputc ( ' ', output );
       }
    }
   if ( ! diff )  fprintf ( output, " %-10d", start - *nr_gaps );
   fprintf ( output, "\n" );
 }

/* Wrap two sequences s1 and s2 including the differences (diff) so that the
   line width is at most len
*/
void wrap ( struct TSeq * s1, struct TSeq * s2, const char * diff, int len, FILE * output )
 {
   int                  m, n, i, j;
   int                  nr_gaps_a;
   int                  nr_gaps_b;
   int                  nr_lines;

   if ( ! len ) return;

   n = strlen ( s1 -> data );
   m = strlen ( s2 -> data );

   if ( ! n && ! m ) return;

   i         = 0;
   j         = 0;
   nr_gaps_a = 0;
   nr_gaps_b = 0;

   //nr_lines = m / len;
   nr_lines = ( n > m ? m : n ) / len;
   for ( i = 0; i < nr_lines; ++ i )
    {
      /* Shortest sequence */
      print_line ( s1 -> data, i * len, ( i + 1 ) * len, &nr_gaps_a, ( i + 1 ) * len, 0, output, s1 -> header );

      /* Difference */
      print_line ( diff, i * len, ( i + 1 ) * len, NULL, ( i + 1 ) * len, 1, output, NULL );

      /* Longest sequence */
      print_line ( s2 -> data, i * len, ( i + 1 ) * len, &nr_gaps_b, ( i + 1 ) * len, 0, output, s2 -> header );
      fprintf ( output, "\n" );
    }

   /* Last line of first sequence and difference */
   j = i * len;
   if ( j < m || j < n ) 
    {
      print_line ( s1 -> data, i * len, min ( m, n ), &nr_gaps_a, ( i + 1 ) * len, 0, output, s1 -> header );
      print_line ( diff, i * len, ( m < n ) ? m : n, NULL, ( i + 1 ) * len, 1, output, NULL );
      print_line ( s2 -> data, i * len, min ( n, m), &nr_gaps_b, ( i + 1 ) * len, 0, output, s2 -> header );
    }
   
 }
