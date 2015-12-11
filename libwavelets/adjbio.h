#include <string.h>
#include "bio_parameters.h"


extern void          adjbioD9_3d();
extern void          adjbioD9_skip_3d();
extern void          adjbioR9_3d();
extern void          adjbioR9_skip_3d();



extern void          adjbioD9_3d_char();
extern void          adjbioR9_3d_char();
extern void          adjbioD9_skip_3d_char();
extern void          adjbioR9_skip_3d_char();

extern void          adjbioD7_3d();
extern void          adjbioD7_skip_3d();
extern void          adjbioR7_3d();
extern void          adjbioR7_skip_3d();

extern void          adjbioD7_3d_char();
extern void          adjbioR7_3d_char();
extern void          adjbioD7_skip_3d_char();
extern void          adjbioR7_skip_3d_char();

extern void          adjbioD5_3d();
extern void          adjbioD5_skip_3d();
extern void          adjbioR5_3d();
extern void          adjbioR5_skip_3d();

extern void          adjbioD5_3d_char();
extern void          adjbioR5_3d_char();
extern void          adjbioD5_skip_3d_char();
extern void          adjbioR5_skip_3d_char();

extern void          adjbioD3_3d();
extern void          adjbioD3_skip_3d();
extern void          adjbioR3_3d();
extern void          adjbioR3_skip_3d();

extern void          adjbioD3_3d_char();
extern void          adjbioR3_3d_char();
extern void          adjbioD3_skip_3d_char();
extern void          adjbioR3_skip_3d_char();

extern void          adjbioD1_3d();
extern void          adjbioD1_skip_3d();
extern void          adjbioR1_3d();
extern void          adjbioR1_skip_3d();

extern void          adjbioD1_3d_char();
extern void          adjbioR1_3d_char();
extern void          adjbioD1_skip_3d_char();
extern void          adjbioR1_skip_3d_char();


/*new in declararion file: */
extern void   adjbioD1__3d();




extern void   adjbioR1__3d();




extern void bio_3d_premult(
					FLOAT * ,   /*Vector*/
					int ,       /* xlength*/
					int ,       /*ylength*/
					int ,       /*zlength*/
					FLOAT,     /*x1*/
					int        /* pairity*/
					);


extern void bio_3d_postmult(FLOAT * , /*Vector*/
					 int ,     /*xlength*/
					 int ,     /*ylength*/
					 int ,     /*zlength*/
					 FLOAT ,   /*x1*/
					 int /*pairity*/
					 );

void bioD1_3dskip(FLOAT*, /* vector*/
				  FLOAT* ,   /*LLL*/
				  int , /*xlength*/
				  int , /*ylength*/
				  int , /*zlength*/
				  FLOAT, /* X */
				  FLOAT /* norm */
				  );

void bioR1_3dskip(FLOAT* ,   /*LLL*/
				  FLOAT*, /* vector*/
				  int , /*xlength*/
				  int , /*ylength*/
				  int , /*zlength*/
				  FLOAT, /* X */
				  FLOAT /* norm */
				  );





extern void      getfilter();
extern void      getadjointfilter();



int wavelet_decompose3( FLOAT*, /*inspacevector*/
					   int, /*Xlength*/
					   int, /* Ylength*/
					   int, /*Zlength, use value 1 in 2dim */
					   char,/* Filterlength*/
					   char, /* Levels */
					   char , /*minZLevels, use  value 0  in 2dim */
					   char , /* MaxZLevels,   use value 0 in 2dim */
					   char , /*minXYLevels, use  value of Levels  in 2dim */
					   char , /* MaxXYLevels,   use value of Levels in 2dim */
					   char , /*Skip */ 
					   FLOAT*,/*covector*/			
                       int * , /*colength_ptr*/ 
						   char /* ifnotSilent*/
					   );


	 int wavelet_reconstruct3(FLOAT* , /*reccovector*/
							  int , /*colength*/
							  FLOAT*, /* outvector,  computation in space  */ 
							  int , /* Xlength*/
							  int , /*Ylength*/
							  int , /*Zlength, use value 1 in 2dim */
							  char, /* Filterlength*/
	       					  char ,/*Levels*/
							  char , /*minZLevels */
							  char , /* maxZLevels,   */
							  char , /*minXYLevels, use  value of Levels  in 2dim */
							  char , /* MaxXYLevels,   use value of Levels in 2dim */
							  char , /* Skip  */
							  char    /* ifnotSilent*/
					   );




















