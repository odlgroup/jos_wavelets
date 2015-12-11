#include <string.h>
#include "bio_parameters.h"



extern void          bioD9_3d();
extern void          bioD9_skip_3d();
extern void          bioR9_3d();
extern void          bioR9_skip_3d();



extern void          bioD9_3d_char();
extern void          bioR9_3d_char();
extern void          bioD9_skip_3d_char();
extern void          bioR9_skip_3d_char();

extern void          bioD7_3d();
extern void          bioD7_skip_3d();
extern void          bioR7_3d();
extern void          bioR7_skip_3d();

extern void          bioD7_3d_char();
extern void          bioR7_3d_char();
extern void          bioD7_skip_3d_char();
extern void          bioR7_skip_3d_char();

extern void          bioD5_3d();
extern void          bioD5_skip_3d();
extern void          bioR5_3d();
extern void          bioR5_skip_3d();

extern void          bioD5_3d_char();
extern void          bioR5_3d_char();
extern void          bioD5_skip_3d_char();
extern void          bioR5_skip_3d_char();

extern void          bioD3_3d();
extern void          bioD3_skip_3d();
extern void          bioR3_3d();
extern void          bioR3_skip_3d();

extern void          bioD3_3d_char();
extern void          bioR3_3d_char();
extern void          bioD3_skip_3d_char();
extern void          bioR3_skip_3d_char();

extern void          bioD1_3d();
extern void          bioD1_skip_3d();
extern void          bioR1_3d();
extern void          bioR1_skip_3d();

extern void          bioD1_3d_char();
extern void          bioR1_3d_char();
extern void          bioD1_skip_3d_char();
extern void          bioR1_skip_3d_char();

/*new in declararion file: */
extern void   bioD1__3d();


extern void   bioR1__3d();





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

int wavelet_invadjoint3( FLOAT*, /*inspacevector*/
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


	 int wavelet_adjoint3(FLOAT* , /*reccovector*/
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

 

 




























void boundarymultiply(FLOAT *,//vector,
		      int, //xlength,
		      int, //ylength,
		      int,// zlength,
		      FLOAT //factor
		      );

void bpairitymultiply(FLOAT * ,//vector,
		      int, // xlength,
		      int, // ylength,
		      int, // zlength,
                      int, // pairity,
		      FLOAT // factor
		      );











