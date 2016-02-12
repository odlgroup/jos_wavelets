// Copyright 2015-2016 Jan-Olov Stromberg <jostromb@kth.se>
//
// jos_wavelets is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// jos_wavelets is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with jos_wavelets.  If not, see <http://www.gnu.org/licenses/>.


#include <string.h>
#include "bio_parameters.h"

extern void          bioD9_2dzz();
extern void          bioD9_skip_2d();
extern void          bioR9_2dzz();
extern void          bioR9_skip_2d();



extern void          bioD9_2dzz_char();
extern void          bioR9_2dzz_char();
extern void          bioD9_skip_2d_char();
extern void          bioR9_skip_2d_char();

extern void          bioD7_2dzz();
extern void          bioD7_skip_2d();
extern void          bioR7_2dzz();
extern void          bioR7_skip_2d();

extern void          bioD7_2dzz_char();
extern void          bioR7_2dzz_char();
extern void          bioD7_skip_2d_char();
extern void          bioR7_skip_2d_char();

extern void          bioD5_2dzz();
extern void          bioD5_skip_2d();
extern void          bioR5_2dzz();
extern void          bioR5_skip_2d();

extern void          bioD5_2dzz_char();
extern void          bioR5_2dzz_char();
extern void          bioD5_skip_2d_char();
extern void          bioR5_skip_2d_char();

extern void          bioD3_2dzz();
extern void          bioD3_skip_2d();
extern void          bioR3_2dzz();
extern void          bioR3_skip_2d();

extern void          bioD3_2dzz_char();
extern void          bioR3_2dzz_char();
extern void          bioD3_skip_2d_char();
extern void          bioR3_skip_2d_char();


extern void          bioD1_2dzz();
extern void          bioD1_skip_2d();
extern void          bioR1_2dzz();
extern void          bioR1_skip_2d();

extern void          bioD1_2dzz_char();
extern void          bioR1_2dzz_char();
extern void          bioD1_skip_2d_char();
extern void          bioR1_skip_2d_char();






extern void          bioD9_2d();
extern void          bioD9_skip_2d();
extern void          bioR9_2d();
extern void          bioR9_skip_2d();



extern void          bioD9_2d_char();
extern void          bioR9_2d_char();
extern void          bioD9_skip_2d_char();
extern void          bioR9_skip_2d_char();

extern void          bioD7_2d();
extern void          bioD7_skip_2d();
extern void          bioR7_2d();
extern void          bioR7_skip_2d();

extern void          bioD7_2d_char();
extern void          bioR7_2d_char();
extern void          bioD7_skip_2d_char();
extern void          bioR7_skip_2d_char();

extern void          bioD5_2d();
extern void          bioD5_skip_2d();
extern void          bioR5_2d();
extern void          bioR5_skip_2d();

extern void          bioD5_2d_char();
extern void          bioR5_2d_char();
extern void          bioD5_skip_2d_char();
extern void          bioR5_skip_2d_char();

extern void          bioD3_2d();
extern void          bioD3_skip_2d();
extern void          bioR3_2d();
extern void          bioR3_skip_2d();

extern void          bioD3_2d_char();
extern void          bioR3_2d_char();
extern void          bioD3_skip_2d_char();
extern void          bioR3_skip_2d_char();

extern void          bioD1_2d();
extern void          bioD1_skip_2d();
extern void          bioR1_2d();
extern void          bioR1_skip_2d();

extern void          bioD1_2d_char();
extern void          bioR1_2d_char();
extern void          bioD1_skip_2d_char();
extern void          bioR1_skip_2d_char();



extern void          bioR5_2d();




extern void          bioR9_skip_2d();


extern void          bioD7_skip_2dzz();

extern void          bioR7_skip_2dzz();

extern void          bioD5_skip_2d();
extern void          bioD5_skip_2dzz();
extern void          bioR5_skip_2d();
extern void          bioR5_skip_2dzz();





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



extern void          bioD9_2d_y();
extern void          bioD9_skip_y();

extern void          bioR9_2d_y();
extern void          bioR9_skip_2d_y();



extern void          bioD9_2d_y_char();
extern void          bioR9_2d_y_char();
extern void          bioD9_skip_2d_y_char();
extern void          bioR9_skip_2d_y_char();
extern void          bioD7_2d_y();
extern void          bioD7_skip_2d_y();
extern void          bioR7_2d_y();
extern void          bioR7_skip_2d_y();

extern void          bioD7_2d_y_char();
extern void          bioR7_2d_y_char();
extern void          bioD7_skip_2d_y_char();
extern void          bioR7_skip_2d_y_char();

extern void          bioD5_2d_y();
extern void          bioD5_skip_2d_y();
extern void          bioR5_2d_y();
extern void          bioR5_skip_2d_y_skip();

extern void          bioD5_2d_y_char();
extern void          bioR5_2d_y_char();
extern void          bioD5_skip_2d_y_char();
extern void          bioR5_skip_2d_y_char();

extern void          bioD3_2d_y();
extern void          bioD3_skip_2d_y();
extern void          bioR3_2d_y();
extern void          bioR3_skip_2d_y();

extern void          bioD3_2d_y_char();
extern void          bioR3_2d_y_char();
extern void          bioD3_skip_2d_y_char();
extern void          bioR3_skip_2d_y_char();

extern void          bioD1_2d_y();
extern void          bioD1_skip_2d_y();
extern void          bioR1_2d_y();
extern void          bioR1_skip_2d_y();

extern void          bioD1_2d_y_char();
extern void          bioR1_2d_y_char();
extern void          bioD1_skip_2d_y_char();
extern void          bioR1_skip_2d_y_char();

/*new in declararion file: */
extern void   bioD1__3d();
extern void   bioD1__2d_mult();
extern void   bioD1__2d_y();
extern void   bioD1_skip_2d_y();


extern void   bioR1__3d();
extern void   bioR1__2d_mult();
extern void   bioR1__2d_y();
extern void   bioR1_skip_2d_y();


extern void bio_2d_y_premult();
extern void bio_2d_y_postmult();


extern  void bio_2d_premult(FLOAT * , /*vector*/
					FLOAT * , /*HH*/
					FLOAT * , /*HL*/
					FLOAT * , /*LH*/
					FLOAT *, /*LL*/
					int ,  /*xlength*/
					int , /*ylength*/
					int , /*rowinc*/
					int , /*HHrowinc*/
					int , /*HLrowinc*/
					int , /*LHrowinc*/
					int , /*LLrowinc */
					int , /*Rotlevels*/
					FLOAT *, /* X*/
					int  /*ifnotAllskip*/
						   );

extern  void bio_2d_postmult(FLOAT * , /*vector*/
					  int   , /*xlength*/
					  int  , /*ylength*/
					  int  , /*rowinc*/
					  int  , /*Rotlevels*/
					  FLOAT * , /*X*/
					  int ,   /*pairity*/
					  int   /*ifnotAllskip*/
	  );


extern void bioD1_2dskip_mult(FLOAT * , /*vector*/
							 FLOAT * , /*LL */
							 int  , /*xlength*/
							  int , /*ylength*/
							  FLOAT  , /*x0*/
							  int  , /*vrowinc*/
							  int   /*LLrowinc*/
						);
extern void bioR1_2dskip_mult( FLOAT * , /*LL */
							  FLOAT * , /*vector*/
							  int  , /*xlength*/
							  int , /*ylength*/
							  FLOAT  , /*x0*/
							  int  , /*vrowinc*/
							  int   /*LLrowinc*/
						);

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




extern void      covectorpeano();
extern void      makepeanolist();
extern void      ZZtopeano();
extern void      bitspread();
extern void      getfilter();
extern void      getadjointfilter();



int wavelet_decompose3( FLOAT*, /*inspacevector*/
			int, /*Xlength*/
			int, /* Ylength*/
			int, /*Zlength, use value 1 in 2dim */
			char,/* Filterlength*/
			char, /* Levels */
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
			 char , /* Skip  */
			 char    /* ifnotSilent*/
			 );

int wavelet_invadjoint3( FLOAT*, /*inspacevector*/
		       int, /*Xlength*/
		       int, /* Ylength*/
		       int, /*Zlength, use value 1 in 2dim */
		       char,/* Filterlength*/
		       char, /* Levels */
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
			      char , /* Skip  */
			      char    /* ifnotSilent*/
			      );



int wavelet_boxrequest3(int * ,//Inbreqvector,
			int ,//Xlength,
			int ,//Ylength,
			int ,//Zlength,
			char, //Filterlength,
			char ,//Levels,
			char ,//minZLevels, 
			char ,//MaxZLevels, 
			char ,//minXYLevels, 
			char ,//MaxXYLevels, 
			char ,//Skip ,
                        int * ,//bsizevector,
			/*	int *colength_ptr,                        
				int *breqcovector,*/			
			int *  ,//request_lengthptr,
			unsigned int **  ,  //reqlist_ptr,
                        int **  ,//reqlistindex_ptr, 
			unsigned char * ,//Tempvector,
			char //ifnotSilent
			);						



int wavelet_recreq3(FLOAT *,//reqcovector,
		    int *,//reqlistindex,
                    int *,//bsizevector,
		    FLOAT *,//outvector,
		    int ,//Xlength,
		    int ,//Ylength,
		    int ,//Zlength, 
		    char ,//Filterlength,
		    char ,//Levels,
			 char , /*minZLevels */
			 char , /* maxZLevels,   */
			 char , /*minXYLevels, use  value of Levels  in 2dim */
			 char , /* MaxXYLevels,   use value of Levels in 2dim */
		    char ,//Skip , 
		    FLOAT* ,// LLLTempvector,
		    char //  ifnotSilent
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











