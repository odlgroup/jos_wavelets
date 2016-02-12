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

extern void          adjbioD9_2dzz();
extern void          adjbioD9_skip_2d();
extern void          adjbioR9_2dzz();
extern void          adjbioR9_skip_2d();



extern void          adjbioD9_2dzz_char();
extern void          adjbioR9_2dzz_char();
extern void          adjbioD9_skip_2d_char();
extern void          adjbioR9_skip_2d_char();

extern void          adjbioD7_2dzz();
extern void          adjbioD7_skip_2d();
extern void          adjbioR7_2dzz();
extern void          adjbioR7_skip_2d();

extern void          adjbioD7_2dzz_char();
extern void          adjbioR7_2dzz_char();
extern void          adjbioD7_skip_2d_char();
extern void          adjbioR7_skip_2d_char();

extern void          adjbioD5_2dzz();
extern void          adjbioD5_skip_2d();
extern void          adjbioR5_2dzz();
extern void          adjbioR5_skip_2d();

extern void          adjbioD5_2dzz_char();
extern void          adjbioR5_2dzz_char();
extern void          adjbioD5_skip_2d_char();
extern void          adjbioR5_skip_2d_char();

extern void          adjbioD3_2dzz();
extern void          adjbioD3_skip_2d();
extern void          adjbioR3_2dzz();
extern void          adjbioR3_skip_2d();

extern void          adjbioD3_2dzz_char();
extern void          adjbioR3_2dzz_char();
extern void          adjbioD3_skip_2d_char();
extern void          adjbioR3_skip_2d_char();


extern void          adjbioD1_2dzz();
extern void          adjbioD1_skip_2d();
extern void          adjbioR1_2dzz();
extern void          adjbioR1_skip_2d();

extern void          adjbioD1_2dzz_char();
extern void          adjbioR1_2dzz_char();
extern void          adjbioD1_skip_2d_char();
extern void          adjbioR1_skip_2d_char();






extern void          adjbioD9_2d();
extern void          adjbioD9_skip_2d();
extern void          adjbioR9_2d();
extern void          adjbioR9_skip_2d();



extern void          adjbioD9_2d_char();
extern void          adjbioR9_2d_char();
extern void          adjbioD9_skip_2d_char();
extern void          adjbioR9_skip_2d_char();

extern void          adjbioD7_2d();
extern void          adjbioD7_skip_2d();
extern void          adjbioR7_2d();
extern void          adjbioR7_skip_2d();

extern void          adjbioD7_2d_char();
extern void          adjbioR7_2d_char();
extern void          adjbioD7_skip_2d_char();
extern void          adjbioR7_skip_2d_char();

extern void          adjbioD5_2d();
extern void          adjbioD5_skip_2d();
extern void          adjbioR5_2d();
extern void          adjbioR5_skip_2d();

extern void          adjbioD5_2d_char();
extern void          adjbioR5_2d_char();
extern void          adjbioD5_skip_2d_char();
extern void          adjbioR5_skip_2d_char();

extern void          adjbioD3_2d();
extern void          adjbioD3_skip_2d();
extern void          adjbioR3_2d();
extern void          adjbioR3_skip_2d();

extern void          adjbioD3_2d_char();
extern void          adjbioR3_2d_char();
extern void          adjbioD3_skip_2d_char();
extern void          adjbioR3_skip_2d_char();

extern void          adjbioD1_2d();
extern void          adjbioD1_skip_2d();
extern void          adjbioR1_2d();
extern void          adjbioR1_skip_2d();

extern void          adjbioD1_2d_char();
extern void          adjbioR1_2d_char();
extern void          adjbioD1_skip_2d_char();
extern void          adjbioR1_skip_2d_char();



extern void          adjbioR5_2d();




extern void          adjbioR9_skip_2d();


extern void          adjbioD7_skip_2dzz();

extern void          adjbioR7_skip_2dzz();

extern void          adjbioD5_skip_2d();
extern void          adjbioD5_skip_2dzz();
extern void          adjbioR5_skip_2d();
extern void          adjbioR5_skip_2dzz();





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



extern void          adjbioD9_2d_y();
extern void          adjbioD9_skip_y();

extern void          adjbioR9_2d_y();
extern void          adjbioR9_skip_2d_y();



extern void          adjbioD9_2d_y_char();
extern void          adjbioR9_2d_y_char();
extern void          adjbioD9_skip_2d_y_char();
extern void          adjbioR9_skip_2d_y_char();
extern void          adjbioD7_2d_y();
extern void          adjbioD7_skip_2d_y();
extern void          adjbioR7_2d_y();
extern void          adjbioR7_skip_2d_y();

extern void          adjbioD7_2d_y_char();
extern void          adjbioR7_2d_y_char();
extern void          adjbioD7_skip_2d_y_char();
extern void          adjbioR7_skip_2d_y_char();

extern void          adjbioD5_2d_y();
extern void          adjbioD5_skip_2d_y();
extern void          adjbioR5_2d_y();
extern void          adjbioR5_skip_2d_y_skip();

extern void          adjbioD5_2d_y_char();
extern void          adjbioR5_2d_y_char();
extern void          adjbioD5_skip_2d_y_char();
extern void          adjbioR5_skip_2d_y_char();

extern void          adjbioD3_2d_y();
extern void          adjbioD3_skip_2d_y();
extern void          adjbioR3_2d_y();
extern void          adjbioR3_skip_2d_y();

extern void          adjbioD3_2d_y_char();
extern void          adjbioR3_2d_y_char();
extern void          adjbioD3_skip_2d_y_char();
extern void          adjbioR3_skip_2d_y_char();

extern void          adjbioD1_2d_y();
extern void          adjbioD1_skip_2d_y();
extern void          adjbioR1_2d_y();
extern void          adjbioR1_skip_2d_y();

extern void          adjbioD1_2d_y_char();
extern void          adjbioR1_2d_y_char();
extern void          adjbioD1_skip_2d_y_char();
extern void          adjbioR1_skip_2d_y_char();

/*new in declararion file: */
extern void   adjbioD1__3d();
extern void   adjbioD1__2d_mult();
extern void   adjbioD1__2d_y();
extern void   adjbioD1_skip_2d_y();


extern void   adjbioR1__3d();
extern void   adjbioR1__2d_mult();
extern void   adjbioR1__2d_y();
extern void   adjbioR1_skip_2d_y();


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
					   char, /*Skip */ 
					   FLOAT* ,/*covector*/			
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
		    char ,//MaxZLevels, 
		    char ,//minXYLevels, 
		    char ,//MaxXYLevels, 
		    char ,//maxXYLevels, 
		    char ,//Skip , 
		    FLOAT* ,// LLLTempvector,
		    char //  ifnotSilent
		    );



 




































