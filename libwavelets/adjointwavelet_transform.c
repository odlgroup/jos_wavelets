# Copyright 2015-2016 Jan-Olov Stromberg <jostromb@kth.se>
#
# jos_wavelets is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# jos_wavelets is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with jos_wavelets.  If not, see <http://www.gnu.org/licenses/>.


#include "bio.h"
#ifndef HIGH_PRECISION
#define FLOAT float
#else
#define FLOAT double
#endif
//#include "../include/wavelet_transform.h"



/* NOTE that  MaxZLevels and MaxXYLevels should be at least as
Levels   in the parameterchoice, this means they have no effect,
and the only the 3d filters will be used.
the 2d and 2dy filters are not updated yet with the adjoint filters  */


int adjointinvwavelet_transform3D(FLOAT *inspacevector,
		       int xlength,
		       int ylength,
		       int zlength,
		       int Filterlength, /* 1,3,5,7 or 9 */
		       int levels,		       
     		       FLOAT *waveletcoefficients
			      			){
int co_length;
 /*
  printf("\ninspacevector5:\n");
  
  for(v=0;v<xlength*ylength*zlength;v++)printf(" %f ",inspacevector[v]);
 */

			
return wavelet_invadjoint3(inspacevector,
			   xlength,
			   ylength,
			   zlength,
			   (char) Filterlength,
			    (char)levels, 
			    (char)0, /*char Skip */
			   waveletcoefficients,
			   &co_length, 
			   0 /* char ifnotSilent: 1 for printing */
			   );
}						



int adjointwavelet_transform3D(FLOAT * waveletcoefficients,
                         int xlength,
			 int ylength,
                         int zlength,
		        int Filterlength, /* 1,3,5,7 or 9 */
		       int levels,		       
                         FLOAT  *outvector){

return  wavelet_adjoint3(waveletcoefficients,
			     xlength*ylength*zlength,/*coeffcient_length*/
                             outvector,
			     xlength,
			     ylength,
			     zlength, 
			    (char)Filterlength,
			     (char)levels,
			     (char)0,  /* Skip: general skiplevel
				    (changes coeff_length)*/ 
			     0 /*ifnotSilent  set 1 for printing */
			     );
}
 


int adjointinvwavelet_transform2D(FLOAT *inspacevector,
		       int xlength,
		       int ylength,
		        int Filterlength, /* 1,3,5,7 or 9 */
		       int levels,		       
		       FLOAT *waveletcoefficients
			      			){
int co_length;
			
return wavelet_invadjoint3(inspacevector,
			   xlength,
			   ylength,
			   1,/*zlength,*/
			   (char)Filterlength,
			   (char)levels,
			   (char)0, /*char Skip */
			   waveletcoefficients,
			   &co_length, 
			   0 /* char ifnotSilent: 1 for printing */
			   );
}						



int adjointwavelet_transform2D(FLOAT * waveletcoefficients,
                         int xlength,
			 int ylength,
		        int Filterlength, /* 1,3,5,7 or 9 */
		       int levels,		       
			 FLOAT  *outvector){

return  wavelet_adjoint3(waveletcoefficients,
			     xlength*ylength,/*coeffcient_length*/
                             outvector,
			     xlength,
			     ylength,
			     1, /* zlength*/ 
			 (char)Filterlength,
			    (char) levels, 
			     (char)0,  /* Skip: general skiplevel
				    (changes coeff_length)*/ 
			     0 /*ifnotSilent  set 1 for printing */
			     );
}
 




int adjointinvwavelet_transform1D(FLOAT *inspacevector,
		       int xlength,
		        int Filterlength, /* 1,3,5,7 or 9 */
		       int levels,		       
			FLOAT *waveletcoefficients
			      			){
int co_length;
			
return wavelet_invadjoint3(inspacevector,
			   xlength,
			  1,/*ylength*/
			   1,/*zlength,*/
			   (char)Filterlength,
			   (char)levels, 
			   (char)0, /*char Skip */
			   waveletcoefficients,
			   &co_length, 
			   0 /* char ifnotSilent: 1 for printing */
			   );
}						


int adjointwavelet_transform1D(FLOAT * waveletcoefficients,
                         int xlength,
		       int Filterlength, /* 1,3,5,7 or 9 */
			       int levels,  
			  FLOAT  *outvector){

return  wavelet_adjoint3(waveletcoefficients,
			     xlength,/*coeffcient_length*/
                             outvector,
			     xlength,
			     1, /*ylength*/
			     1, /* zlength*/ 
			    (char) Filterlength,
			    (char) levels, /*Levels of transforms*/
			     (char)0,  /* Skip: general skiplevel
				    (changes coeff_length)*/ 
			     0 /*ifnotSilent  set 1 for printing */
			     );
}
 
