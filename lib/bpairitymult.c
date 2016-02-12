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


#include "bio_parameters.h"

void bpairitymultiply(FLOAT *vector,
		      int xlength,
		      int ylength,
		      int zlength,
                      int pairity,
		      FLOAT factor){
  int j;
  int k;


  if((pairity==0)&&(xlength>1))
    for(k=0;k<zlength;k++)
      for(j=0;j<ylength;j++) vector[xlength*(k*ylength+j)] *=factor;

  if(((xlength&1)==(1-pairity))&&(xlength>1)) 
    for(k=0;k<zlength;k++)
      for(j=0;j<ylength;j++) vector[xlength*(k*ylength+j+1)-1] *=factor;

  if((pairity==0)&&(ylength>1))
     for(k=0;k<zlength;k++)
      for(j=0;j<xlength;j++) vector[xlength*k*ylength +j] *=factor;

  if(((ylength&1)==(1-pairity))&&(ylength>1)) 
    for(k=0;k<zlength;k++)
      for(j=0;j<xlength;j++) vector[j+xlength*((k+1)*ylength-1)] *=factor;

  if((pairity==0)&&(zlength>1))
    for(k=0;k<ylength;k++)
      for(j=0;j<xlength;j++) vector[(xlength*k)+j] *=factor;

  if(((zlength&1)==(1-pairity))&&(zlength>1)) 
    for(k=0;k<ylength;k++)
      for(j=0;j<xlength;j++) vector[xlength*(k+(zlength-1)*ylength)+j] *=factor;
}




void boundarymultiply(FLOAT *vector,
		      int xlength,
		      int ylength,
		      int zlength,
               		FLOAT factor){

  int k;
  int j;

  for(k=0;k<ylength;k++)
    for(j=0;j<xlength;j++) vector[k*xlength+j] *=factor;
  for(k=0;k<zlength;k++){
    for(j=0;j<xlength;j++) vector[xlength*k*ylength +j] *=factor;   
    for(j=0;j<ylength;j++){
      vector[xlength*(k*ylength+j)] *=factor;
      vector[xlength*(k*ylength+j+1)-1] *=factor;
    }
    for(j=0;j<xlength;j++) vector[xlength*((k+1)*ylength -1) +j] *=factor;     }
  for(k=0;k<ylength;k++)
    for(j=0;j<xlength;j++) vector[((zlength-1)*ylength+k)*xlength+j] *=factor;
}



