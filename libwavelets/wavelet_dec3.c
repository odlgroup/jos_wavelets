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



/*file decompose3m.c */
/* file defining procedure for decomposing image into a bitstream
 using biorthogonal filters and bitcoding.


 
                                   Algorithm and code developped 
                                               by
                                      Jan-Olov Stromberg

                                KTH, Stockholm, Sweden
                     Fast Mathematical Algorithms&Hardware, Hamden, CT, USA


                             Preliminay version by Octobler 19, 1998

*/
#include <stdio.h> 
#include <stdlib.h>
#include <math.h> 
#include <time.h>
#include <string.h> 
#define TRIG_UNSIGNEDCHAR
#include "bio_parameters.h" 
/*#include "coding_Color.h"*/
#include "bio.h"
/*
#ifndef FLOAT
#define FLOAT float
#endif
*/

/*#define FULLIMAGE   0,0,0,0,0*/
#define FULLIMAGE   0 ,0 ,0, 0,0
#define FLLSKIPIMAGE  0,0
#define XPARTIALLENGTH  (xlength)
#define YPARTIALLENGTH  (ylength)
#define XFULLENGTH   (xlength)
#define YFULLENGTH   (ylength)


#define PRINT   0
#define NEW 1
/************************/




int wavelet_decompose3(FLOAT *inspacevector,
		       int Xlength,
		       int Ylength,
		       int Zlength,
		       char Filterlength,
		       char Levels,
		       char Skip ,
		       FLOAT *covector,			                      
		       int *colength_ptr, 
		       char ifnotSilent
				   )						



	 
{
int ifnotallskip=1;
FLOAT  *vector=NULL;
int Size;
int mxlength,mylength,mzlength;
int scale,xlength,ylength,zlength,maxscale; 
/*int finalxlength,finalylength,finalzlength;*/
int LLLsubsize;
int LLHsubsize;
int LHLsubsize;
int LHHsubsize;
int HLHsubsize;
int HHLsubsize;
int HHHsubsize;
int HLLsubsize;
FLOAT *LLL=NULL,*LHH=NULL,*LLH=NULL,*LHL=NULL;
FLOAT *HLL=NULL,*HHH=NULL,*HLH=NULL,*HHL=NULL;

FLOAT *LLLv=NULL;


int skipsize=0;
void (*bioD_3d)();
void (*bioD_skip_3d)();
void (*bioD_3d_char)();
void (*bioD_skip_3d_char)();

/*   Following restrictions to the code are made

     2015-12-09    Jan-Olov Stömberg     */

 /*******************/


getfilter(&bioD_3d,&bioD_skip_3d,
	  &bioD_3d_char,&bioD_skip_3d_char,Filterlength,7);



xlength = ((Xlength-1)>>Skip)+1;
ylength = ((Ylength-1)>>Skip)+1;
zlength = ((Zlength-1)>>Skip)+1;





maxscale=Levels-1;
scale=0;

 Size=Xlength*Ylength*Zlength;



LLL=inspacevector;
LLLsubsize =  zlength*ylength*xlength;

while((scale<=maxscale)&&(scale<Skip)){
  if(scale==maxscale) LLL=covector;
      xlength=1+((Xlength-1)  >>scale);
      ylength=1+((Ylength-1) >>scale);
      zlength=1+((Zlength-1) >>scale);
                 
      if(!scale)   bioD_skip_3d_char
	(inspacevector,HHH,HHL,HLH,HLL,LHH,LHL,LLH,LLL,xlength,ylength,zlength,
	 0);
      else   bioD_skip_3d
	(vector,HHH,HHL,HLH,HLL,LHH,LHL,LLH,LLL,xlength,ylength,zlength,
	 0);
      vector=LLL;
      scale++;
    }
    skipsize=xlength*ylength*zlength;
 


  LLLv=covector;    

while(scale<=maxscale){

  xlength=1+((Xlength-1)  >>scale);
  ylength=1+((Ylength-1) >>scale);
  zlength=1+((Zlength-1) >>scale);

  mxlength=(xlength+1)>>1;
  mylength=(ylength+1)>>1;
  mzlength=(zlength+1)/2;


  vector=LLL;



  LLLsubsize =  mzlength*mylength*mxlength;
  LLHsubsize =  mzlength*mylength*(xlength-mxlength);
  LHLsubsize =mzlength*(ylength-mylength)*mxlength;
  LHHsubsize = mzlength*(ylength-mylength)*(xlength-mxlength);
  HLLsubsize =(zlength-mzlength)*mylength*mxlength;
  HLHsubsize =(zlength-mzlength)*mylength*(xlength-mxlength);
  HHLsubsize =(zlength-mzlength)*(ylength-mylength)*mxlength;
  HHHsubsize =(zlength-mzlength)*(ylength-mylength)*(xlength-mxlength);

#ifdef COARSEEND
  HHH=LLLv;
  HHL=HHH+HHHsubsize;
  HLH=HHL+HHLsubsize;
  HLL=HLH+HLHsubsize;
  LHH=HLL+HLLsubsize;
  LHL=LHH+LHHsubsize;
  LLH=LHL+LHLsubsize;
  LLLv=LLH+LLHsubsize;
 if(scale==maxscale)LLL=LLLv;
#else
 if(scale==maxscale)LLL=covector;
  LLH =covector+LLLsubsize;
  LHL=LLH+LLHsubsize;
  LHH=LHL+LHLsubsize;
  HLL=LHH+LHHsubsize;
  HLH=HLL+HLLsubsize;
  HHL=HLH+HLHsubsize;
  HHH=HHL+HHLsubsize;
#endif



 
if((!scale)){
    bioD_3d_char  (inspacevector,HHH,HHL,HLH,HLL,LHH,LHL,LLH,LLL,
                    xlength,ylength,zlength,ifnotallskip);
  }else{
    bioD_3d  (vector,HHH,HHL,HLH,HLL,LHH,LHL,LLH,LLL,
	      xlength,ylength,zlength,ifnotallskip);
  }
 

 scale++;
 vector=LLL;
   }      
 if(LLL!=LLLv)memcpy(LLLv,LLL,LLLsubsize*sizeof(FLOAT));

/***********************************************/
*colength_ptr=skipsize;

return 0;
 }   


/************************/

























































































