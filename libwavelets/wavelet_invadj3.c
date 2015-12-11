
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




int wavelet_invadjoint3(FLOAT *inspacevector,
		       int Xlength,
		       int Ylength,
		       int Zlength,
		       char Filterlength,
		       char Levels,
		       char minZLevels, 
		       char MaxZLevels, 
		       char minXYLevels, 
		       char MaxXYLevels, 
		       char Skip ,
		       FLOAT *covector,			                      
		       int *colength_ptr, 
		       char ifnotSilent
				   )						



	 
{

int ifnotallskip=1;
//long time0;
//long time1;
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

int LLHlength;
int LHHlength;
int HLLlength;
int zscale;
int xyscale;
int maxzscale;
FLOAT *LLL=NULL,*LHH=NULL,*LLH=NULL,*LHL=NULL;
FLOAT *HLL=NULL,*HHH=NULL,*HLH=NULL,*HHL=NULL;
FLOAT *LLLv=NULL,*LHHv=NULL,*LLHv=NULL,*LHLv=NULL;
FLOAT *HLLv=NULL,*HHHv=NULL,*HLHv=NULL,*HHLv=NULL;
FLOAT *H, *L;
FLOAT * HH0;
int skipsize=0;
//FLOAT Thresh;
int xylength;
int z0length;
int x0length;
int y0length;


void (*adjbioR_3d)();
void (*adjbioR_skip_3d)();
void (*adjbioR_3d_char)();
void (*adjbioR_skip_3d_char)();



 

/*   Following restrictions to the code are made                               \
                                                                                
                                                                               \
                                                                                
									       2015-12-09    Jan-Olov Stömberg     */

 MaxXYLevels=Levels;
 MaxZLevels=Levels;
 minZLevels=0;
 minXYLevels=0;

 /*******************/





/*only using first pointer the other are redefined here */
getadjointfilter(&adjbioR_3d,&adjbioR_skip_3d,
	  &adjbioR_3d_char,&adjbioR_skip_3d_char,-Filterlength,7);



xlength = ((Xlength-1)>>Skip)+1;
ylength = ((Ylength-1)>>Skip)+1;
zlength = ((Zlength-1)>>Skip)+1;





maxscale=Levels-1;
maxzscale=MaxZLevels-1;
scale=0;

 Size=Xlength*Ylength*Zlength;



LLL=inspacevector;


while((scale<=maxscale)&&(scale<Skip)){
if((scale==maxscale)&&(scale>=minZLevels-1))LLL=covector;
      xlength=1+((Xlength-1)  >>scale);
      ylength=1+((Ylength-1) >>scale);
   zlength=1+((Zlength-1) >>scale);

if(scale<=maxzscale){
      xlength=1+((Xlength-1)  >>scale);
      ylength=1+((Ylength-1) >>scale);
      zlength=1+((Zlength-1) >>scale);
       
      if(!scale)   adjbioR_skip_3d_char
	(inspacevector,HHH,HHL,HLH,HLL,LHH,LHL,LLH,LLL,xlength,ylength,zlength,
	 0);
      else   adjbioR_skip_3d
	(vector,HHH,HHL,HLH,HLL,LHH,LHL,LLH,LLL,xlength,ylength,zlength,
	 0);
      vector=LLL;
 }
      scale++;
 }


      xlength=1+((Xlength-1)  >>scale);
      ylength=1+((Ylength-1) >>scale);
   zlength=1+((Zlength-1) >>scale);

skipsize=xlength*ylength*zlength;




LLLv=covector;
    
while(scale<=maxscale){

  xlength=1+((Xlength-1)  >>scale);
  ylength=1+((Ylength-1) >>scale);
  if(scale>maxzscale) zlength=1+((Zlength-1) >>(maxzscale+1));
  else  zlength=1+((Zlength-1) >>scale);
  
  mxlength=(xlength+1)>>1;
  mylength=(ylength+1)>>1;
  if(scale>maxzscale) mzlength=zlength;
  else mzlength=(zlength+1)/2;


  vector=LLL;
  LLLsubsize =  mzlength*mylength*mxlength;
  LLHsubsize =  mzlength*mylength*(xlength-mxlength);
  LHLsubsize =mzlength*(ylength-mylength)*mxlength;
  LHHsubsize = mzlength*(ylength-mylength)*(xlength-mxlength);
  HLLsubsize =(zlength-mzlength)*mylength*mxlength;
  HLHsubsize =(zlength-mzlength)*mylength*(xlength-mxlength);
  HHLsubsize =(zlength-mzlength)*(ylength-mylength)*mxlength;
  HHHsubsize =(zlength-mzlength)*(ylength-mylength)*(xlength-mxlength);

  HHHv=LLLv;
  HHLv=HHHv+HHHsubsize;
  HLHv=HHLv+HHLsubsize;
  if(minZLevels>=minXYLevels){
  HLLv=HLHv+HLHsubsize;
  LHHv=HLLv+HLLsubsize;
  LHLv=LHHv+LHHsubsize;
  LLHv=LHLv+LHLsubsize;
  LLLv=LLHv+LLHsubsize;
  }


  HLH=HLHv;
  HHL=HHLv;
  HHH=HHHv;
  HLL=HLLv;


  LLH=LLHv;
  LHL=LHLv;
  LHH=LHHv;

if((scale==maxscale)&&(scale>=minZLevels-1)&&
     ((scale>=minXYLevels-1)||(HLLsubsize==0)))LLL=LLLv;


 if((scale<=maxzscale)&&(scale <MaxXYLevels)){  //always true
 zlength=1+((Zlength-1) >>scale); 
  z0length=((zlength+1)>>1);
  // printf("Doing scale %d here \n",scale);  
if((!scale)){
    adjbioR_3d_char  (inspacevector,HHH,HHL,HLH,HLL,LHH,LHL,LLH,LLL,
                    xlength,ylength,zlength,ifnotallskip);
  }else{
    adjbioR_3d  (vector,HHH,HHL,HLH,HLL,LHH,LHL,LLH,LLL,
	      xlength,ylength,zlength,ifnotallskip);
  }

 }


 scale++;

if(ifnotallskip){


if(minXYLevels<=minZLevels){
x0length=(xlength+1)>>1;
y0length=(ylength+1)>>1;
z0length=zlength>>1;
HLLlength=x0length*y0length*z0length;
if(HLLlength>0){
xyscale=scale;


HH0=HLLv;

 }
 }
 

	     z0length=((zlength+1)>>1);
	     xylength=((ylength)>>1)*((xlength)>>1);   
	     LHHlength =z0length*xylength;
		 if(LHHlength>0){
		   /*******************************************/
		   zscale=scale;
		   H=LHHv;
		   L=LHH; 
		   vector=LHH;		 
/***********************************************/
		 }


 z0length=((zlength+1)>>1);  
	     xylength=((ylength+1)>>1)*(xlength>>1);   
	     LLHlength =z0length*xylength;
		 if(LLHlength>0){
		   /*******************************************/
		   zscale=scale;
		   H=LLHv;
		   L=LLH;
		   vector=LLH;
		 
		   /***********************************************/
		 }


 }

 vector=LLL;
 }

zscale=scale;


 if(LLL!=LLLv)memcpy(LLLv,LLL,LLLsubsize*sizeof(FLOAT));

/***********************************************/
*colength_ptr=skipsize;

return 0;
}   


/************************/

























































































