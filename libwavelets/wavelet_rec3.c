#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#define TRIG_UNSIGNEDCHAR
#include "bio_parameters.h"
#include "bio.h"
#ifndef FLOAT
#define FLOAT float
#endif

#define max(s,t) (s>t?s:t)
#define min(s,t) (s<t?s:t)
#define FULLIMAGE 0,0,0,0,0


int wavelet_reconstruct3(FLOAT *reccovector,
                         int colength,
                         FLOAT  *outvector,
			 int Xlength,
			 int Ylength,
			 int Zlength, 
			 char Filterlength,
			 char Levels,
			 char Skip , 
			 char     ifnotSilent
						 )




{
  //long time0=0;
  //long time1=1;
FLOAT *vector=NULL;
int scale,xlength,ylength,zlength,maxscale;
int mxlength,mylength,mzlength;

FLOAT *LLL=NULL,*LHH=NULL,*LLH=NULL,*LHL=NULL;
FLOAT *HLL=NULL,*HHH=NULL,*HLH=NULL,*HHL=NULL;
 FLOAT  *LLLv=NULL;
int skipsize;
int Size;
int LLLsubsize;
int LLHsubsize;
int LHLsubsize;
int LHHsubsize;
int HLHsubsize;
int HHLsubsize;
int HHHsubsize;
int HLLsubsize;
int finalxlength,finalylength,finalzlength;
int finalLLLsize;




int ifnotallskip=1;
 
void  (*bioR_3d)();
void  (*bioR_skip_3d)();
void  (*bioR_3d_char)();
void  (*bioR_skip_3d_char)();

/*   Following restrictions to the code are made

     2015-12-09    Jan-Olov Stömberg     */


 /*******************/

getfilter(&bioR_3d,&bioR_skip_3d,
	  &bioR_3d_char,&bioR_skip_3d_char,-Filterlength,7);

 Skip=(Skip<Levels?Skip:Levels);
 Size=Xlength*Ylength*Zlength;
xlength=1+((Xlength-1)>>Skip);
ylength=1+((Ylength-1)>>Skip);
zlength=1+((Zlength-1)>>Skip);
skipsize= xlength*ylength*zlength;

maxscale=Levels-1;

finalxlength= 1+((Xlength-1)>>(maxscale+1));
finalylength= 1+((Ylength-1)>>(maxscale+1));
finalzlength=1+((Zlength-1) >>(maxscale+1));
finalLLLsize=finalxlength*finalylength*finalzlength;




if(Levels==0)skipsize=Xlength*Ylength*Zlength;

vector = outvector;
 if(maxscale< 0){
 memcpy(outvector,reccovector,Size*sizeof(FLOAT));

return 0;
 }
scale=maxscale;

LLLv =reccovector + skipsize ;



xlength=Xlength;
ylength=Ylength;
zlength=Zlength;
 LLLsubsize =  zlength*ylength*xlength;
/********************************************/





/************************************************/

 while(scale >=Skip){  /* over scales */

  xlength=1+((Xlength-1)  >>scale);
  ylength=1+((Ylength-1) >>scale);
  zlength=1+((Zlength-1) >>scale);
   
  mxlength=(xlength+1)>>1;
  mylength=(ylength+1)>>1;
  mzlength=(zlength+1)>>1;

  vector=outvector;
 
  LLLsubsize =  mzlength*mylength*mxlength;
  LLHsubsize =  mzlength*mylength*(xlength-mxlength);
  LHLsubsize =mzlength*(ylength-mylength)*mxlength;
  LHHsubsize = mzlength*(ylength-mylength)*(xlength-mxlength);
  HLLsubsize =(zlength-mzlength)*mylength*mxlength;
  HLHsubsize =(zlength-mzlength)*mylength*(xlength-mxlength);
  HHLsubsize =(zlength-mzlength)*(ylength-mylength)*mxlength;
  HHHsubsize =(zlength-mzlength)*(ylength-mylength)*(xlength-mxlength);

#ifdef COARSEEND
  if(scale==maxscale)LLL=LLLv - finalLLLsize;
  LLH = LLLv - LLLsubsize - LLHsubsize;
  LHL = LLH - LHLsubsize;
  LHH = LHL - LHHsubsize;
  HLL = LHH -    HLLsubsize;
  HLH = HLL -  HLHsubsize;
  HHL = HLH- HHLsubsize;
  HHH = HHL-   HHHsubsize;
#else
  if(scale==maxscale)LLL=reccovector;
  LLH = reccovector + LLLsubsize;
  LHL = LLH + LLHsubsize;
  LHH = LHL + LHLsubsize;
  HLL = LHH + LHHsubsize;
  HLH = HLL + HLLsubsize;
  HHL = HLH + HLHsubsize;
  HHH = HHL + HHLsubsize;
#endif  
 

   
  ifnotallskip=(scale<Skip?0:1);


  if(scale<maxscale)LLL=outvector;

  zlength=1+((Zlength-1) >>scale);

 
/************************************************/   
 vector=outvector;

if(scale==0){
bioR_3d_char(HHH,HHL,HLH,HLL,LHH,LHL,LLH,LLL,
                   outvector,xlength,ylength,zlength,
                   ifnotallskip);
}else{
 bioR_3d(HHH,HHL,HLH,HLL,LHH,LHL,LLH,LLL,
                   vector,xlength,ylength,zlength,
                   ifnotallskip);
}
 LLL=vector;
  scale--;
}
 LLL=outvector;

 while(scale>=0){
   ifnotallskip=0;
  xlength=1+((Xlength-1)  >>scale);
  ylength=1+((Ylength-1) >>scale);
  zlength=1+((Zlength-1) >>scale); 
	 if(scale==0){
	   bioR_3d_char(HHH,HHL,HLH,HLL,LHH,LHL,LLH,LLL,
					outvector,xlength,ylength,zlength,
					ifnotallskip);
	 }else{
	   bioR_3d(HHH,HHL,HLH,HLL,LHH,LHL,LLH,LLL,
			   vector,xlength,ylength,zlength,
			   ifnotallskip);
	 }
  
 
LLL=vector;	 
	 scale--;
 }

 
 return 0;
}

/*********************/





















