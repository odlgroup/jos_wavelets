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


      


int wavelet_adjoint3(FLOAT *reccovector,
                         int colength,
                         FLOAT  *outvector,
			 int Xlength,
			 int Ylength,
			 int Zlength, 
			 char Filterlength,
			 char Levels,
			 char minZLevels, 
			 char maxZLevels, 
			 char minXYLevels, 
			 char maxXYLevels, 
			 char Skip , 
			 char     ifnotSilent
						 )




{
  //long time0=0;
  //long time1=1;
FLOAT *vector=NULL;
int scale,xlength,ylength,zlength,maxscale;
int  xylength,maxzscale; 
int mxlength,mylength,mzlength;
FLOAT *LLL=NULL,*LHH=NULL,*LLH=NULL,*LHL=NULL;
FLOAT *HLL=NULL,*HHH=NULL,*HLH=NULL,*HHL=NULL;
FLOAT  *LLLv=NULL,*LHHv=NULL,*LLHv=NULL,*LHLv=NULL;
FLOAT *HLLv=NULL;
FLOAT *LLLv0=NULL;
int skipsize;
int Size;
int finalxlength,finalylength,finalzlength;
int finalLLLsize;
int LLLsubsize;
int LLHsubsize;
int LHLsubsize;
int LHHsubsize;
int HLHsubsize;
int HHLsubsize;
int HHHsubsize;
int HLLsubsize;
int ifnotallskip=1;
void  (*adjbioD_3d)();
void  (*adjbioD_skip_3d)();
void  (*adjbioD_3d_char)();
void  (*adjbioD_skip_3d_char)();




 

/*   Following restrictions to the code are made                                                                  
                                                                                                                  
     2015-12-09    Jan-Olov Stömberg     */

 maxXYLevels=Levels;
 maxZLevels=Levels;
 minZLevels=0;
 minXYLevels=0;

 /*******************/



getadjointfilter(&adjbioD_3d,&adjbioD_skip_3d,
	  &adjbioD_3d_char,&adjbioD_skip_3d_char,Filterlength,7);



 Size=Xlength*Ylength*Zlength;
xlength=1+((Xlength-1)>>Skip);
ylength=1+((Ylength-1)>>Skip);
if(Skip<maxZLevels)zlength=1+((Zlength-1)>>Skip);
else zlength=1+((Zlength-1)>>maxZLevels);
skipsize= xlength*ylength*zlength;

maxscale=Levels-1;
maxzscale= maxZLevels-1;

if(Levels==0)skipsize=Xlength*Ylength*Zlength;

vector = outvector;
if(((maxscale< 0)&&(minZLevels==0))||((maxZLevels==0)&&(maxXYLevels==0))) return 0;
scale=maxscale;

finalxlength= 1+((Xlength-1)>>(maxscale+1));
finalylength= 1+((Ylength-1)>>(maxscale+1));
 if(maxscale>maxzscale) finalzlength=1+((Zlength-1) >>(maxzscale+1));
 else  finalzlength=1+((Zlength-1) >>(maxscale+1));
finalLLLsize=finalxlength*finalylength*finalzlength;


LLLv =reccovector + skipsize ;
LLL=LLLv - finalLLLsize;
LLLv0=outvector+finalLLLsize;




xlength=Xlength;
ylength=Ylength;
zlength=Zlength;


/************************************************/

 while(scale >=0*Skip){  /* over scales */

  xlength=1+((Xlength-1)  >>scale);
  ylength=1+((Ylength-1) >>scale);

  zlength=1+((Zlength-1) >>scale);
   
  mxlength=(xlength+1)>>1;
  mylength=(ylength+1)>>1;
 
mzlength=(zlength+1)/2;


  vector=outvector;
  LLLsubsize =  mzlength*mylength*mxlength;
  LLHsubsize =  mzlength*mylength*(xlength-mxlength);
  LHLsubsize =mzlength*(ylength-mylength)*mxlength;
  LHHsubsize = mzlength*(ylength-mylength)*(xlength-mxlength);
  HLLsubsize =(zlength-mzlength)*mylength*mxlength;
  HLHsubsize =(zlength-mzlength)*mylength*(xlength-mxlength);
  HHLsubsize =(zlength-mzlength)*(ylength-mylength)*mxlength;
  HHHsubsize =(zlength-mzlength)*(ylength-mylength)*(xlength-mxlength);

 



  if(minXYLevels<=minZLevels){ 
   LLH = LLLv - LLLsubsize - LLHsubsize;
  LHL = LLH - LHLsubsize;
  LHH = LHL - LHHsubsize;
	HLL  = LHH -    HLLsubsize;
	HLH   = HLL -  HLHsubsize;
  }

    HHL= HLH- HHLsubsize;
	HHH =  HHL-   HHHsubsize;

 
 HLLv = HLL;
   LLHv = LLH ;
   LHLv = LHL ;
   LHHv = LHH  ;

 LLLv0 = outvector+LLLsubsize;

   
  ifnotallskip=(scale<Skip?0:1);
    
  xylength=mylength*(xlength-mxlength);
/*********************************
	

 ******************************/

if(scale<maxscale)LLL=outvector;

 
 
    xylength=((ylength+1)>>1)*((xlength)>>1);
 



/************************************************/
    
/********************************************/
  xylength=((ylength)>>1)*((xlength+1)>>1);  
   

   
/************************************************/
    
/********************************************/
	     xylength=((ylength)>>1)*(xlength>>1);
 		
	

zlength=1+((Zlength-1) >>scale);
 


 
 
/************************************************/
   

 

 vector=outvector;
if((scale<=maxzscale)&&(scale <maxXYLevels)){
 zlength=1+((Zlength-1) >>scale); 
  
if(scale==0){
adjbioD_3d_char(HHH,HHL,HLH,HLL,LHH,LHL,LLH,LLL,
                   outvector,xlength,ylength,zlength,
                   ifnotallskip);
}else{
 adjbioD_3d(HHH,HHL,HLH,HLL,LHH,LHL,LLH,LLL,
                   vector,xlength,ylength,zlength,
                   ifnotallskip);


}

 }



 LLL=vector;
  scale--;
 }

 if((scale==maxscale)&&(scale >=  minZLevels-1))LLL=reccovector;
 else LLL=outvector;


 while(scale>=0){
   ifnotallskip=0;
  xlength=1+((Xlength-1)  >>scale);
  ylength=1+((Ylength-1) >>scale);

   

 zlength=1+((Zlength-1) >>scale);
   /************************************************/
   
   vector=outvector;
   LLH=NULL;
   LHL=NULL;
   LHH=NULL;
   HLL=NULL;
   HLH=NULL;
   HHL=NULL;
   HHH=NULL;
   if((scale<=maxzscale)&&(scale <maxXYLevels)){
	 zlength=1+((Zlength-1) >>scale); 

	 if(scale==0){
	   adjbioD_3d_char(HHH,HHL,HLH,HLL,LHH,LHL,LLH,LLL,
					outvector,xlength,ylength,zlength,
					ifnotallskip);
	 }else{
	   adjbioD_3d(HHH,HHL,HLH,HLL,LHH,LHL,LLH,LLL,
			   vector,xlength,ylength,zlength,
			   ifnotallskip);
	 }
   }
 
LLL=vector;	 
	 scale--;
 }

 
 return 0;
}

/*********************/





















