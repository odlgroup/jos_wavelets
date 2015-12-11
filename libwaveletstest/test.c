#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>



#define XLength 470
#define YLength 500   // has to be 1 for 1D test
#define ZLength 400   // has to be 1 for 1D  and 2D  test

 /* makefile for test  has an include file;  ../build/makefile.setting
     where "HIGH_PRECISION"   (is/is not)  defined  */

#define HIGH_PRECISION 1

#ifndef HIGH_PRECISION 
#define FLOAT float
#else 
#define FLOAT double
#endif

#include "../libwavelets/wavelet_transform.h"
#include "test.h"


int main()
{
  int Xlength=XLength;
  int Ylength=YLength;
  int Zlength=ZLength;
  clock_t  time0, time1,time2,time3,time4,time5;

 

    FLOAT *body,*body_out;
  FLOAT *waveletcoeff;
  unsigned int a=0;
  unsigned int mask12=0xfff;  /*  mask with 12  "ones"  */     
  FLOAT * ptr;
int size;

 FLOAT *diffvector;
FLOAT  *body_save;
double maxvalue=0;
double poserror=0;
double negerror=0;
double error2=0;
double norm2=0;
double error1=0;
double norm1=0;
double dsum=0;
double errorlog10=0;
int nonexactnumber=0;
int Size=0;
FLOAT SNR;

/* for 2D test and  1D test : */
//Zlength=1; 
/* for 1D test:  */
//Ylength=1;
size = Xlength*Ylength*Zlength;

/* Allocation of arrays used in test */

 body=(FLOAT*)malloc((unsigned int)(size*sizeof(FLOAT)));
waveletcoeff=(FLOAT*)malloc((unsigned int)(size*sizeof(FLOAT)));
body_out=(FLOAT*)malloc((unsigned int)(size*sizeof(FLOAT)));
 body_save=(FLOAT*)malloc((unsigned int)(size*sizeof(FLOAT))); 
diffvector=(FLOAT*)malloc((unsigned int)(size*sizeof(FLOAT))); 
ptr=body;

 time0=clock();    /* uses c-call clock()  for timing process */ 
 
for(int k=0;k<size;k++){
a=mask12&random();       /*generetin 12 bits unsigned integer */
//a=1.0;
 *ptr++=(FLOAT)a;       /* converting and saving on FLOAT array  */
 }
time1=clock(); 

/*  saving a backup  of the generated arry */
 memcpy(body_save,body,(size_t)size*sizeof(FLOAT));
 //for(int k=0;k<size;k++)printf("indata  %f , ",body[k]);
 printf("\n");
 time2=clock(); 

 
 /* CHOOSE dimension  1D  2D or 3D  also TRANSFORM or ADJOINT */
 /*  ERASE  comment sign  "//" for one only of following six lines
 */
  wavelet_transform3D(body,Xlength,Ylength,Zlength,waveletcoeff);
 //wavelet_transform2D(body,Xlength,Ylength,waveletcoeff);
 //wavelet_transform1D(body,Xlength,waveletcoeff);
 //adjointinvwavelet_transform3D(body,Xlength,Ylength,Zlength,waveletcoeff);
 //adjointinvwavelet_transform2D(body,Xlength,Ylength,waveletcoeff);
 //adjointinvwavelet_transform1D(body,Xlength,waveletcoeff);

 time3=clock();
 
 //for(int k=0;k<size;k++)printf("wcoeff  %f , ",waveletcoeff[k]);
 printf("\n");
 

 /* INVERSE OPERATORS  has to match operator chosen above  */
 /*  ERASE  comment sign  "//" for one only of following six lines
 */ 
 /* CHOOSE dimension  1D  2D or 3D  also TRANSFORM or ADJOINT */

 invwavelet_transform3D(waveletcoeff,Xlength,Ylength,Zlength,body_out);
 //invwavelet_transform2D(waveletcoeff,Xlength,Ylength,body_out);
 //invwavelet_transform1D(waveletcoeff,Xlength,body_out);
 //adjointwavelet_transform3D(waveletcoeff,Xlength,Ylength,Zlength,body_out);
 //adjointwavelet_transform2D(waveletcoeff,Xlength,Ylength,body_out);
 //adjointwavelet_transform1D(waveletcoeff,Xlength,body_out);

 time4=clock();
 // for(int k=0;k<size;k++)printf("outdata  %f , ",body_out[k]);
 printf("\n");

 /* collect statistics on arrays */
compimages(body_save,body_out,diffvector,size,
	      &maxvalue,&poserror,&negerror,&error2,&norm2,&error1,&norm1,
              &errorlog10, &nonexactnumber,&Size,&dsum);

 time5=clock();

 /* Printing out timings of processes done above  */

 printf("Size %d x %d x %d, random: %f sek,copy %f ,  WT  %f sek, invWT %f sek, compare %f sek\n",
	
	Xlength,Ylength,Zlength, 0.000001*(time1-time0),0.000001*(time2-time1),0.000001*(time3-time2),0.000001*(time4-time3),0.000001*(time5-time4));



 if(body!=NULL)free(body);
 if(body_save!=NULL)free(body_save);
 if(body_out!=NULL)free(body_out);
 if(waveletcoeff!=NULL)free(waveletcoeff);

 /* show statistics from the test */
show_compresults(maxvalue,poserror,negerror,error2,norm2,error1,norm1,
                   errorlog10,nonexactnumber,size,&SNR,sizeof(FLOAT),dsum);

 return 0;
}
