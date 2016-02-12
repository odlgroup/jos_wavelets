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


/*  uncomment next line for checking reconstruction  vith inverse adjoint */
//#define  ADJOINT_TRANSF

#define DIM   3           /*  dimension 1, 2 or 3 */
#define FILTERLENGTH    9     /* 1, 3, 5, 7 or 9 */
/* FILTERLENGtH 1: decimation only  L is even decimation; 
                                    H is odd decimation */  
#define SCALES 10


#define XLength 512
#define YLength 513    /*   used only in dimension 2 and 3 */
#define ZLength 511      /*  used only in dimension 3  */


/***********DEFINITINS ABOVE MAY BE CHANGES BY USER     *****************/

 #include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>




 /* makefile for test  has an include file;  ../build/makefile.setting
     where "HIGH_PRECISION"   (is/is not)  defined  */

#define HIGH_PRECISION 1

#ifndef HIGH_PRECISION 
#define FLOAT float
#else 
#define FLOAT double
#endif

#include "../lib/wavelet_transform.h"
#include "test.h"


int main()
{

#if(DIM==3)
  int Xlength=XLength;
  int Ylength=YLength;
  int Zlength=ZLength;
#endif
#if(DIM==2)
  int Xlength=XLength;
  int Ylength=YLength;
  int Zlength=1;
#endif
#if(DIM==1)
  int Xlength=XLength;
  int Ylength=1;
  int Zlength=1;
#endif

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
 int k;
/* for 2D test and  1D test : */
//Zlength=1; 
/* for 1D test:  *///Ylength=1;

#ifndef ADJOINT_TRANSF
  printf("Test of reconstruction with symmetric biorthogonal filters \n \
         Filterlength=%d, Number of scales=%d\n\n",FILTERLENGTH,SCALES);

#else
        printf("Test of reconstruction with the adjoint filters of \n \
\t symmetric  biorthogonal filters with filterlength=%d.\n\
\t Number of scales=%d\n\n",FILTERLENGTH,SCALES);
#endif

size = Xlength*Ylength*Zlength;

/* Allocation of arrays used in test */

 body=(FLOAT*)malloc((unsigned int)(size*sizeof(FLOAT)));
waveletcoeff=(FLOAT*)malloc((unsigned int)(size*sizeof(FLOAT)));
body_out=(FLOAT*)malloc((unsigned int)(size*sizeof(FLOAT)));
 body_save=(FLOAT*)malloc((unsigned int)(size*sizeof(FLOAT))); 
 diffvector=(FLOAT*)malloc((unsigned int)(size*sizeof(FLOAT))); 
ptr=body;

 time0=clock();    /* uses c-call clock()  for timing process */ 
 
for(k=0;k<size;k++){
  a=mask12&random();       /*generetin 12 bits unsigned integer */
  //a=1.0;
 *ptr++=(FLOAT)a;       /* converting and saving on FLOAT array  */
 }
time1=clock(); 

/*  saving a backup  of the generated arry */
 memcpy(body_save,body,(size_t)size*sizeof(FLOAT));
 //for(k=0;k<size;k++)printf("indata  %f , ",body[k]);
 printf("\n");
 time2=clock(); 

 

 /* CHOOSE dimension  1D  2D or 3D  also TRANSFORM or ADJOINT */
 /*  ERASE  comment sign  "//" for one only of following six lines
 */
#ifndef ADJOINT_TRANSF

#if(DIM==3) 
wavelet_transform3D(body,Xlength,Ylength,Zlength,
 FILTERLENGTH,SCALES,waveletcoeff);
#endif

#if(DIM==2) 
 wavelet_transform2D(body,Xlength,Ylength,
		     FILTERLENGTH,SCALES,waveletcoeff);
#endif

#if(DIM==1) 
wavelet_transform1D(body,Xlength,
		    FILTERLENGTH,SCALES,waveletcoeff);
#endif

#else

#if(DIM==3) 
adjointinvwavelet_transform3D(body,Xlength,Ylength,Zlength,
				FILTERLENGTH,SCALES,waveletcoeff);

#endif

#if(DIM==2) 
adjointinvwavelet_transform2D(body,Xlength,Ylength,
			      FILTERLENGTH,SCALES,waveletcoeff);
#endif

#if(DIM==1) 
adjointinvwavelet_transform1D(body,Xlength,
  FILTERLENGTH,SCALES,waveletcoeff);
#endif

#endif

 time3=clock();
 
 //for(k=0;k<size;k++)printf("wcoeff  %f , ",waveletcoeff[k]);
 printf("\n");
 

 /* INVERSE OPERATORS  has to match operator chosen above  */
 /*  ERASE  comment sign  "//" for one only of following six lines
 */ 
 /* CHOOSE dimension  1D  2D or 3D  also TRANSFORM or ADJOINT */
#ifndef ADJOINT_TRANSF

#if(DIM==3) 
invwavelet_transform3D(waveletcoeff,Xlength,Ylength,Zlength,
		       FILTERLENGTH,SCALES,body_out);
#endif

#if(DIM==2) 
 invwavelet_transform2D(waveletcoeff,Xlength,Ylength,
			FILTERLENGTH,SCALES,body_out);
#endif

#if(DIM==1) 
invwavelet_transform1D(waveletcoeff,Xlength,
		       FILTERLENGTH,SCALES,body_out);
#endif

#else

#if(DIM==3) 
 adjointwavelet_transform3D(waveletcoeff,Xlength,Ylength,Zlength,
			    FILTERLENGTH,SCALES,body_out);
#endif

#if(DIM==2) 
adjointwavelet_transform2D(waveletcoeff,Xlength,Ylength,
			   FILTERLENGTH,SCALES,body_out);
#endif

#if(DIM==1) 
adjointwavelet_transform1D(waveletcoeff,Xlength,
			   FILTERLENGTH,SCALES,body_out);
#endif

#endif
 time4=clock();
 // for(k=0;k<size;k++)printf("outdata  %f , ",body_out[k]);
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
