#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <math.h>


#define XLength 470
#define YLength 500
#define ZLength 400


 /* makefile for test  has an include file;  ../build/makefile.setting
     where "HIGH_PRECISION"   (is/is not)  defined  */


#ifndef HIGH_PRECISION
#define FLOAT float
#else 
#define FLOAT double
#endif

#include "wavelet_transform.h"
#include "test.h"


int main()
{
  int Xlength=XLength;
  int Ylength=YLength;
  int Zlength=ZLength;
  clock_t  time0, time1,time2,time3,time4,time5;
  int k,l,j;
 

    FLOAT *randbody,*randcoeff;
    //  FLOAT *waveletcoeff,*adjointvector;
  FLOAT *invwaveletvector,*adjointvector;
  FLOAT *coeff_save, *body_save;
  unsigned int a=0;
  unsigned int mask12=0xfff;  /*  mask with 12  "ones"  */     
  FLOAT * ptr;
int size;


 double   sum1,sum2,sum3;
  double xsum1,ysum1,zsum1;
 double xsum2,ysum2,zsum2;
 double xsum3,ysum3,zsum3;

 FLOAT * Dptr,*Wptr,*Cptr,*Aptr;


/* for 2D test and  1D test : */
//Zlength=1; 
/* for 1D test:  */
//Ylength=1;
size = Xlength*Ylength*Zlength;
 printf("Xlength=%d, Ylength=%d, Zlength=%d\n",Xlength,Ylength,Zlength);
 printf("size=%d,\n",size);
 //return 1;
/* Allocation of arrays used in test */

 randbody=(FLOAT*)calloc((unsigned int)size,sizeof(FLOAT));
 body_save=(FLOAT*)calloc((unsigned int)size,sizeof(FLOAT));
 randcoeff=(FLOAT*)calloc((unsigned int)size,sizeof(FLOAT)); 
 coeff_save=(FLOAT*)calloc((unsigned int)size,sizeof(FLOAT));

 //waveletcoeff=(FLOAT*)calloc((unsigned int)size,sizeof(FLOAT));
invwaveletvector=(FLOAT*)calloc((unsigned int)size,sizeof(FLOAT)); 
adjointvector=(FLOAT*)calloc((unsigned int)size,sizeof(FLOAT));

 time0=clock();    /* uses c-call clock()  for timing process */ 
 ptr=randbody;
for(int k=0;k<size;k++){
a=mask12&random();
//a=1.0;       /*generetin 12 bits unsigned integer */
 *ptr++ =(FLOAT)a;       /* converting and saving on FLOAT array  */
 }

 ptr=randcoeff;
for(int k=0;k<size;k++){
  a=mask12&random();
  //  a=1.0;
       /*generetin 12 bits unsigned integer */
   *ptr++ =(FLOAT)a;       /* converting and saving on FLOAT array  */
  //  *ptr++ =(FLOAT)1.0;
 }

time1=clock(); 
/*copying files */
 for(int k=0;k<size;k++)body_save[k]=randbody[k];
 for(int k=0;k<size;k++)coeff_save[k]=randcoeff[k];
 time2=clock();
 /*  calling the wavelet transform*/
 adjointinvwavelet_transform3D(randbody,Xlength,Ylength,Zlength,adjointvector);
 //  wavelet_transform3D(randbody,Xlength,Ylength,Zlength,waveletcoeff);
 // wavelet_transform2D(body,Xlength,Ylength,waveletcoeff);
 //wavelet_transform1D(body,Xlength,waveletcoeff);

 time3=clock();

 /* calling the inverse wavelettransform */
 /*  printf("\nrandcoeff:\n");
     for(v=0;v<Xlength*Ylength*Zlength;v++)printf(" %f ",randcoeff[v]);
*/

 //adjointwavelet_transform3D(randcoeff,Xlength,Ylength,Zlength,adjointvector);
 invwavelet_transform3D(randcoeff,Xlength,Ylength,Zlength,invwaveletvector );
 //invwavelet_transform2D(waveletcoeff,Xlength,Ylength,body_out);
 //invwavelet_transform1D(Waveletcoeff,Xlength,body_out);
 time4=clock();

 sum1=0.0;
 // 
 printf("\n");
 //for(k=0;k<size;k++)   sum1 +=body_save[k]*body_save[k];
 //for(k=0;k<size;k++){
     //sum1 +=waveletcoeff[k]*coeff_save[k];
 Dptr=body_save;
 Wptr=invwaveletvector;
 Cptr=coeff_save;
 Aptr=adjointvector;  
 zsum1=0.0;zsum2=0.0;zsum3=0.0;
 for(k=0;k<Zlength;k++){
   ysum1 = 0.0; ysum2 = 0.0; ysum3 = 0.0;
   for(j=0;j<Ylength;j++){
     xsum1=0.0; xsum2=0.0; xsum3=0.0;
     for(l=0;l<Xlength;l++){
       /*
       xsum1 += (*Wptr)*(*Cptr);
       xsum2 += (*Dptr)*(*Aptr);
       xsum3 += (*Wptr++)*(*Cptr++) - (*Dptr++)*(*Aptr++);
       */
       xsum1 += (*Dptr)*(*Dptr);
       xsum2 += (*Cptr)*(*Cptr);
       xsum3 += (*Wptr++)*(*Dptr++) - (*Cptr++)*(*Aptr++);

     }
   ysum1 += xsum1; ysum2 += xsum2; ysum3 += xsum3;
   }
   zsum1 += ysum1; zsum2 += ysum2;  zsum3 += ysum3;
 }
 sum1=zsum1;sum2=zsum2; sum3=zsum3;
printf("\n");

 time5=clock();

 printf("Inner produkt1=%10.4f, Inner produkt2 = %10.4f, DIFFERENS = %e\n", sum1,sum2,sum3); 

printf("SignaltoError rate %f decibel\n",
       (-10.0*log10((double)(fabs(sum3)/sqrt(sum1*sum2))

		    )));
/* Printing out timings of processes done above  */

 printf("Size %d x %d x %d,Gen Rand:%f sek,Cpy File  %fe Sek;  WT;  %f sek, \
 AT:  %f sek, In prod %f sek\n",
	Xlength,Ylength,Zlength, 0.000001*(time1-time0),0.000001*(time2-time1),0.000001*(time3-time2),0.000001*(time4-time3),0.000001*(time5-time4));


#if(0)
 sum1=0.0;
 // for(k=0;k<size;k++)sum1 +=waveletcoeff[k]*randcoeff[k];
 // for(k=0;k<size;k++)sum1 +=body_save[k]*body_save[k];
 sum2=0.0;
 // for(k=0;k<size;k++)sum2 +=waveletcoeff[k]*waveletcoeff[k];
 //for(k=0;k<size;k++)sum2 +=adjointvector[k]*adjointvector[k];
 time5=clock();

 // printf("Inner produkt3=%e, Inner produkt4 = %e, differens = %e\n",sum1,sum2,sum1-sum2); 
/* Printing out timings of processes done above  */

 /*
printf("Size %d x %d x %d,Gen Rand:%f sek,Cpy File  %fe Sek;  WT;  %f sek, \
 AT:  %f sek, In prod %f sek\n",
		Xlength,Ylength,Zlength, 0.000001*(time1-time0),0.000001*(time2-time1),0.000001*(time3-time2),0.000001*(time4-time3),0.000001*(time5-time4));
 */

#endif 
/*
 ptr=coeff_save;
 printf("\nCoefficientvector\n\n");
 for(k=0;k<Zlength;k++){
   for(j=0;j<Ylength;j++){
     for(m=0;m<Xlength;m++){
       printf(" %f ",*ptr++);
}   

     printf("\n");
   }
   printf("\n-------\n");
}
 ptr=adjointvector;
 printf("\nAdjointvector\n\n");
 for(k=0;k<Zlength;k++){
   for(j=0;j<Ylength;j++){
     for(m=0;m<Xlength;m++){
       printf(" %f ",*ptr++);
}   
     printf("\n");
}
   printf("\n-------\n");
}
 */


 if(randbody!=NULL)free(randbody);
 if(randcoeff!=NULL)free(randcoeff);
 if(adjointvector!=NULL)free(adjointvector);
 //if(waveletcoeff!=NULL)free(waveletcoeff);
if(invwaveletvector!=NULL)free(invwaveletvector);

 return 0;
}
