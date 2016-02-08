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


#define DIM   3           /*  dimension 1, 2 or 3 */
#define FILTERLENGTH    3     /* 1, 3, 5, 7 or 9 */
/* FILTERLENGtH 1: decimation only  L is even decimation; 
                                    H is odd decimation */  
#define SCALES 2

#define XLength 400
#define YLength 100    /*   used only in dimension 2 and 3 */
#define ZLength 400      /*  used only in dimension 3  */



/***********DEFINITINS ABOVE MAY BE CHANGES BY USER     *****************/

 /* makefile for test  has an include file;  ../build/makefile.setting
     where "HIGH_PRECISION"   (is/is not)  defined  */



#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <math.h>



#ifndef HIGH_PRECISION
#define FLOAT float
#else 
#define FLOAT double
#endif

#include "../libwavelets/wavelet_transform.h"
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



  clock_t  time0, time1,time2,time3,time4,time5,time6,time7;
  clock_t  time3a, time5a;
  int k,l,j;
 

    FLOAT *randbody,*randcoeff;
    FLOAT *waveletcoeff;
    FLOAT *invwaveletvector,*adjointvector,*invadjointvector;
  FLOAT *coeff_save, *body_save;
  unsigned int a=0;
  unsigned int mask12=0xfff;  /*  mask with 12  "ones"  */     
  FLOAT * ptr;
int size;


 double   sum1,sum2,sum3,sum4;
  double xsum1,ysum1,zsum1;
 double xsum2,ysum2,zsum2;
 double xsum3,ysum3,zsum3;
 double xsum4,ysum4,zsum4;

 FLOAT * Dptr,*iWptr,*Wptr, *Cptr,*Aptr,*iAptr;;


 printf("Test of adjoint with symmetric biorthogonal filters \n \
         Filterlength=%d, Number of scales=%d\n\n",FILTERLENGTH,SCALES);


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

 waveletcoeff=(FLOAT*)calloc((unsigned int)size,sizeof(FLOAT));
invwaveletvector=(FLOAT*)calloc((unsigned int)size,sizeof(FLOAT)); 
adjointvector=(FLOAT*)calloc((unsigned int)size,sizeof(FLOAT));
invadjointvector=(FLOAT*)calloc((unsigned int)size,sizeof(FLOAT));

 time0=clock();    /* uses c-call clock()  for timing process */ 
 ptr=randbody;
for(k=0;k<size;k++){
a=mask12&random();
//a=1.0;       /*generetin 12 bits unsigned integer */
 *ptr++ =(FLOAT)a;       /* converting and saving on FLOAT array  */
 }

 ptr=randcoeff;
for(k=0;k<size;k++){
  a=mask12&random();
  //  a=1.0;
       /*generetin 12 bits unsigned integer */
   *ptr++ =(FLOAT)a;       /* converting and saving on FLOAT array  */
  //  *ptr++ =(FLOAT)1.0;
 }

time1=clock(); 
/*copying files */
 for(k=0;k<size;k++)body_save[k]=randbody[k];
 for(k=0;k<size;k++)coeff_save[k]=randcoeff[k];
 time2=clock();
 /*  calling the wavelet transform*/
#if(DIM==3)
 adjointinvwavelet_transform3D(randbody,Xlength,Ylength,Zlength,FILTERLENGTH,
			       SCALES,invadjointvector);
 #endif
#if(DIM==2)
 adjointinvwavelet_transform2D(randbody,Xlength,Ylength,FILTERLENGTH,
			       SCALES,invadjointvector);
 #endif
#if(DIM==1)
 adjointinvwavelet_transform1D(randbody,Xlength,FILTERLENGTH,
			       SCALES,invadjointvector);
 #endif
 time3=clock();
 for(k=0;k<size;k++)randbody[k]=body_save[k];
 time3a=clock();
 #if(DIM==3)
 wavelet_transform3D(randbody,Xlength,Ylength,Zlength,FILTERLENGTH,
		    SCALES,waveletcoeff);
#endif
 #if(DIM==2)
 wavelet_transform2D(randbody,Xlength,Ylength,FILTERLENGTH,
		    SCALES,waveletcoeff);
#endif
 #if(DIM==1)
 wavelet_transform1D(randbody,Xlength,FILTERLENGTH,
		    SCALES,waveletcoeff);
#endif

time4=clock(); 

#if(DIM==3)
 adjointwavelet_transform3D(randcoeff,Xlength,Ylength,Zlength,FILTERLENGTH,
                           SCALES,adjointvector);
#endif
#if(DIM==2)
 adjointwavelet_transform2D(randcoeff,Xlength,Ylength,FILTERLENGTH,
                           SCALES,adjointvector);
#endif
#if(DIM==1)
 adjointwavelet_transform1D(randcoeff,Xlength,FILTERLENGTH,
                           SCALES,adjointvector);
#endif

time5=clock();
 for(k=0;k<size;k++)randcoeff[k]=coeff_save[k];
 time5a=clock();
#if(DIM==3)
invwavelet_transform3D(randcoeff,Xlength,Ylength,Zlength,FILTERLENGTH,
		       SCALES,invwaveletvector );
 #endif
#if(DIM==2)
invwavelet_transform2D(randcoeff,Xlength,Ylength,FILTERLENGTH,
		       SCALES,invwaveletvector );
 #endif
#if(DIM==1)
invwavelet_transform1D(randcoeff,Xlength,FILTERLENGTH,
		       SCALES,invwaveletvector );
 #endif
 time6=clock();

 sum1=0.0;
 // 
 printf("\n");

 Dptr=body_save;
 Cptr=coeff_save;
 Aptr=adjointvector;  
 iWptr=invwaveletvector;
 Wptr=waveletcoeff;
 iAptr=invadjointvector;  


 zsum1=0.0;zsum2=0.0;zsum3=0.0;zsum4=0.0;
 for(k=0;k<Zlength;k++){ysum1 = 0.0; ysum2 = 0.0; ysum3 = 0.0;ysum4=0.0;
   for(j=0;j<Ylength;j++){
     xsum1=0.0; xsum2=0.0; xsum3=0.0;xsum4=0.0;
     for(l=0;l<Xlength;l++){
       
       
       xsum1 += (*Dptr)*(*Dptr);
       xsum2 += (*Cptr)*(*Cptr);
       xsum3 += (*Aptr++)*(*Dptr) - (*Cptr)*(*Wptr++);
       xsum4 += (*iAptr++)*(*Cptr++) - (*Dptr++)*(*iWptr++);

     }
     ysum1 += xsum1; ysum2 += xsum2; ysum3 += xsum3; ysum4 +=xsum4;
   }
   zsum1 += ysum1; zsum2 += ysum2;  zsum3 += ysum3; zsum4 +=ysum4;
 }
 sum1=zsum1;sum2=zsum2; sum3=zsum3; sum4=zsum4;
printf("\n");

 time7=clock();

 printf("Inner produkt1=%e, Inner produkt2 = %e, adjointdiff = %e,invadjointdiff = %e\n",
 sum1,sum2,sum3,sum4); 

printf("SignaltoError rate adjoint  %f decibel\n",
       (-10.0*log10((double)(fabs(sum3)/sqrt(sum1*sum2))  )));
printf("SignaltoError rate invadjoint  %f decibel\n",
       (-10.0*log10((double)(fabs(sum4)/sqrt(sum1*sum2) ))));

/* Printing out timings of processes done above  */

 printf("Size %d x %d x %d,Gen Rand:%f sek,Cpy File  %fe Sek;\n  \
  iAT  %f sek, WT;  %f sek, AT;  %f sek,  iWT;  %f sek, \n \
, Inner prod %f sek\n",
	Xlength,Ylength,Zlength, 0.000001*(time1-time0),
	0.000001*((time2-time1)+(time3a-time3)+(time5a-time5)),
	0.000001*(time3-time2),0.000001*(time4-time3a),
0.000001*(time5-time4),0.000001*(time6-time5a),0.000001*(time7-time6));


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
 if(waveletcoeff!=NULL)free(waveletcoeff);
if(invwaveletvector!=NULL)free(invwaveletvector);

 return 0;
}
