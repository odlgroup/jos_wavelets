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


/* file bio13579.3d.c */


/*       Fast algorithm for a symmetric biorthorgonal Low and High
          pass filter  of length  1,3,5,7 and 9
         i dimension two.
         This is a control file. The used procedures are defined
         in the files bio_3.c  (the local operations)
                      bio1.3d.c ordering in HH.LH.HL and LL subbands

                                 Algorithm and code developped 
                                               by
                                      Jan-Olov Stromberg

				      KTH, Stockholm, Sweden
                     Fast Mathematical Algorithms&Hardware, Hamden, CT, USA


                             Preliminary  version by October 18 , 1997

*/


#include "bio_parameters.h"
//#include "boundary_macro.h"    old stuff not use

#ifndef FLOAT
#define FLOAT float
#endif

#ifndef INTYPE
#define INTYPE FLOAT
#endif

#ifndef OUTTYPE
#define OUTTYPE FLOAT
#endif
#include "bio.h"



#ifndef NORMALIZATION9
#define NORMALIZATION9   (0.763)
#endif


#ifndef NORMALIZATION7
#define NORMALIZATION7   (0.98)
#endif


#ifndef NORMALIZATION5
#define NORMALIZATION5   (0.5)
#endif
#ifndef NORMALIZATION3
#define NORMALIZATION3   (1.0)
#endif

/**********************************************/

void bioD9_3d(invector, HHH,HHL,HLH,HLL,
	      LHH,LHL,LLH,LLL,xlength,ylength,zlength,
	      ifnotAllskip)
INTYPE *invector;
FLOAT *HLL,*HLH,*HHL,*HHH;
FLOAT *LLL,*LLH,*LHL,*LHH;
int  xlength,ylength,zlength; 
int ifnotAllskip;

{ 





	FLOAT x0=    -1.586134342059;
	FLOAT x1=    -0.052980118573;
	FLOAT x2=     0.882911075529;
	FLOAT x3=     0.443506852045;

	/* Rotlevels=4;*/

FLOAT X0,X1,X2,X3,X4;
int pairity; 
FLOAT Norm;
 int dim;

 dim=(xlength>1)+(ylength>1)+(zlength>1);

 if(0||ifnotAllskip){
 X0=x0;
 X1=x0*x1;
 X2=x1*x2;
 X3=x2*x3;
 X4=1.0/x3;
 }else{
 X0=x0;
 X1=x0*x1;
 X2=x1*x2;
 X3=1.0;
 
 X4=1.0/x2;
 }

 if(dim==3)Norm =(x0*x1*x2);
 if(dim==2)Norm=1.0; 
 if(dim==1)Norm= 1.0/(x0*x1*x2);


 if(xlength>1) Norm /= NORMALIZATION9;
if(ylength>1) Norm  /=  NORMALIZATION9; 
 if(zlength>1) Norm  /=  NORMALIZATION9;



if(dim==0){Norm=1.0;X0=1.0;X1=1.0;X2=1.0;X3=1.0;X4=1.0;}



/*warning invector will be overwritten */
pairity=1;

bio_3d_premult(invector,xlength,ylength,zlength,X0,pairity);
pairity=0;

bio_3d_premult(invector,xlength,ylength,zlength,X1,pairity);
pairity=1;

bio_3d_premult(invector,xlength,ylength,zlength,X2,pairity);
pairity=0;

bio_3d_premult(invector,xlength,ylength,zlength,X3,pairity);

if(ifnotAllskip){
bioD1__3d(invector,HHH,HHL,HLH,HLL,LHH,LHL,LLH,LLL,
			  xlength,ylength,zlength,X4,Norm);
}else  bioD1_3dskip(invector,LLL,xlength,ylength,zlength,X4,Norm);

      }

/*********************************************/
 void bioR9_3d
(HHH,HHL,HLH,HLL,LHH,LHL,LLH,LLL,
outvector,xlength,ylength,zlength,
 ifnotAllskip)
FLOAT *HLL,*HLH,*HHL,*HHH;
FLOAT *LLL,*LLH,*LHL,*LHH;
OUTTYPE *outvector;
int  xlength,ylength; 
int ifnotAllskip;
{
	FLOAT ix0=    1.586134342059;
	FLOAT ix1=     0.052980118573;
	FLOAT ix2=    -0.882911075529;
	FLOAT ix3=     -0.443506852045;



FLOAT iX0,iX1,iX2,iX3,iX4;
int pairity; 
FLOAT Norm;
int dim;

 dim=(xlength>1)+(ylength>1)+(zlength>1);

if(ifnotAllskip){
  iX4=ix3;
  iX3=1.0/(ix3*ix2);
  iX2=1.0/(ix2*ix1);
  iX1=1.0/(ix1*ix0);
  iX0=1.0/ix0;
}else{
 iX4=ix2;
  iX3=1.0;
  iX2=1.0/(ix2*ix1);
  iX1=1.0/(ix1*ix0);
  iX0=1.0/ix0;
}





 if(dim==3)Norm =1.0/(ix0*ix1*ix2);
 if(dim==2)Norm=1.0; 
 if(dim==1)Norm= (ix0*ix1*ix2);



 if(xlength>1) Norm *= NORMALIZATION9;
if(ylength>1) Norm  *=  NORMALIZATION9; 
 if(zlength>1) Norm  *=  NORMALIZATION9;

if(dim==0){Norm=1.0;iX0=1.0;iX1=1.0;iX2=1.0;iX3=1.0;iX4=1.0;}


 
if(ifnotAllskip){
  bioR1__3d( HHH,HHL,HLH,HLL,LHH,LHL,LLH,LLL,
	    outvector,xlength,ylength,zlength,iX4,Norm);
pairity=0;

bio_3d_postmult(outvector,xlength,ylength,zlength,iX3,pairity);
}

else  bioR1_3dskip(LLL,outvector,xlength,ylength,zlength,iX4,Norm);
pairity=1;
bio_3d_postmult(outvector,xlength,ylength,zlength,iX2,pairity);
pairity=0;
bio_3d_postmult(outvector,xlength,ylength,zlength,iX1,pairity);
pairity=1;
bio_3d_postmult(outvector,xlength,ylength,zlength,iX0,pairity);
   }

/*********************************************/
void bioD7_3d
           (invector, HHH,HHL,HLH,HLL,
	    LHH,LHL,LLH,LLL,xlength,ylength,zlength,
           ifnotAllskip)
INTYPE *invector;
FLOAT *HLL,*HLH,*HHL,*HHH;
FLOAT *LLL,*LLH,*LHL,*LHH;
int  xlength,ylength,zlength; 
int ifnotAllskip;
{ 


 FLOAT x0=    0.2;
 FLOAT x1=   -0.357142857136;
 FLOAT x2=     0.21;

FLOAT X0,X1,X2,X3;
int pairity; 
FLOAT Norm;
int dim;

 dim=(xlength>1)+(ylength>1)+(zlength>1);

 if(0||ifnotAllskip){
X0=x0;
X1=x0*x1;
X2=x1*x2;
X3=1.0/x2;
 }else{
X0=x0;
X1=x0*x1;
X2=1;
X3=1.0/x1;

}


 if(dim==3)Norm =(x0*x1);
 if(dim==2)Norm=1.0; 
 if(dim==1)Norm= 1.0/(x0*x1);



if(xlength>1) Norm /= NORMALIZATION7;
if(ylength>1) Norm  /=  NORMALIZATION7; 
if(zlength>1)Norm  /=  NORMALIZATION7;


if(dim==0){Norm=1.0;X0=1.0;X1=1.0;X2=1.0;X3=1.0;}





/*warning invector will be overwritten */

pairity=0;
bio_3d_premult(invector,xlength,ylength,zlength,X0,pairity);
pairity=1;
bio_3d_premult(invector,xlength,ylength,zlength,X1,pairity);
pairity=0;
bio_3d_premult(invector,xlength,ylength,zlength,X2,pairity);


if(ifnotAllskip)bioD1__3d(invector,HHH,HHL,HLH,HLL,LHH,LHL,LLH,LLL,
			  xlength,ylength,zlength,X3,Norm);
else  bioD1_3dskip(invector,LLL,xlength,ylength,zlength,X3,Norm);
      
}
/**********************************************/
 void bioR7_3d(HHH,HHL,HLH,HLL,LHH,LHL,LLH,LLL,
	       outvector,xlength,ylength,zlength,
 ifnotAllskip)
FLOAT *HLL,*HLH,*HHL,*HHH;
FLOAT *LLL,*LLH,*LHL,*LHH;
OUTTYPE *outvector;
int  xlength,ylength; 
int ifnotAllskip;
{


 FLOAT ix0=    -0.2;
 FLOAT ix1=     0.357142857136;
 FLOAT ix2=    -0.21;

FLOAT iX0,iX1,iX2,iX3;
int pairity; 
FLOAT Norm;
int dim;

 dim=(xlength>1)+(ylength>1)+(zlength>1);

 if(ifnotAllskip){
  iX3=ix2;
  iX2=1.0/(ix2*ix1);
  iX1=1.0/(ix1*ix0);
  iX0=1.0/ix0;
 }else{
  iX3=ix1;
  iX2=1.0;
  iX1=1.0/(ix1*ix0);
  iX0=1.0/ix0;
 }


 /*
Norm = 1.0 /(iX0*iX1*iX2*iX3) ;
if(xlength>1) Norm *= iX1;
if(ylength>1) Norm *= iX1;
if(zlength>1) Norm *= iX1;
 */

 if(dim==3)Norm =1.0/(ix0*ix1);
 if(dim==2)Norm=1.0; 
 if(dim==1)Norm= (ix0*ix1);



 if(xlength>1) Norm *= NORMALIZATION7;
if(ylength>1) Norm  *=  NORMALIZATION7; 
 if(zlength>1) Norm  *=  NORMALIZATION7;



if(dim==0){Norm=1.0;iX0=1.0;iX1=1.0;iX2=1.0;iX3=1.0;}



 
 

if(ifnotAllskip){
bioR1__3d(HHH,HHL,HLH,HLL,LHH,LHL,LLH,LLL,
	  outvector,xlength,ylength,zlength,iX3,Norm);
pairity=0;

bio_3d_postmult(outvector,xlength,ylength,zlength,iX2,pairity);
}
else  bioR1_3dskip(LLL,outvector,xlength,ylength,zlength,iX3,Norm);
pairity=1;
bio_3d_postmult(outvector,xlength,ylength,zlength,iX1,pairity);
pairity=0;
bio_3d_postmult(outvector,xlength,ylength,zlength,iX0,pairity);
   }


/*****************************************************/
void bioD5_3d
           (invector, HHH,HHL,HLH,HLL,
	    LHH,LHL,LLH,LLL,xlength,ylength,zlength,
           ifnotAllskip)
INTYPE *invector;
FLOAT *HLL,*HLH,*HHL,*HHH;
FLOAT *LLL,*LLH,*LHL,*LHH;
int  xlength,ylength,zlength; 
int ifnotAllskip;

{ 

  FLOAT x0=    -0.50;
  FLOAT x1=	0.25;

FLOAT X0,X1,X2;
int pairity; 
FLOAT Norm;
int dim;

 dim=(xlength>1)+(ylength>1)+(zlength>1);

  if(0||ifnotAllskip){
 X0=x0;
X1=x0*x1;
 X2=1.0/x1;
  }else{
 X0=x0;
 X1=1.0;
 X2=1.0/x0;
  }


 if(dim==3)Norm =x0;
 if(dim==2)Norm=1.0; 
 if(dim==1)Norm= 1.0/x0;


 if(xlength>1) Norm /= NORMALIZATION5;
if(ylength>1) Norm  /=  NORMALIZATION5; 
 if(zlength>1) Norm  /=  NORMALIZATION5;


if(dim==0){Norm=1.0;X0=1.0;X1=1.0;X2=1.0;}



/*warning invector will be overwritten */
pairity=1;
bio_3d_premult(invector,xlength,ylength,zlength,X0,pairity);
pairity=0;
bio_3d_premult(invector,xlength,ylength,zlength,X1,pairity);
if(ifnotAllskip)bioD1__3d(invector,HHH,HHL,HLH,HLL,LHH,LHL,LLH,LLL,
			  xlength,ylength,zlength,X2,Norm);
else  {
bioD1_3dskip(invector,LLL,xlength,ylength,zlength,X2,Norm);
}

      }
/******************************************************/
 void bioR5_3d
(HHH,HHL,HLH,HLL,LHH,LHL,LLH,LLL,
outvector,xlength,ylength,zlength,
 ifnotAllskip)
FLOAT *HLL,*HLH,*HHL,*HHH;
FLOAT *LLL,*LLH,*LHL,*LHH;
OUTTYPE *outvector;
int  xlength,ylength; 
int ifnotAllskip;
{

  FLOAT ix0=      0.50;
  FLOAT ix1=     -0.25;


FLOAT iX0,iX1,iX2;
int pairity;
FLOAT Norm;
int dim;

dim=(xlength>1)+(ylength>1)+(zlength>1);

  if(ifnotAllskip){
 iX2=ix1;
 iX1=1.0/(ix1*ix0);
 iX0=1.0/ix0;
  }else{
 iX2=ix0;
 iX1=1.0;
 iX0=1.0/ix0;
  }



 if(dim==3)Norm =1.0/ix0;
 if(dim==2)Norm=1.0; 
 if(dim==1)Norm= ix0;


if(xlength>1) Norm *= NORMALIZATION5;
if(ylength>1) Norm  *=  NORMALIZATION5; 
 if(zlength>1) Norm  *=  NORMALIZATION5;

if(dim==0){Norm=1.0;iX0=1.0;iX1=1.0;iX2=1.0;}



if(ifnotAllskip){
bioR1__3d(HHH,HHL,HLH,HLL,LHH,LHL,LLH,LLL,
	  outvector,xlength,ylength,zlength,iX2,Norm);
pairity=0;
bio_3d_postmult(outvector,xlength,ylength,zlength,iX1,pairity);
}
else  bioR1_3dskip(LLL,outvector,xlength,ylength,zlength,iX2,Norm);
pairity=1;
bio_3d_postmult(outvector,xlength,ylength,zlength,iX0,pairity);
   }


/********************************/
void bioD3_3d
           (invector, HHH,HHL,HLH,HLL,
	    LHH,LHL,LLH,LLL,xlength,ylength,zlength,
           ifnotAllskip)
INTYPE *invector;
FLOAT *HLL,*HLH,*HHL,*HHH;
FLOAT *LLL,*LLH,*LHL,*LHH;
int  xlength,ylength,zlength; 
int ifnotAllskip;
{ 
  FLOAT x0=    -0.5;

 

FLOAT X0=x0; 
FLOAT X1=x0;


int pairity; 

FLOAT Norm;
 int dim;






dim = (xlength > 1)+(ylength > 1)+(zlength > 1);
;
if(dim==3)Norm=x0;
if(dim==2)Norm=1.0;
if(dim==1)Norm=1.0/x0;


 if(xlength>1) Norm /= NORMALIZATION3;
if(ylength>1) Norm  /=  NORMALIZATION3; 
 if(zlength>1) Norm  /=  NORMALIZATION3;

 if(dim==0){Norm=1.0;X0=1.0,X1=1.0;}



/*warning invector will be overwritten */
 pairity=1; /* NOTE: A DIFFERENT PAIRITY HERE THE ON LONGER FILTERS */ 
  bio_3d_premult(invector,xlength,ylength,zlength,X0,pairity);
  /*bio_3d(invector,xlength,ylength,zlength,tempvector,x0,pairity);*/
  /* switching needed since we have pairith=1: */

if(ifnotAllskip)bioD1__3d(invector,HHH,HHL,HLH,HLL,LHH,LHL,LLH,LLL,
			  xlength,ylength,zlength,X1,Norm);
else  bioD1_3dskip(invector,LLL,xlength,ylength,zlength,X1,Norm);

      }
/**********************************************/
 void bioR3_3d
(HHH,HHL,HLH,HLL,LHH,LHL,LLH,LLL,
outvector,xlength,ylength,zlength,
 ifnotAllskip)
FLOAT *HLL,*HLH,*HHL,*HHH;
FLOAT *LLL,*LLH,*LHL,*LHH;
OUTTYPE *outvector;
int  xlength,ylength; 
int ifnotAllskip;
{
 FLOAT ix0=    0.5;
 
 


 FLOAT  iX1;
 FLOAT iX0;

int pairity;
FLOAT Norm;
 int dim;

dim = (xlength > 1)+(ylength > 1)+(zlength > 1);




  if(ifnotAllskip){
 iX1=ix0;
 iX0=1.0/ix0;
  }else{
 iX1=1.0;
 iX0=1.0;
  }

  iX1=1/iX1;

if(dim==3)Norm=1.0/ix0;
if(dim==2)Norm=1.0;
if(dim==1)Norm=ix0;


 if(xlength>1) Norm *= NORMALIZATION3;
if(ylength>1) Norm  *=  NORMALIZATION3; 
 if(zlength>1) Norm  *=  NORMALIZATION3;



 if(dim==0){Norm=1.0;iX0=1.0;iX1=1.0;}
 
if(ifnotAllskip){
bioR1__3d(HHH,HHL,HLH,HLL,LHH,LHL,LLH,LLL,
	  outvector,xlength,ylength,zlength,iX1,Norm);
 pairity=1; /* NOTE: A DIFFERENT PAIRITY HERE THE ON LONGER FILTERS */ 
bio_3d_postmult(outvector,xlength,ylength,zlength,iX0,pairity);
}else  bioR1_3dskip(LLL,outvector,xlength,ylength,zlength,iX1,Norm);
   }

/***************************************/ 
void bioD1_3d
           (invector, HHH,HHL,HLH,HLL,
	    LHH,LHL,LLH,LLL,xlength,ylength,zlength,
           ifnotAllskip)
INTYPE *invector;
FLOAT *HLL,*HLH,*HHL,*HHH;
FLOAT *LLL,*LLH,*LHL,*LHH;
int  xlength,ylength,zlength; 
int ifnotAllskip;
{ 
FLOAT  X0=1.0; 
FLOAT Norm=1.0;

if(ifnotAllskip)bioD1__3d(invector,HHH,HHL,HLH,HLL,LHH,LHL,LLH,LLL,
			  xlength,ylength,zlength,X0,Norm);
else  bioD1_3dskip(invector,LLL,xlength,ylength,zlength,X0,Norm);
      }
/*************************************************/
 void bioR1_3d
(HHH,HHL,HLH,HLL,LHH,LHL,LLH,LLL,
outvector,xlength,ylength,zlength,
 ifnotAllskip)
FLOAT *HLL,*HLH,*HHL,*HHH;
FLOAT *LLL,*LLH,*LHL,*LHH;
OUTTYPE *outvector;
int  xlength,ylength; 
int ifnotAllskip;
{
FLOAT iX0 =1.0;

FLOAT Norm=1.0;


 
if(ifnotAllskip)bioR1__3d(HHH,HHL,HLH,HLL,LHH,LHL,LLH,LLL,
			  outvector,xlength,ylength,zlength,iX0,Norm);
else  bioR1_3dskip(LLL,outvector,xlength,ylength,zlength,iX0,Norm);
         }
































































































































































































