#include <stdio.h>
#include "bio.h"
#include "adjbio.h"



void getadjointfilter(p_bio_xd,p_bio_skip_xd,
	  p_bio_xd_char,p_bio_skip_xd_char,Filterlength,flag)

void (**p_bio_xd)();
void (**p_bio_skip_xd)();	  
void (**p_bio_xd_char)();
void (**p_bio_skip_xd_char)();
int flag;
int Filterlength;
{
 
if(flag==3){   /* 2d  */
switch(Filterlength){
 case 1:
*p_bio_xd=adjbioD1_2d;
*p_bio_skip_xd=adjbioD1_2d;	  
*p_bio_xd_char=adjbioD1_2d;
*p_bio_skip_xd_char=adjbioD1_2d;
 break;
 case 3:
*p_bio_xd=adjbioD3_2d;
*p_bio_skip_xd=adjbioD3_2d;	  
*p_bio_xd_char=adjbioD3_2d;
*p_bio_skip_xd_char=adjbioD3_2d;
 break;
 case 5:
*p_bio_xd=adjbioD5_2d;
*p_bio_skip_xd=adjbioD5_2d;	  
*p_bio_xd_char=adjbioD5_2d;
*p_bio_skip_xd_char=adjbioD5_2d;
 break;
 case 7:
*p_bio_xd=adjbioD7_2d;
*p_bio_skip_xd=adjbioD7_2d;	  
*p_bio_xd_char=adjbioD7_2d;
*p_bio_skip_xd_char=adjbioD7_2d;
 break;
 case 9:
*p_bio_xd=adjbioD9_2d;
*p_bio_skip_xd=adjbioD9_2d;	  
*p_bio_xd_char=adjbioD9_2d;
*p_bio_skip_xd_char=adjbioD9_2d;
 break;
 case -1:
*p_bio_xd=adjbioR1_2d;
*p_bio_skip_xd=adjbioR1_2d;	  
*p_bio_xd_char=adjbioR1_2d;
*p_bio_skip_xd_char=adjbioR1_2d;
 break;
 case -3:
*p_bio_xd=adjbioR3_2d;
*p_bio_skip_xd=adjbioR3_2d;	  
*p_bio_xd_char=adjbioR3_2d;
*p_bio_skip_xd_char=adjbioR3_2d;
 break;
 case -5:
*p_bio_xd=adjbioR5_2d;
*p_bio_skip_xd=adjbioR5_2d;	  
*p_bio_xd_char=adjbioR5_2d;
*p_bio_skip_xd_char=adjbioR5_2d;
 break;
 case -7:
*p_bio_xd=adjbioR7_2d;
*p_bio_skip_xd=adjbioR7_2d;	  
*p_bio_xd_char=adjbioR7_2d;
*p_bio_skip_xd_char=adjbioR7_2d;
 break;
 case -9:
*p_bio_xd=adjbioR9_2d;
*p_bio_skip_xd=adjbioR9_2d;	  
*p_bio_xd_char=adjbioR9_2d;
*p_bio_skip_xd_char=adjbioR9_2d;
  break;
 default:
   printf("Error: Non-valid filterlength\n");
   break;
}
}


 
if(flag==2) {  /* 2dy  */
switch(Filterlength){
 case 1:
*p_bio_xd=adjbioD1_2d_y;
*p_bio_skip_xd=adjbioD1_2d_y;	  
*p_bio_xd_char=adjbioD1_2d_y;
*p_bio_skip_xd_char=adjbioD1_2d_y;
 break;
 case 3:
*p_bio_xd=adjbioD3_2d_y;
*p_bio_skip_xd=adjbioD3_2d_y;	  
*p_bio_xd_char=adjbioD3_2d_y;
*p_bio_skip_xd_char=adjbioD3_2d_y;
 break;
 case 5:
*p_bio_xd=adjbioD5_2d_y;
*p_bio_skip_xd=adjbioD5_2d_y;	  
*p_bio_xd_char=adjbioD5_2d_y;
*p_bio_skip_xd_char=adjbioD5_2d_y;
 break;
 case 7:
*p_bio_xd=adjbioD7_2d_y;
*p_bio_skip_xd=adjbioD7_2d_y;	  
*p_bio_xd_char=adjbioD7_2d_y;
*p_bio_skip_xd_char=adjbioD7_2d_y;
 break;
 case 9:
*p_bio_xd=adjbioD9_2d_y;
*p_bio_skip_xd=adjbioD9_2d_y;	  
*p_bio_xd_char=adjbioD9_2d_y;
*p_bio_skip_xd_char=adjbioD9_2d_y;
 break;
 case -1:
*p_bio_xd=adjbioR1_2d_y;
*p_bio_skip_xd=adjbioR1_2d_y;	  
*p_bio_xd_char=adjbioR1_2d_y;
*p_bio_skip_xd_char=adjbioR1_2d_y;
 break;
 case -3:
*p_bio_xd=adjbioR3_2d_y;
*p_bio_skip_xd=adjbioR3_2d_y;	  
*p_bio_xd_char=adjbioR3_2d_y;
*p_bio_skip_xd_char=adjbioR3_2d_y;
 break;
 case -5:
*p_bio_xd=adjbioR5_2d_y;
*p_bio_skip_xd=adjbioR5_2d_y;	  
*p_bio_xd_char=adjbioR5_2d_y;
*p_bio_skip_xd_char=adjbioR5_2d_y;
 break;
 case -7:
*p_bio_xd=adjbioR7_2d_y;
*p_bio_skip_xd=adjbioR7_2d_y;	  
*p_bio_xd_char=adjbioR7_2d_y;
*p_bio_skip_xd_char=adjbioR7_2d_y;
 break;
 case -9:
*p_bio_xd=adjbioR9_2d_y;
*p_bio_skip_xd=adjbioR9_2d_y;	  
*p_bio_xd_char=adjbioR9_2d_y;
*p_bio_skip_xd_char=adjbioR9_2d_y;
  break;
 default:
   printf("Error: Non-valid filterlength\n");
   break;
}
}


if(flag==7)   /*3d */
   {

switch(Filterlength){
 case 1:
*p_bio_xd=adjbioD1_3d;
*p_bio_skip_xd=adjbioD1_3d;	  
*p_bio_xd_char=adjbioD1_3d;
*p_bio_skip_xd_char=adjbioD1_3d;
 break;
 case 3:
*p_bio_xd=adjbioD3_3d;
*p_bio_skip_xd=adjbioD3_3d;	  
*p_bio_xd_char=adjbioD3_3d;
*p_bio_skip_xd_char=adjbioD3_3d;
 break;
 case 5:
*p_bio_xd=adjbioD5_3d;
*p_bio_skip_xd=adjbioD5_3d;	  
*p_bio_xd_char=adjbioD5_3d;
*p_bio_skip_xd_char=adjbioD5_3d;
 break;
 case 7:
*p_bio_xd=adjbioD7_3d;
*p_bio_skip_xd=adjbioD7_3d;	  
*p_bio_xd_char=adjbioD7_3d;
*p_bio_skip_xd_char=adjbioD7_3d;
 break;
 case 9:
*p_bio_xd=adjbioD9_3d;
*p_bio_skip_xd=adjbioD9_3d;	  
*p_bio_xd_char=adjbioD9_3d;
*p_bio_skip_xd_char=adjbioD9_3d;
 break;
 case -1:
*p_bio_xd=adjbioR1_3d;
*p_bio_skip_xd=adjbioR1_3d;	  
*p_bio_xd_char=adjbioR1_3d;
*p_bio_skip_xd_char=adjbioR1_3d;
 break;
 case -3:
*p_bio_xd=adjbioR3_3d;
*p_bio_skip_xd=adjbioR3_3d;	  
*p_bio_xd_char=adjbioR3_3d;
*p_bio_skip_xd_char=adjbioR3_3d;
 break;
 case -5:
*p_bio_xd=adjbioR5_3d;
*p_bio_skip_xd=adjbioR5_3d;	  
*p_bio_xd_char=adjbioR5_3d;
*p_bio_skip_xd_char=adjbioR5_3d;
 break;
 case -7:
*p_bio_xd=adjbioR7_3d;
*p_bio_skip_xd=adjbioR7_3d;	  
*p_bio_xd_char=adjbioR7_3d;
*p_bio_skip_xd_char=adjbioR7_3d;
 break;
 case -9:
*p_bio_xd=adjbioR9_3d;
*p_bio_skip_xd=adjbioR9_3d;	  
*p_bio_xd_char=adjbioR9_3d;
*p_bio_skip_xd_char=adjbioR9_3d;
  break;
 default:
   printf("Error: Non-valid filterlength\n");
   break;
}
}


}













