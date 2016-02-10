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


#include <stdio.h>
#include "bio.h"



void getfilter(p_bio_xd,p_bio_skip_xd,
	  p_bio_xd_char,p_bio_skip_xd_char,Filterlength,flag)

void (**p_bio_xd)();
void (**p_bio_skip_xd)();	  
void (**p_bio_xd_char)();
void (**p_bio_skip_xd_char)();
int flag;
int Filterlength;
{
 


if(flag==7)   /*3d */
   {

switch(Filterlength){
 case 1:
*p_bio_xd=bioD1_3d;
*p_bio_skip_xd=bioD1_3d;	  
*p_bio_xd_char=bioD1_3d;
*p_bio_skip_xd_char=bioD1_3d;
 break;
 case 3:
*p_bio_xd=bioD3_3d;
*p_bio_skip_xd=bioD3_3d;	  
*p_bio_xd_char=bioD3_3d;
*p_bio_skip_xd_char=bioD3_3d;
 break;
 case 5:
*p_bio_xd=bioD5_3d;
*p_bio_skip_xd=bioD5_3d;	  
*p_bio_xd_char=bioD5_3d;
*p_bio_skip_xd_char=bioD5_3d;
 break;
 case 7:
*p_bio_xd=bioD7_3d;
*p_bio_skip_xd=bioD7_3d;	  
*p_bio_xd_char=bioD7_3d;
*p_bio_skip_xd_char=bioD7_3d;
 break;
 case 9:
*p_bio_xd=bioD9_3d;
*p_bio_skip_xd=bioD9_3d;	  
*p_bio_xd_char=bioD9_3d;
*p_bio_skip_xd_char=bioD9_3d;
 break;
 case -1:
*p_bio_xd=bioR1_3d;
*p_bio_skip_xd=bioR1_3d;	  
*p_bio_xd_char=bioR1_3d;
*p_bio_skip_xd_char=bioR1_3d;
 break;
 case -3:
*p_bio_xd=bioR3_3d;
*p_bio_skip_xd=bioR3_3d;	  
*p_bio_xd_char=bioR3_3d;
*p_bio_skip_xd_char=bioR3_3d;
 break;
 case -5:
*p_bio_xd=bioR5_3d;
*p_bio_skip_xd=bioR5_3d;	  
*p_bio_xd_char=bioR5_3d;
*p_bio_skip_xd_char=bioR5_3d;
 break;
 case -7:
*p_bio_xd=bioR7_3d;
*p_bio_skip_xd=bioR7_3d;	  
*p_bio_xd_char=bioR7_3d;
*p_bio_skip_xd_char=bioR7_3d;
 break;
 case -9:
*p_bio_xd=bioR9_3d;
*p_bio_skip_xd=bioR9_3d;	  
*p_bio_xd_char=bioR9_3d;
*p_bio_skip_xd_char=bioR9_3d;
  break;
 default:
   printf("Error: Non-valid filterlength\n");
   break;
}
}


}













