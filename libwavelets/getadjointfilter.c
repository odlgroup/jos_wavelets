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













