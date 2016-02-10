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


int wavelet_transform3D(FLOAT * ,//invector,
			int, //xlength,
			int ,//ylength,
			int ,//zlength,
                        int,//levels,
                        char, // Filterlength can be  1,3,5,7 or 9	
		FLOAT * //waveletcoefficients
			);

int invwavelet_transform3D(FLOAT * ,//waveltcoefficient,
			   int ,//xlength,
			   int ,//ylength,
			   int ,//zlength,
                        char, // Filterlength can be  1,3,5,7 or 9	
                        int,//levels,
			   FLOAT *//outvector
			   );

int wavelet_transform2D(FLOAT * ,//invector,
			int, //xlength,
			int ,//ylength,
                        char, // Filterlength can be  1,3,5,7 or 9	
                        int,//levels,
			FLOAT * //waveletcoefficients
			);

int invwavelet_transform2D(FLOAT * ,//waveltcoefficient,
			   int ,//xlength,
			   int ,//ylength,
                        char, // Filterlength can be  1,3,5,7 or 9	
                        int,//levels,
			   FLOAT *//outector
			   );

int wavelet_transform1D(FLOAT * ,//invector,
			int, //xlength,
                        char, // Filterlength can be  1,3,5,7 or 9	
                        int,//levels,
			FLOAT * //waveletcoefficients
			);

int invwavelet_transform1D(FLOAT * ,//waveltcoefficient,
			   int ,//xlength,
                        char, // Filterlength can be  1,3,5,7 or 9	
                        int,//levels,
			   FLOAT *//outvector
			   );

int adjointinvwavelet_transform3D(FLOAT * ,//invector,
			int, //xlength,
			int ,//ylength,
			int ,//zlength,
                        char, // Filterlength can be  1,3,5,7 or 9	
                        int,//levels,
			FLOAT * //waveletcoefficients
			);

int adjointwavelet_transform3D(FLOAT * ,//waveltcoefficient,
			   int ,//xlength,
			   int ,//ylength,
			   int ,//zlength,
                        char, // Filterlength can be  1,3,5,7 or 9	
                        int,//levels,
			   FLOAT *//outvector
			   );

int adjointinvwavelet_transform2D(FLOAT * ,//invector,
			int, //xlength,
			int ,//ylength,
                        char, // Filterlength can be  1,3,5,7 or 9	
                        int,//levels,
			FLOAT * //waveletcoefficients
			);

int adjointwavelet_transform2D(FLOAT * ,//waveltcoefficient,
			   int ,//xlength,
			   int ,//ylength,
                        char, // Filterlength can be  1,3,5,7 or 9	
                        int,//levels,
			   FLOAT *//outector
			   );

int adjointinvwavelet_transform1D(FLOAT * ,//invector,
			int, //xlength,
                        char, // Filterlength can be  1,3,5,7 or 9	
                        int,//levels,
			FLOAT * //waveletcoefficients
			);

int adjointwavelet_transform1D(FLOAT * ,//waveltcoefficient,
			   int ,//xlength,
                        char, // Filterlength can be  1,3,5,7 or 9	
                        int,//levels,
			   FLOAT *//outvector
			   );
