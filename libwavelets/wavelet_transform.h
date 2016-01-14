
int wavelet_transform3D(FLOAT * ,//invector,
			int, //xlength,
			int ,//ylength,
			int ,//zlength,
                        int,//levels,
                        int, // Filterlength can be  1,3,5,7 or 9	
		FLOAT * //waveletcoefficients
			);

int invwavelet_transform3D(FLOAT * ,//waveltcoefficient,
			   int ,//xlength,
			   int ,//ylength,
			   int ,//zlength,
                        int, // Filterlength can be  1,3,5,7 or 9	
                        int,//levels,
			   FLOAT *//outvector
			   );

int wavelet_transform2D(FLOAT * ,//invector,
			int, //xlength,
			int ,//ylength,
                        int, // Filterlength can be  1,3,5,7 or 9	
                        int,//levels,
			FLOAT * //waveletcoefficients
			);

int invwavelet_transform2D(FLOAT * ,//waveltcoefficient,
			   int ,//xlength,
			   int ,//ylength,
                        int, // Filterlength can be  1,3,5,7 or 9	
                        int,//levels,
			   FLOAT *//outector
			   );

int wavelet_transform1D(FLOAT * ,//invector,
			int, //xlength,
                        int, // Filterlength can be  1,3,5,7 or 9	
                        int,//levels,
			FLOAT * //waveletcoefficients
			);

int invwavelet_transform1D(FLOAT * ,//waveltcoefficient,
			   int ,//xlength,
                        int, // Filterlength can be  1,3,5,7 or 9	
                        int,//levels,
			   FLOAT *//outvector
			   );

int adjointinvwavelet_transform3D(FLOAT * ,//invector,
			int, //xlength,
			int ,//ylength,
			int ,//zlength,
                        int, // Filterlength can be  1,3,5,7 or 9	
                        int,//levels,
			FLOAT * //waveletcoefficients
			);

int adjointwavelet_transform3D(FLOAT * ,//waveltcoefficient,
			   int ,//xlength,
			   int ,//ylength,
			   int ,//zlength,
                        int, // Filterlength can be  1,3,5,7 or 9	
                        int,//levels,
			   FLOAT *//outvector
			   );

int adjointinvwavelet_transform2D(FLOAT * ,//invector,
			int, //xlength,
			int ,//ylength,
                        int, // Filterlength can be  1,3,5,7 or 9	
                        int,//levels,
			FLOAT * //waveletcoefficients
			);

int adjointwavelet_transform2D(FLOAT * ,//waveltcoefficient,
			   int ,//xlength,
			   int ,//ylength,
                        int, // Filterlength can be  1,3,5,7 or 9	
                        int,//levels,
			   FLOAT *//outector
			   );

int adjointinvwavelet_transform1D(FLOAT * ,//invector,
			int, //xlength,
                        int, // Filterlength can be  1,3,5,7 or 9	
                        int,//levels,
			FLOAT * //waveletcoefficients
			);

int adjointwavelet_transform1D(FLOAT * ,//waveltcoefficient,
			   int ,//xlength,
                        int, // Filterlength can be  1,3,5,7 or 9	
                        int,//levels,
			   FLOAT *//outvector
			   );
