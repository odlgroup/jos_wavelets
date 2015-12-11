
int wavelet_transform3D(FLOAT * ,//invector,
			int, //xlength,
			int ,//ylength,
			int ,//zlength,
			FLOAT * //waveletcoefficients
			);

int invwavelet_transform3D(FLOAT * ,//waveltcoefficient,
			   int ,//xlength,
			   int ,//ylength,
			   int ,//zlength,
			   FLOAT *//outvector
			   );

int wavelet_transform2D(FLOAT * ,//invector,
			int, //xlength,
			int ,//ylength,
			FLOAT * //waveletcoefficients
			);

int invwavelet_transform2D(FLOAT * ,//waveltcoefficient,
			   int ,//xlength,
			   int ,//ylength,
			   FLOAT *//outector
			   );

int wavelet_transform1D(FLOAT * ,//invector,
			int, //xlength,
			FLOAT * //waveletcoefficients
			);

int invwavelet_transform1D(FLOAT * ,//waveltcoefficient,
			   int ,//xlength,
			   FLOAT *//outvector
			   );

int adjointinvwavelet_transform3D(FLOAT * ,//invector,
			int, //xlength,
			int ,//ylength,
			int ,//zlength,
			FLOAT * //waveletcoefficients
			);

int adjointwavelet_transform3D(FLOAT * ,//waveltcoefficient,
			   int ,//xlength,
			   int ,//ylength,
			   int ,//zlength,
			   FLOAT *//outvector
			   );

int adjointinvwavelet_transform2D(FLOAT * ,//invector,
			int, //xlength,
			int ,//ylength,
			FLOAT * //waveletcoefficients
			);

int adjointwavelet_transform2D(FLOAT * ,//waveltcoefficient,
			   int ,//xlength,
			   int ,//ylength,
			   FLOAT *//outector
			   );

int adjointinvwavelet_transform1D(FLOAT * ,//invector,
			int, //xlength,
			FLOAT * //waveletcoefficients
			);

int adjointwavelet_transform1D(FLOAT * ,//waveltcoefficient,
			   int ,//xlength,
			   FLOAT *//outvector
			   );
