set(Headers
	bio.h
	bio_parameters.h
	wavelet_transform.h
 	adjbio.h
        )

set(Sources
	bio_3d.c
	bio1.3d.c
	bio13579.3d.c
	getfilter.c
	wavelet_dec3.c
	wavelet_rec3.c
	wavelet_transform.c
        adjointwavelet_transform.c
        wavelet_adj3.c
	wavelet_invadj3.c
 	getadjointfilter.c
	adjbio13579.3d.c
	bpairitymult.c
	)


#include_directories ("${PROJECT_SOURCE_DIR}/DEPENDECY")

add_library(libwavelets ${Headers} ${Sources})
#target_link_libraries(libwavelets DEPENDECY)
