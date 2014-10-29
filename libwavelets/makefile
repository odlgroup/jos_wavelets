
include  ../build/makefile.setting

RM	= rm -f
RM_CMD 	= $(RM)  *.o ,*~ *.a,


WL_SRCS=  	wavelet_transform.c \
		Decompose3.c  wavelet_dec3.c \
		bio13579.2d.c bio1.2d.c bio_d.c \
	  	bio13579.2dy.c bio1.2dy.c bio_dy.c \
		bio13579.3d.c bio1.3d.c bio_3d.c bio_3.c \
		wavelet_rec3.c   getfilter.c  


WL_OBJS=  	wavelet_transform.o \
		wavelet_dec3.o \
		bio13579.3d.o bio1.3d.o bio_3d.o \
		bio13579.2d.o bio1.2d.o bio_2d.o \
	 	bio13579.2dy.o bio1.2dy.o bio_2dy.o \
		wavelet_rec3.o   getfilter.o 



CC             = gcc
CCOPTS 	        = $(OFLAG)
INCLUDE	= 
LDOPTS	=	

library:  libwavelets.a


TARGET1=libwavelets.a



$(TARGET1):  $(WL_OBJS)
	ar -rcv $(TARGET1)   $(WL_OBJS)

.c.o:  
	$(CC) $(PRECISION)  $(CCOPTS) -c $< 



clean:
	$(RM_CMD) 

depend:
	makedepend -s "# DO NOT DELETE" -- $(SRC)


























