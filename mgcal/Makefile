CC			= gcc
AR			= ar
RM			= rm -f
MAKE		= make
INSTALL		= install

OPENMP_FLG	= -fopenmp

DESTDIR		= ../lib
LIBS		= -lm $(OPENMP_FLG)
CFLAGS		= -O3
CPPFLAGS	= -I./include

LIBSRC_OBJS	= src/calc.o src/grid.o src/kernel.o src/scattered.o src/vector3d.o\
			  src/data_array.o src/io.o src/mgcal.o src/source.o src/private/util.o\
			  src/calc_bhattacharya.o

all	:		libmgcal

libmgcal:	$(LIBSRC_OBJS)
			$(AR) r $(DESTDIR)/$@.a $(LIBSRC_OBJS)

.c.o:
			$(CC) $(CFLAGS) -o $*.o -c $(CPPFLAGS) $< $(OPENMP_FLG)

install:

clean-objs:
			@ for i in $(LIBSRC_OBJS) ; do \
				$(RM) $$i ; \
			done
			$(RM) *~ */*~

clean:		clean-objs
			$(RM) $(DESTDIR)/libmgcal.a

