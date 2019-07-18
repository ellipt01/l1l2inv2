
include ../make.config

DESTLIBDIR	= ../lib
CFLAGS		= -O3
CPPFLAGS	= -I./include

LIBSRC_OBJS	= src/cdescent.o src/linregmodel.o src/regression.o src/update.o\
			  src/cyclic.o src/mmio.o src/stepsize.o\
			  src/io.o src/mmreal.o src/stochastic.o\
			  src/private/atomic.o src/private/private.o

all	:		libcdescent

libcdescent:	$(LIBSRC_OBJS)
			$(AR) r $(DESTLIBDIR)/$@.a $(LIBSRC_OBJS)

.c.o:
			$(CC) $(CFLAGS) -o $*.o -c $(CPPFLAGS) $< $(OPENMP_FLG)

install:

clean-objs:
			@ for i in $(LIBSRC_OBJS) ; do \
				$(RM) $$i ; \
			done
			$(RM) *~

clean:		clean-objs
			$(RM) $(DESTLIBDIR)/libcdescent.a

