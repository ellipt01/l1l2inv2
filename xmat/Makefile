
include ../make.config

DESTDIR		= ../bin
DESTLIBDIR	= ./lib

LOCALLIBS	= -L./lib -ll1l2inv_xmat -lcdescent_xmat -L../lib -lmgcal 
LIBS		= $(BLAS_LIB) $(GSL_LIB) -lm $(OPENMP_FLG)
CPPFLAGS	= -I./include -I../include -I../mgcal/include -I./cdescent/include

LIBSRC_OBJS	= src/l1l2inv.o src/simeq.o ../src/smooth.o src/utils.o ../src/settings.o

L1L2INV_OBJS= src/main.o

RECOV_OBJS	= tools/recover.o
RESOL_OBJS	= tools/resolution.o

OBJS		= $(LIBSRC_OBJS) $(L1L2INV_OBJS) $(RECOV_OBJS) $(RESOL_OBJS)

SUBDIRS		= cdescent scripts

PROGRAMS	= l1l2inv_xmat
TOOLS		= recover_xmat resolution_xmat

all	:		libl1l2inv_xmat $(SUBDIRS) $(PROGRAMS) $(TOOLS)

libl1l2inv_xmat:	$(LIBSRC_OBJS)
			$(AR) r $(DESTLIBDIR)/$@.a $(LIBSRC_OBJS)

# PROGRAMS
l1l2inv_xmat:	$(L1L2INV_OBJS) libl1l2inv_xmat
			$(CC) $(CFLAGS) -o src/$@ $(L1L2INV_OBJS) $(CPPFLAGS) $(LOCALLIBS) $(LIBS)

# TOOLS
recover_xmat:	$(RECOV_OBJS)
			$(CC) $(CFLAGS) -o tools/$@ $(RECOV_OBJS) $(CPPFLAGS) $(LOCALLIBS) $(LIBS)

resolution_xmat:	$(RESOL_OBJS)
			$(CC) $(CFLAGS) -o tools/$@ $(RESOL_OBJS) $(CPPFLAGS) $(LOCALLIBS) $(LIBS)

$(SUBDIRS):	FORCE
			$(MAKE) -C $@

FORCE:

.c.o:
			$(CC) $(CFLAGS) -o $*.o -c $(CPPFLAGS) $< $(OPENMP_FLG)

install:
			@ for i in $(PROGRAMS) ; do \
				$(INSTALL) -m 755 src/$$i $(DESTDIR)/ ; \
			done
			@ for i in $(TOOLS) ; do \
				$(INSTALL) -m 755 tools/$$i $(DESTDIR)/ ; \
			done
			@ for i in $(SUBDIRS) ; do\
				$(MAKE) install -C $$i ; \
			done

clean-objs:
			$(RM) $(OBJS)
			@ for i in $(OBJS) ; do \
				$(RM) $$i ; \
			done
			@ for i in $(SUBDIRS) ; do \
				$(MAKE) clean-objs -C $$i ; \
			done
			$(RM) *~ */*~


clean:		clean-objs
			$(RM) $(DESTLIBDIR)/libl1l2inv_xmat.a
			@ for i in $(PROGRAMS) ; do \
				$(RM) src/$$i ; \
				$(RM) $(DESTDIR)/$$i ; \
			done
			@ for i in $(TOOLS) ; do \
				$(RM) tools/$$i ; \
				$(RM) $(DESTDIR)/$$i ; \
			done
			@ for i in $(SUBDIRS) ; do \
				echo $(MAKE) clean -C $$i ; \
				$(MAKE) clean -C $$i ; \
			done
			$(RM) *~ */*~

