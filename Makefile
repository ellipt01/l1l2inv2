
include make.config

DESTDIR		= ./bin
DESTLIBDIR	= ./lib

LOCALLIBS	= -L./lib -ll1l2inv -lcdescent -lmgcal
LIBS		= $(BLAS_LIB) $(GSL_LIB) $(EXTRA_LIBS)
CPPFLAGS	= -I./include -I./mgcal/include -I./cdescent/include $(OPENMP_FLG)

LIBSRC_OBJS	= src/l1l2inv.o src/simeq.o src/smooth.o src/utils.o src/settings.o

L1L2INV_OBJS= src/main.o
LCVINTP_OBJS= src/lcurve_interp.o
OPTLAM_OBJS	= src/optimal_lambda.o

RECOV_OBJS	= tools/recover.o
EXTR_OBJS	= tools/extract.o
CROSS_OBJS	= tools/cross_sect.o
MAKEIN_OBJS	= demo/src/makeinput.o

OBJS		= $(LIBSRC_OBJS) $(L1L2INV_OBJS) $(LCV_OBJS) $(LCVINTP_OBJS) $(OPTLAM_OBJS)\
			  $(RECOV_OBJS) $(EXTR_OBJS) $(CROSS_OBJS) $(MAKEIN_OBJS)

SUBDIRS		= mgcal cdescent scripts xmat

PROGRAMS	= l1l2inv lcurve_interp optimal_lambda
TOOLS		= recover extract cross_sect
DEMO		= makeinput

all	:	libl1l2inv $(SUBDIRS) $(PROGRAMS) $(TOOLS) $(DEMO)

libl1l2inv:	$(LIBSRC_OBJS)
			$(AR) r $(DESTLIBDIR)/$@.a $(LIBSRC_OBJS)

# PROGRAMS
l1l2inv:	$(L1L2INV_OBJS) libl1l2inv
			$(CC) $(CFLAGS) -o src/$@ $(L1L2INV_OBJS) $(CPPFLAGS) $(LOCALLIBS) $(LIBS)

lcurve_interp:	$(LCVINTP_OBJS)
			$(CC) $(CFLAGS) -o src/$@ $(LCVINTP_OBJS) $(CPPFLAGS) $(LOCALLIBS) $(LIBS)

optimal_lambda:	$(OPTLAM_OBJS)
			$(CC) $(CFLAGS) -o src/$@ $(OPTLAM_OBJS) $(CPPFLAGS) $(LOCALLIBS) $(LIBS)

# TOOLS
recover:	$(RECOV_OBJS)
			$(CC) $(CFLAGS) -o tools/$@ $(RECOV_OBJS) $(CPPFLAGS) $(LOCALLIBS) $(LIBS)

extract:	$(EXTR_OBJS)
			$(CC) $(CFLAGS) -o tools/$@ $(EXTR_OBJS) $(CPPFLAGS) $(LOCALLIBS) $(LIBS)

cross_sect:	$(CROSS_OBJS)
			$(CC) $(CFLAGS) -o tools/$@ $(CROSS_OBJS) $(CPPFLAGS) $(LOCALLIBS) $(LIBS)

# DEMO
makeinput:	$(MAKEIN_OBJS)
			$(CC) $(CFLAGS) -o demo/src/$@ $(MAKEIN_OBJS) $(CPPFLAGS) $(LOCALLIBS) $(LIBS)


$(SUBDIRS):	FORCE
			$(MAKE) -C $@

FORCE:

.c.o:
			$(CC) $(CFLAGS) -o $*.o -c $(CPPFLAGS) $< 

install:
			@ for i in $(PROGRAMS) ; do \
				$(INSTALL) -m 755 src/$$i $(DESTDIR)/ ; \
			done
			@ for i in $(TOOLS) ; do \
				$(INSTALL) -m 755 tools/$$i $(DESTDIR)/ ; \
			done
			@ for i in $(DEMO) ; do \
				$(INSTALL) -m 755 demo/src/$$i demo/bin/ ; \
			done
			@ for i in $(SUBDIRS) ; do\
				$(MAKE) install -C $$i ; \
			done

clean-objs:
			$(RM) $(OBJS)
			@ for i in $(OBJS) ; do \
				$(RM) $$i ; \
			done
			@ for i in $(PROGRAMS) ; do \
				$(RM) src/$$i ; \
			done
			@ for i in $(TOOLS) ; do \
				$(RM) tools/$$i ; \
			done
			@ for i in $(DEMO) ; do \
				$(RM) demo/src/$$i ; \
			done
			@ for i in $(SUBDIRS) ; do \
				$(MAKE) clean-objs -C $$i ; \
			done
			$(RM) *~ */*~


clean:		clean-objs
			$(RM) $(DESTLIBDIR)/libl1l2inv.a
			@ for i in $(PROGRAMS) ; do \
				$(RM) $(DESTDIR)/$$i ; \
			done
			@ for i in $(TOOLS) ; do \
				$(RM) $(DESTDIR)/$$i ; \
			done
			@ for i in $(DEMO) ; do \
				$(RM) demo/$(DESTDIR)/$$i ; \
			done
			@ for i in $(SUBDIRS) ; do \
				echo $(MAKE) clean -C $$i ; \
				$(MAKE) clean -C $$i ; \
			done
			$(RM) *~ */*~

