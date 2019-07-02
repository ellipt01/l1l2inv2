CC		= gcc
RM		= rm -f
INSTALL	= install

DESTDIR	= ../bin

SCRIPTS	= calc.sh

all:

install:
	@ for i in $(SCRIPTS) ; do \
		$(INSTALL) -m 755 $$i $(DESTDIR)/ ; \
	done

clean-objs:
	$(RM) *~ */*~

clean:	clean-objs
	@ for i in $(SCRIPTS) ; do \
		$(RM) $(DESTDIR)/$$i ; \
	done


