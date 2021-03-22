# Linux makefile

LDLIBS=-lgdal
CPPFLAGS = -O2

default:	mbuilder

clean:
	$(RM) mbuilder
