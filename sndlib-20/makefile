SHELL = /bin/sh
top_srcdir = .

INSTALL = /usr/bin/install -c
prefix = /usr/local
exec_prefix = ${prefix}
datarootdir = @datarootdir@
bindir = ${exec_prefix}/bin
srcdir = .
mandir = ${prefix}/man
libdir = ${exec_prefix}/lib
includedir = ${prefix}/include

mkinstalldirs = $(SHELL) $(srcdir)/mkinstalldirs

CC = gcc
DEFS = -DHAVE_CONFIG_H
CFLAGS = -g -O2
LDFLAGS = 
LIBS = -lm 

GUILE_LIBS = @GUILE_LIBS@
GUILE_CFLAGS = @GUILE_CFLAGS@
AUDIO_LIB = @AUDIO_LIB@
GSL_LIBS = @GSL_LIBS@
GSL_FLAGS = @GSL_CFLAGS@

LDSO_FLAGS = @LDSO_FLAGS@
SO_FLAGS = @SO_FLAGS@
SO_INSTALL = @SO_INSTALL@
A_INSTALL = @A_INSTALL@
SO_LD = @SO_LD@
A_LD = @A_LD@
A_LD_FLAGS = @A_LD_FLAGS@
LD_FLAGS = @LD_FLAGS@

.c.o:
	$(CC) -c $(DEFS) $(CFLAGS) $(SO_FLAGS) $(GUILE_CFLAGS) $<

SNDLIB_HEADERS = sndlib.h sndlib-strings.h vct.h clm.h xen.h sndlib2xen.h clm2xen.h
SNDLIB_O_FILES = headers.o audio.o io.o sound.o xen.o vct.o clm.o sndlib2xen.o clm2xen.o midi.o
SNDLIB_SIMPLE_O_FILES = headers.o audio.o io.o sound.o


sndlib: $(SNDLIB_HEADERS) $(SNDLIB_O_FILES)
	$(SO_LD) $(LDFLAGS) $(SNDLIB_O_FILES) -o sndlib.so $(LDSO_FLAGS) $(AUDIO_LIB) $(GSL_LIBS) $(GUILE_LIBS) $(LIBS)
#	$(A_LD) $(LD_FLAGS) $(LDFLAGS) $(SNDLIB_O_FILES) -o sndlib.a $(A_LD_FLAGS) $(AUDIO_LIB) $(GSL_LIBS) $(GUILE_LIBS) $(LIBS)
	$(A_LD) $(A_LD_FLAGS) sndlib.a $(SNDLIB_O_FILES)

$(SNDLIB_O_FILES): $(SNDLIB_HEADERS)

clean:
	rm -f $(SNDLIB_O_FILES)

sndplay: $(SNDLIB_HEADERS) $(SNDLIB_O_FILES) sndplay.o
	$(CC) sndplay.c -o sndplay $(SNDLIB_SIMPLE_O_FILES) $(AUDIO_LIB) $(LIBS)

sndrecord: $(SNDLIB_HEADERS) $(SNDLIB_O_FILES) sndrecord.o
	$(CC) sndrecord.c -o sndrecord $(SNDLIB_SIMPLE_O_FILES) $(AUDIO_LIB) $(LIBS)

sndinfo: $(SNDLIB_HEADERS) $(SNDLIB_O_FILES) sndinfo.o
	$(CC) sndinfo.c -o sndinfo $(SNDLIB_SIMPLE_O_FILES) $(AUDIO_LIB) $(LIBS)

audinfo: $(SNDLIB_HEADERS) $(SNDLIB_O_FILES) audinfo.o
	$(CC) audinfo.c -o audinfo $(SNDLIB_SIMPLE_O_FILES) $(AUDIO_LIB) $(LIBS)

install: sndlib
	$(mkinstalldirs) $(bindir)
	$(mkinstalldirs) $(libdir)
	$(mkinstalldirs) $(includedir)
	$(SO_INSTALL) sndlib.so $(libdir)/libsndlib.so
	$(A_INSTALL) sndlib.a $(libdir)/libsndlib.a
	$(INSTALL) sndlib-config $(bindir)/sndlib-config
	$(INSTALL) sndlib.h $(includedir)/sndlib.h

uninstall:
	rm -f $(libdir)/libsndlib.so
	rm -f $(libdir)/libsndlib.a

Makefile: Makefile.in config.status
	./config.status

config.status: configure
	./config.status --recheck

configure: configure.in
	cd $(srcdir); autoconf

