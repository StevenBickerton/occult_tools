# Makefile for cdiffracsim tools
#    make -W detect.c detect DEBUG=1
#

ifdef PREC
	CFLAGS_PREC = -DPRECISION
endif

# if we're on my Mac, use gcc 4.2
CC = gcc
UNAME := $(shell uname -s)
ifeq ($(UNAME),Darwin)
	CC = clang
endif

HEADS = detection.h statistics.h fft-gsl.h fft-fftw.h
#OBJS  = libdetection.o libstatistics.o libfft-gsl.o libfft-fftw.o
OBJS  = libdetection.o libstatistics.o libfft-fftw.o

CFLAGS   = -g -Wall -O2 -std=c99 $(CFLAGS_PREC)
#CPPFLAGS = -I$(HOME)/usr/include -I/opt/local/include
CPPFLAGS = -I/opt/local/include
#LDFLAGS  = $(OBJS) -L$(HOME)/usr/lib -L/opt/local/lib \
#	-lcfitsio -lfftw3 -lfftw3f -lm -lgsl -lgslcblas

#LDFLAGS  = $(OBJS) -L$(HOME)/usr/lib -L/opt/local/lib 
LDFLAGS  = $(OBJS) -L/opt/local/lib \
	-lcfitsio -lfftw3 -lm -lgsl -lgslcblas

INSTALL = /usr/bin/install
BINDIR  = ${HOME}/usr/bin


########################################################
BINARIES = detect addKBO komplete makeKernel offsetPattern xcorrelate \
	xchi basis fresnelBox fresnelT hideKBO buildParams fresnelTA \
	hideKBO_TA fits2ascii

all: $(OBJS) $(BINARIES)

########################################################
# link the final executibles

detect: detect.o $(HEADS) $(OBJS)
	$(CC) $(CFLAGS) $(CPPFLAGS) $< $(LDFLAGS) -o $@
komplete: komplete.o $(HEADS) $(OBJS)
	$(CC) $(CFLAGS) $(CPPFLAGS) $< $(LDFLAGS) -o $@
basis: basis.o $(HEADS) $(OBJS)
	$(CC) $(CFLAGS) $(CPPFLAGS) $< $(LDFLAGS) -o $@
addKBO: addKBO.o $(HEADS) $(OBJS)
	$(CC) $(CFLAGS) $(CPPFLAGS) $< $(LDFLAGS) -o $@
makeKernel: makeKernel.o $(HEADS) $(OBJS)
	$(CC) $(CFLAGS) $(CPPFLAGS) $< $(LDFLAGS) -o $@
offsetPattern: offsetPattern.o $(HEADS) $(OBJS)
	$(CC) $(CFLAGS) $(CPPFLAGS) $< $(LDFLAGS) -o $@
xcorrelate: xcorrelate.o $(HEADS) $(OBJS)
	$(CC) $(CFLAGS) $(CPPFLAGS) $< $(LDFLAGS) -o $@
xchi: xchi.o $(HEADS) $(OBJS)
	$(CC) $(CFLAGS) $(CPPFLAGS) $< $(LDFLAGS) -o $@
hideKBO: hideKBO.o $(HEADS) $(OBJS)
	$(CC) $(CFLAGS) $(CPPFLAGS) $< $(LDFLAGS) -o $@
hideKBO_TA: hideKBO_TA.o $(HEADS) $(OBJS)
	$(CC) $(CFLAGS) $(CPPFLAGS) $< $(LDFLAGS) -o $@
buildParams: buildParams.o $(HEADS) $(OBJS)
	$(CC) $(CFLAGS) $(CPPFLAGS) $< $(LDFLAGS) -o $@
fits2ascii: fits2ascii.o $(HEADS) $(OBJS)
	$(CC) $(CFLAGS) $(CPPFLAGS) $< $(LDFLAGS) -o $@

build the objects:
%.o : %.c $(HEADS)
	$(CC) $(CFLAGS) $(CPPFLAGS) -c $<


#########################################################
#fresnelBox:  fresnelBox.c detection.h libdetection.o
#	$(CC) $(CFLAGS) $(CPPFLAGS) $< -o $@
#fresnelT:  fresnelT.c detection.h libdetection.o
#	$(CC) $(CFLAGS) $(CPPFLAGS) $< -o $@
#fresnelTA:  fresnelTA.c detection.h libdetection.o
#	$(CC) $(CFLAGS) $(CPPFLAGS) $< -o $@




#########################################################
.PHONY: clean

clean: 
	$(RM) *.o $(BINARIES)
	rm -rf test
	rm -rf *.dSYM

#########################################################
.PHONY: install

install:
	$(INSTALL) $(BINARIES)        $(BINDIR)/
	$(INSTALL) fresplot.pl        $(BINDIR)/fresplot
	$(INSTALL) plotBoxes.pl       $(BINDIR)/plotBoxes
	$(INSTALL) animFresPatt.pl    $(BINDIR)/animFresPatt
	$(INSTALL) plotStars.pl       $(BINDIR)/plotStars
	$(INSTALL) hitSummary.pl      $(BINDIR)/hitSummary	
	$(INSTALL) eventplot.pl       $(BINDIR)/eventplot
	$(INSTALL) writeFresParams.pl $(BINDIR)/writeFresParams
	$(INSTALL) plotHitSummary.pl  $(BINDIR)/plotHitSummary
	$(INSTALL) hitSort.pl         $(BINDIR)/hitSort
	$(INSTALL) rmin_bmax_from_stats.py $(BINDIR)/rmin_bmax_from_stats
	$(INSTALL) rip_stat.pl        $(BINDIR)/rip_stat
	$(INSTALL) rip_snr.pl	      $(BINDIR)/rip_snr
	$(INSTALL) rip_diff.pl        $(BINDIR)/rip_diff
	$(INSTALL) getRateFromStats.pl $(BINDIR)/getRateFromStats

#########################################################
.PHONY: uninstall

uninstall:

# compiled C code
	$(RM) $(BINDIR)/{fresnelBox,fresnelT,addKBO,detect,komplete}
	$(RM) $(BINDIR)/{makeKernel,offsetPattern,basis,hideKBO}
	$(RM) $(BINDIR)/{xcorrelate,xchi,fresnelTA,hideKBO_TA}
	$(RM) $(BINDIR)/{fits2ascii}

# perl scripts
	$(RM) $(BINDIR)/{fresplot,plotBoxes,animFresPatt,plotStars}
	$(RM) $(BINDIR)/{hitSort,plotHitSummary}
	$(RM) $(BINDIR)/{hitSummary,eventplot,writeFresParams}
	$(RM) $(BINDIR)/{rmin_bmax_from_stats,rip_stats,rip_snr,rip_diff}

#########################################################
.PHONY: test

test:
	./test.pl

