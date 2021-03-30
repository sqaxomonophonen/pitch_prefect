#OPT=-g -O0
OPT=-O2
pkgs=gl sdl2 flac
CFLAGS=$(OPT) --std=c99 -Wall $(shell pkg-config --cflags $(pkgs)) -DBUILD_LINUX -Inanovg/src -I.. -I../3rd

all: pptool

pptool.o: pptool.c ../pitch_prefect.h

pptool: pptool.o miniaudio.o
	$(CC) $(CFLAGS) \
		miniaudio.o \
		$< -o $@ \
		$(shell pkg-config --libs $(pkgs)) \
		-Lnanovg/build \
		-pthread \
		-lnanovg \
		-lm -ldl

clean:
	rm -f pptool *.o
