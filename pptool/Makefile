#OPT=-g -O0
OPT=-O2
CFLAGS=$(OPT) --std=c99 -Wall

all: pptool

pptool: pptool.c miniaudio.o
	$(CC) $(CFLAGS) \
		-DBUILD_LINUX \
		-Inanovg/src \
		-I.. -I../3rd \
		$(shell pkg-config --cflags gl sdl2 flac) \
		miniaudio.o \
		$< -o $@ \
		$(shell pkg-config --libs gl sdl2 flac) \
		-Lnanovg/build \
		-pthread \
		-lnanovg \
		-lm -ldl

clean:
	rm -f pptool *.o
