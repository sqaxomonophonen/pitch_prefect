#OPT=-g -O0
OPT=-O2
CFLAGS=$(OPT) --std=c99 -Wall

all: pptool

pptool: pptool.c
	$(CC) $(CFLAGS) \
		-DBUILD_LINUX \
		-Ilibsoundio \
		-Inanovg/src \
		-I.. -I../3rd \
		$(shell pkg-config --cflags gl sdl2 flac) \
		$< -o $@ \
		$(shell pkg-config --libs gl sdl2 flac) \
		-Lnanovg/build \
		libsoundio/build/libsoundio.a \
		-pthread \
		-lpulse -ljack -lasound \
		-lnanovg \
		-lm

clean:
	rm -f pptool
