#OPT=-g -O0
OPT=-O2
CFLAGS=$(OPT) --std=c99 -Wall

all: fungame

fungame: fungame.c sound.h
	$(CC) $(CFLAGS) \
		-DBUILD_LINUX \
		-Ilibsoundio \
		-Inanovg/src \
		-I.. -I../3rd \
		$(shell pkg-config --cflags gl sdl2) \
		$< -o $@ \
		$(shell pkg-config --libs gl sdl2) \
		-Lnanovg/build \
		libsoundio/build/libsoundio.a \
		-pthread \
		-lpulse -ljack -lasound \
		-lnanovg \
		-lm

clean:
	rm -f fungame
