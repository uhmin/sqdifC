CFLAGS= -Wall -lm -O -D_FILE_OFFSET_BITS=64
OBJECT= -O -D_FILE_OFFSET_BITS=64

all: sqdifC.exe

sqdifC.exe: sqdif_c.o
        gcc $(CFLAGS) -o $@ $+

clean:
        rm sqdifC.exe *.o

%.o: %.c
        gcc $(OBJECT) -c -o $@ $+
