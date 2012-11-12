SRC = gapmis.c functions.c output.c
OBJ = $(SRC:.c=.o)
CC  = gcc

CFLAGS = -g -Wall -msse3 -O3 -fomit-frame-pointer -funroll-loops

all: gapmis

gapmis: $(OBJ)
	$(CC) $(CFLAGS) -o $@ $(OBJ) -lm

clean:
	rm -f gapmis $(OBJ) *~ gapmis.out

$(OBJ) : functions.h EDNAFULL.h EBLOSUM62.h output.h types.h makefile
