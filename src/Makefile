CXX=g++
CC=gcc
LD=${CC}

DEBUG=-O1 -ggdb -fsanitize=address -DDEBUG -fno-omit-frame-pointer
OPTIM=-O3 -march=native
OPTIM+= -DDEBUG

RELEASE=0

CFLAGS+= -Wall -Wextra -Wpedantic
ifneq ($(RELEASE), 1)
CFLAGS+= ${DEBUG}
else
CFLAGS+= ${OPTIM}
endif

# only for main target
CFLAGS+= -Iglad/include

CFLAGS+= $(shell pkg-config --cflags glfw3 gl)
LDFLAGS+= $(CFLAGS) -lm -ldl
LDFLAGS+= $(shell pkg-config --libs glfw3 gl)

.phony: clean run

main: main.o glad/src/gl.c

clean:
	rm -f *.o main

run: clean main
	ASAN_OPTIONS=fast_unwind_on_malloc=0 ./main

grind: main
	valgrind --leak-check=full --show-reachable=yes ./main
