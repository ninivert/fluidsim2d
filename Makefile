CXX=g++
CC=gcc
LD=${CC}

DEBUG=-O1 -ggdb -fsanitize=address -DDEBUG -fno-omit-frame-pointer
OPTIM=-O3 -march=native
OPTIM+= -DDEBUG

RELEASE=0

CFLAGS+= -Wall -Wextra -Wpedantic
CFLAGS+= -fopenmp
ifneq ($(RELEASE), 1)
CFLAGS+= ${DEBUG}
else
CFLAGS+= ${OPTIM}
endif
LDFLAGS+= $(CFLAGS) -lm

main_perf: main_perf.o
main_gui: CFLAGS+= -Iglad/include $(shell pkg-config --cflags glfw3 gl)
main_gui: LDFLAGS+= -ldl $(shell pkg-config --libs glfw3 gl)
main_gui: main_gui.o glad/src/gl.c

.phony: clean run_gui

clean:
	rm -f *.o main

run_gui: clean main_gui
	ASAN_OPTIONS=fast_unwind_on_malloc=0 ./main_gui

# grind: main
# 	valgrind --leak-check=full --show-reachable=yes ./main
