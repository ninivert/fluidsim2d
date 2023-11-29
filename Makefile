CXX=g++
# CC=gcc
CC=mpicc
LD=${CC}

DEBUG=-O1 -ggdb -fsanitize=address -DDEBUG -fno-omit-frame-pointer
OPTIM=-O3 -march=native
OPTIM+= -DDEBUG

USE_OMP=0
ifeq ($(USE_OMP), 1)
CFLAGS+= -fopenmp
endif

# use gauss-seidel algo (race condition !!)
USE_GS=0
ifeq ($(USE_GS), 1)
CFLAGS+= -DUSE_GS=1
endif

RELEASE=0
CFLAGS+= -Wall -Wextra -Wpedantic
ifneq ($(RELEASE), 1)
CFLAGS+= ${DEBUG}
else
CFLAGS+= ${OPTIM}
endif
LDFLAGS+= $(CFLAGS) -lm

main_perf: main_perf.o solver.o
main_perf_mpi: main_perf_mpi.o solver.o

main_gui: CFLAGS+= -Iglad/include $(shell pkg-config --cflags glfw3 gl)
main_gui: LDFLAGS+= -ldl $(shell pkg-config --libs glfw3 gl)
main_gui: main_gui.o solver.o glad/src/gl.c

.phony: clean run_gui

clean:
	rm -f *.o main_gui main_perf main_perf_mpi

# run_gui: clean main_gui
# 	ASAN_OPTIONS=fast_unwind_on_malloc=0 ./main_gui

# grind: main
# 	valgrind --leak-check=full --show-reachable=yes ./main
