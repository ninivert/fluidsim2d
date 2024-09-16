# fluidsim2d

a rushed project for [parallel programming](https://edu.epfl.ch/coursebook/en/parallel-programming-PHYS-743)

interactive gui (requires glfw3 gl libs probably):

```sh
make -B RELEASE=1 main_gui && ./main_gui
```

make perf (requires openmp):

```sh
make -B USE_OMP=1 USE_GS=1 RELEASE=1 main_perf && ./main_perf 100 100 1
```

presentation: [pres.pdf](./pres.pdf)
