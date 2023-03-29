# ising

```
scl enable devtoolset-11 bash
```

Compilation instructions:

``` bash
git clone [the fork]
cd ising
bash config/bootstrap
mkdir build
cd build
../configure
make
```

## Configuration options:

```
  CXXFLAGS=-O3
	                      Optimization turned on
  --enable-measurement-cache
                          Enable cache for measurement
  --enable-local-energy-change
                          Enable local energy change measurement
  --enable-lookup-table-accept-reject
                          Enable lookup table for accept/reject step
```

try making different build and compare times:

``` bash
# no optimization
mkdir buildNoOpt
cd buildNoOpt
../configure
make
cd ..

#optimization
mkdir buildOpt
cd buildOpt
../configure CXXFLAGS=-O3
make
cd ..

#optimization and local energy change estimates
mkdir buildOptLocalChange
cd buildOptLocalChange
../configure CXXFLAGS=-O3 --enable-local-energy-change
make
cd ..

...

```


## Input parameters

To make things simple, the input can be passed from screen:
```
L? 30
beta? 0.54
nEvol? 1000
seed? 23634542
```

One can read from an input file:
```bash
bin/ising < input
```
where input contains simply:

```
30
0.54
1000
23634542
```

Outputs two files, `magnetization.dat` and `energy.dat`, containing for each line the evolution and measurement.


## Adidng parallelization

```
#pragma omp parallel for reduction(+:energy)
```

Showing the number of threads employed
```
#include <omp.h>

...

/// Setup the simulation
void setup()
{
  printf("nthreads: %d\n",omp_get_max_threads());
...
```


Choose number of threads from terminal, before launching the code:
```
$ export OM_NUM_THREADS=2
```
