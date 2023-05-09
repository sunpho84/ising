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

## Adding parallelization of the measurements

To parallelize the energy calculation:
```
#pragma omp parallel for reduction(+:energy)
  for(Site site=0;site<V;site++)
    for(Dir dir=0;dir<2;dir++)
      energy-=conf[site]*conf[neighs[site][dir]];

```

Showing the number of threads employed
```
#include <omp.h>

...

/// Setup the simulation
void setup()
{
  printf("nthreads: %d\n",omp_get_max_threads());

Choose number of threads from terminal, before launching the code:
```
$ export OMP_NUM_THREADS=2
```

## Gnuplot commands

Save a file with gnuplot

```
set terminal png 
set output "file.png"
plot "magnetization.dat" w l
```

Plot multiple files

```
plot "magnetization1.dat" w l, "magnetization2.dat" w l
```

## Parallelization of the algorithm

The global generator needs to become local (per-site):

``` c++
mt19937 gen;
```

must be changed into:

``` c++
vector<mt19937> gen;
```

This must be properly resized and initialized:

``` c++
  gen.seed(inputSeed);

```

into

``` c++
  gen.resize(V);
  for(Site site=0;site<V;site++)
    gen[site].seed(inputSeed+site);
```


The site generator must be passed explicitly to the routine which uses it:

``` c++
double drawUniformNumber()
{
...
```

into:

``` c++
double drawUniformNumber(mt19937& gen)
{
...
```

A structure holding the update result must be created:

``` c++
/// Holds the results of the update of a site
struct SiteUpdRes
{
  bool acc;
  Spin dM;
  int dE;
};
```

To be returned from the accept/reject:
