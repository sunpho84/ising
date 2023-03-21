# ising

Compilation instructions:

``` bash
git clone [the fork]
cd ising
bash config/bootstrap
mkdir build
cd build
../configure CXXFLAGS=-O3
make
```

Configuration options:

``` bash
  --enable-measurement-cache
                          Enable cache for measurement
  --enable-local-energy-change
                          Enable local energy change measurement
  --enable-lookup-table-accept-reject
                          Enable lookup table for accept/reject step
```


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
