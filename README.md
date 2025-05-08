# GCoM \[ISCA '22\]

This repository contains the source code for [GCoM \[ISCA '22\]](https://doi.org/10.1145/3470496.3527384), a detailed GPU core model for accurate analytical modeling of modern GPUs. Please cite the following paper if you utilize GCoM in your research.

```bibtex
@inproceedings{lee2022gcom,
  author    = {Jounghoo Lee and Yeonan Ha and Suhyun Lee and Jinyoung Woo and Jinho Lee and Hanhwi Jang and Youngsok Kim},
  title     = {{GCoM: A Detailed GPU Core Model for Accurate Analytical Modeling of Modern GPUs}},
  booktitle = {Proc. 49th IEEE/ACM International Symposium on Computer Architecture (ISCA)},
  year      = {2022},
}
```

## Dependencies

* Ubuntu packages
```
$ sudo apt install xutils-dev bison zlib1g-dev flex libglu1-mesa-dev libssl-dev libxml2-dev libxml2-dev
```
* gcc/g++ 8.4.0 or higher

* zlib & gzstream
```bash
third-party/zlib-1.2.12$ prefix=. ./configure
third-party/zlib-1.2.12$ make test
third-party/zlib-1.2.12$ make install prefix=.

third-party/gzstream$ make test
third-party/gzstream$ make
```

* boost serialization library
```
third-party/boost_1_86_0$ ./bootstrap.sh
third-party/boost_1_86_0$ ./b2 --with-serialization
```

## Building & Running GCoM
### Collect SASS trace of an application
Please refer to [Accel-sim Tracer](https://github.com/accel-sim/accel-sim-framework/tree/2260456ea5e6a1420f5734f145a4b7d8ab1d4737) for the SASS trace collection.
Then, compress resulting `*.traceg` files to `*.traceg.gz` using `gzip`.

Otherwise, we provide some sample traces, which can be downloaded using `./get_trace_samples.sh`.

### Build and run model

```
$ make -j [debug]
$ ./bin/GCoM -h
usage : ./bin/GCoM 
 -h    help
 -w [INT] work ID
 -r [path] representative warp index trace path
 -C [path] configuration file path
 -t [path] benchmark trace directory

[work 1] RunGCoMWithSimpleCacheModel options
 -O [path] (Optional) export representative warp selector input to [path]

# run Rodinia NN
$ ./bin/GCoM -w 1 -r repwarp_selector/nn-rodinia-3.1/__data_filelist_4__r_5__lat_30__lng_90/rep_warp_out.bin -t trace_samples/rtx2060/11.0/nn-rodinia-3.1/__data_filelist_4__r_5__lat_30__lng_90/traces/ -C configs/RTX2060.config
...
K1 Result       ,numThreadInst  ,numCycle       ,cpi    ,base   ,comData        ,comStruct      ,memData        ,memStruct      ,idle
Value   ,1327636        ,2182.24        ,0.00164371     ,0.000270707    ,8.44526e-05    ,0.00018843     ,0.00102486     ,3.03547e-05    ,4.49025e-05
Total Result    ,numThreadInst  ,numCycle       ,cpi    ,base   ,comData        ,comStruct      ,memData        ,memStruct      ,idle
Value   ,1327636        ,2182.24        ,0.00164371     ,0.000270707    ,8.44526e-05    ,0.00018843     ,0.00102486     ,3.03547e-05    ,4.49025e-05
```

## Code Organization

```
configs/: hardware configuration file
src/: source code of performance model and cache simulator
repwarp_selector/: place-holder for external representative warp selector
third-party/ : code from other projects. Modified based on our needs.
get_trace_samples.sh: download trace samples
```

