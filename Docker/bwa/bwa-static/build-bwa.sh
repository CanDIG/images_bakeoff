#!/bin/bash
##
## Builds the static bwa executable for the docker.
## Needs to be run on a linux box.
##

## requires: build-essential git

# get musl for static libgc
wget -nv https://www.musl-libc.org/releases/musl-1.1.15.tar.gz
tar -xzf musl-1.1.15.tar.gz
rm musl-1.1.15.tar.gz 
cd musl-1.1.15 && ./configure --prefix=.. && make
cd ..

# compile zlib with musl
wget -nv http://zlib.net/zlib-1.2.8.tar.gz
tar -xzf zlib-1.2.8.tar.gz
rm zlib-1.2.8.tar.gz
cd zlib-1.2.8
CC=../bin/musl-gcc ./configure --prefix=../zlib/1.2.8 --static
make install
cd ..
rm -rf zlib-1.2.8

# download a fixed version of BWA and make
export VERSION=0.7.15
wget -nv https://github.com/lh3/bwa/archive/v$VERSION.tar.gz
tar -xjf v$VERSION.tar.gz
rm v$VERSION.tar.gz
cd bwa-$VERSION
ex -sc '1i|#include <inttypes.h>' -cx kthread.c
ex -sc '1i|#include <sys/types.h>' -cx bwtgap.h
ex -sc '1i|#include <sys/types.h>' -cx bwt_lite.c
make bwa CC=../bin/musl-gcc CFLAGS="-g -Wall -Wno-unused-function -O2 -static -I../zlib/1.2.8/include" LIBS="-lm -L../zlib/1.2.8/lib -lz -lpthread -lrt"
cp bwa ..
cd ..
