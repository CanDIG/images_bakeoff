FROM alpine:3.4
MAINTAINER Jonathan Dursi <jonathan@dursi.ca>
LABEL Description="Alpine install of bedtools"

## get needed packages

RUN apk update \
    && apk add wget build-base bash zlib-dev python

## Download and build
ENV VERSION 2.26.0

RUN cd /tmp \
    && wget -nv --no-check-certificate https://github.com/arq5x/bedtools2/releases/download/v2.26.0/bedtools-$VERSION.tar.gz \
    && gunzip bedtools-$VERSION.tar.gz \
    && tar xf bedtools-$VERSION.tar \
    && rm bedtools-$VERSION.tar \
    && cd bedtools2 \
    && make \
    && mv bin/* /usr/local/bin 
    
## cleanup 

RUN apk del wget build-base python \
    && apk cache clean

RUN rm -rf /tmp/bedtools2

ENTRYPOINT ["/usr/local/bin/bedtools"]
