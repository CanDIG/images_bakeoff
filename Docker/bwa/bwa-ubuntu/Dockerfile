FROM ubuntu:16.04
MAINTAINER Jonathan Dursi <jonathan@dursi.ca>
LABEL Description="Full-ubuntu install w/ bwa for comparison"

RUN apt-get update \
    && apt-get install -y wget build-essential zlib1g-dev

ENV VERSION 0.7.15

RUN cd /tmp \
    && wget -nv https://github.com/lh3/bwa/archive/v$VERSION.tar.gz \
    && tar -xzf v$VERSION.tar.gz \
    && rm v$VERSION.tar.gz \
    && cd bwa-$VERSION \
    && CFLAGS="-static -g -Wall -Wno-unused-function -O2" make bwa \
    && mv bwa /usr/local/bin 
    
RUN rm -rf /tmp/bwa-$VERSION

RUN apt-get remove -y build-essential git 

RUN apt autoremove -y

ENTRYPOINT ["/usr/local/bin/bwa"]
