FROM alpine:3.4
MAINTAINER Jonathan Dursi <jonathan@dursi.ca>
LABEL Description="Alpine install of samtools"

## get needed packages

RUN apk update \
    && apk add wget build-base bash bzip2 zlib-dev ncurses-static ncurses-dev

## Download and build
ENV VERSION 1.3.1

RUN cd /tmp \
    && wget -nv --no-check-certificate -O samtools-$VERSION.tar.bz2 http://github.com/samtools/samtools/releases/download/$VERSION/samtools-$VERSION.tar.bz2 \
    && bunzip2 samtools-$VERSION.tar.bz2 \
    && tar xf samtools-$VERSION.tar \
    && rm samtools-$VERSION.tar \
    && cd samtools-$VERSION \
    && sed -i -e 's/-g -Wall -O2/-g -Wall -O2 -fPIC -static/' Makefile \
    && sed -i -e 's/-g -Wall -O2/-g -Wall -O2 -fPIC -static/' htslib-1.3.1/Makefile \
    && sed -i -e 's/-lz/ -static -static-libgcc -lz/' Makefile \
    && sed -i -e 's/-lz/ -static -static-libgcc -lz/' htslib-1.3.1/Makefile \
    && sed -i -e 's/-rdynamic//' -e 's/-lcurses/-lncursesw/' Makefile \
    && sed -i -e 's/-rdynamic//' -e 's/-lcurses/-lncursesw/' htslib-1.3.1/Makefile \
    && make \
    && make install
    
## cleanup 

RUN apk del wget bzip2

## add wrapper

RUN mkdir /usr/wrapper

COPY samtools /usr/wrapper
COPY samtools_random /usr/wrapper

ENTRYPOINT ["/usr/wrapper/samtools"]
