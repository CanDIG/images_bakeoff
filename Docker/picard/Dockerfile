FROM alpine:3.4
MAINTAINER Jonathan Dursi <jonathan@dursi.ca>

## Using as a starting point:
### http://sirile.github.io/2015/01/01/java-8-minicontainer-with-busybox.html  - busybox w/ Java
### https://github.com/broadinstitute/picard/blob/master/Dockerfile - picard   

## need wget to get files
RUN apk update && apk add wget

## get git to download picard
RUN apk add git 
RUN git config --global http.sslVerify false

RUN mkdir /usr/picard

RUN cd /usr \
    && git clone https://github.com/broadinstitute/picard.git \
    && cd picard \
    && mkdir /usr/picard/lib 

RUN cd /usr/picard \
    && git clone https://github.com/samtools/htsjdk.git 

RUN apk add bash openjdk8 openjdk8-jre
RUN apk add --update --no-cache libstdc++

RUN cd /usr/picard \
    && ./gradlew -g /tmp shadowJar  \
    && mv build/libs/picard.jar picard.jar \
    && rm -rf htsjdk \
    && rm -rf src \
    && rm -rf lib \
    && rm build.xml \
    && ./gradlew clean

RUN mkdir /usr/working
WORKDIR /usr/working
