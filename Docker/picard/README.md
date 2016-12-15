# Picard Dockerfile

Note that the dockerfile at https://github.com/broadinstitute/picard doesn't work, and
likely hasn't for some time.

In an effort to make the image as small as possible, used lightweight alpine linux
( https://alpinelinux.org/ ) rather than ubuntu; result is still quite sizable:

```
REPOSITORY          TAG                 IMAGE ID            CREATED             SIZE
picard              latest              6d5ff528ef1c        16 hours ago        682.8 MB
```

although the image could still benefit from removing gradle, the absurdly heavyweight 
Java build system used.
