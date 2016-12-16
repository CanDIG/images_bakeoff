# bedtools in Docker

Two ways of getting bedtools into a docker - a small-ish image based on 
Alpine linux that builds statically, and then a very small image containing
only those static executables, provdied here (and obtained from building the
bedtools-alpine docker and copying the files out).

```
REPOSITORY                             TAG                 IMAGE ID            CREATED             SIZE
bedtools-static                        latest              1a0b948fb6df        10 hours ago        8.47 MB
bedtools-alpine                        latest              de2e4f971255        12 hours ago        263 MB
```
