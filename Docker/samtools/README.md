# samtools in Docker

Two ways of getting samtools into a docker - a small-ish image based on 
Alpine linux that builds statically, and then a very small image containing
only those static executables, provdied here (and obtained from building the
samtools-alpine docker and copying the files out).

Samtools-alpine has a wrapper that times the executable so that you can get
timings from outside the docker (which will include startup/takedown) and
the executable within the docker.  The static version just has the single 
executable.

```
REPOSITORY                             TAG                 IMAGE ID            CREATED             SIZE
samtools-static                        latest              a4e9707ad374        15 seconds ago       4.87 MB
samtools-alpine                        latest              f8e671f0b91a        17 minutes ago       232 MB
```
