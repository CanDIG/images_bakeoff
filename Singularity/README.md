# Singularity images

Based on the Singularity package (http://singularity.lbl.gov) out of Lawrence Berkeley labs. 
Non-docker-compatibile, less "containing" of process/network/IO, no root-owned daemon processes,
all processes run as user.

The singularity definition scripts given assume you are building on the same distro as you
are building for, but once built, the images should work on any linux system.
