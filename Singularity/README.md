# Singularity images

Based on the Singularity package (http://singularity.lbl.gov) out of Lawrence Berkeley labs. 
Non-docker-compatibile, less "containing" of process/network/IO, no root-owned daemon processes,
all processes run as user.

The singularity definition scripts given assume you are building on the same distro as you
are building for, but once built, the images should work on any linux system.

Root is required to create and build the images, but any user can run the images, and they
are executed as those images.  e.g.:

```
$ sudo singularity create /tmp/Ubuntu-BWA.img
$ sudo singularity bootstrap /tmp/Ubuntu-BWA.img ubuntu16-bwa-0.7.15.def
$ singularity run /tmp/Ubuntu-BWA.img

Program: bwa (alignment via Burrows-Wheeler transformation)
Version: 0.7.15-r1140
Contact: Heng Li <lh3@sanger.ac.uk>

Usage:   bwa <command> [options]

Command: index         index sequences in the FASTA format
         mem           BWA-MEM algorithm
         fastmap       identify super-maximal exact matches
...

$ id
uid=1000(vagrant) gid=1000(vagrant) groups=1000(vagrant),4(adm),24(cdrom),27(sudo),30(dip),46(plugdev),110(lxd),115(lpadmin),116(sambashare)

$ singularity exec /tmp/Ubuntu-BWA.img id
uid=1000(vagrant) gid=1000(vagrant) groups=1000(vagrant),4(adm),24(cdrom),27(sudo),30(dip),46(plugdev),110,115,116
```
