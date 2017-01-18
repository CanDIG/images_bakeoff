# Images Bakeoff

A benchmarking exercise to examine the suitability for the CanDIG project of various 
ways of bundling applications.  Criteria will be:

- Performance
    - startup/teardown time
    - I/O latency
    - I/O throughput
    - memory latency
- Maintainability
    - automatability of updated builds
    - Ease of installing: amount of needed s/w infrastructure at sites
- Generality
    - Will it cover all use cases
- Sustainability
    - Will it exist/break backward compatibility over the course of 4 years
- Security

We will be examining:
- Application packagers
    - [AppImage](http://appimage.org)
    - [Snaps](http://snapcraft.io)
- Containers on bare metal:
    - [Docker](http://docker.com)
    - [rkt+CoreOS](https://coreos.com/rkt/)
    - [LXD](https://www.ubuntu.com/cloud/lxd)
    - [Singularity](http://singularity.lbl.gov)
    - [Shifter](https://github.com/NERSC/shifter) (maybe?)
- VMs
- The above containers, packages within a VM
- VM within a VM

And we will use bare-metal native executables as a baseline.

The benchmark is a remapping benchmark, where GRCh37-aligned BAMs
are converted to FASTQ (using Picard - I/O throughput), GRCh38 is
indexed (samtools - I/O throughput), and the FASTQs are re-aligned
(using BWA mem - memory throughput), and then coverage is measured
at an unsorted list of random positions (bedtools - I/O random
access).  This involves two native executables (bedtools and bwa) and
one large java application (Picard).  We will also examine tools for
encapsulating the entire pipeline.
