# Images Bakeoff

A benchmarking exercise to examine the suitability for the CanDIG project of various 
ways of bundling applications.  Criteria will be:

- Performance
    - startup/teardown time
    - I/O latency
    - I/O throughput
    - memory latency
- Maintainability
    - automatiability of updated builds
    - Ease of installing need s/w infrastructure at sites
- Generality
    - Will it cover all use cases
- Security

We will be examining:
- VMs
- Containers:
    - [Docker](http://docker.com)
    - [rkt+CoreOS](https://coreos.com/rkt/)
    - [LXD](https://www.ubuntu.com/cloud/lxd)
    - [Singularity](http://singularity.lbl.gov)
    - [Shifter](https://github.com/NERSC/shifter)
- Application packagers
    - [AppImage](http://appimage.org)
    - [Snaps](http://snapcraft.io)

The benchmark is a remapping benchmark, where GRCh37-aligned BAMs
are converted to FASTQ (using Picard - I/O throughput), GRCh38 is
indexed (samtools - I/O throughput), and the FASTQs are re-aligned
(using BWA mem - memory throughput), and then coverage is measured
at an unsorted list of random positions (bedtools - I/O random
access).  This involves two native executables (bedtools and bwa) and
one large java application (Picard).  We will also examine tools for
encapsulating the entire pipeline.

The steps in the benchmark suite are (for native executables):

```
$ wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz \
    && gunzip hg38.fa.gz \
    && bgzip hg38.fa

# Time for IO throughput
$ time samtools faidx hg38.fa.gz

# Time for memory access (and I/O; reduce blocksize with -b N, N < 10000000, for more I/O ops)
$ bwa index hg38.fa.gz

# Remove unpaired reads (test for I/O throughput):
$ samtools view -hf 0x2 input.bam  > test.bam

# Time for IO throughput
$ java -Xmx4g -jar ${PATH_TO_PICARD}/picard.jar SamToFastq VALIDATION_STRINGENCY=LENIENT I=input.bam F=file_1.fastq F2="file_2.fastq"

# Time for memory, IO; test for composability (can also samtools sort standalone with very small memory (-m 10K) for high-IOPS/lots of files)
$ bwa mem hg38.fa.gz file_1.fastq file_2.fastq | samtools view -bS - | samtools sort - > file.bam
```

Benchmark datasets are thousand genomes BAMs:

* Test sets (~400 MB)
    - [Low coverage, phase1, chr 20 HG00119](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/data/HG00119/alignment/HG00119.chrom20.ILLUMINA.bwa.GBR.low_coverage.20101123.bam)
    - [Low coverage, phase1, chr 20 NA19700](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/data/NA19700/alignment/NA19700.chrom20.ILLUMINA.bwa.ASW.low_coverage.20101123.bam)
    - [Low coverage, phase1, chr 20 HG00513](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/data/HG00513/alignment/HG00513.chrom20.ILLUMINA.bwa.CHS.low_coverage.20101123.bam)
* Small (~ 1GB)
    - [Low coverage, phase1, unmapped reads HG00119](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/data/HG00119/alignment/HG00119.unmapped.ILLUMINA.bwa.GBR.low_coverage.20101123.bam)
    - [Low coverage, phase1, unmapped reads NA19700](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/data/NA19700/alignment/NA19700.unmapped.ILLUMINA.bwa.ASW.low_coverage.20101123.bam)
    - [Low coverage, phase1, unmapped reads HG00513](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/data/HG00513/alignment/HG00513.unmapped.ILLUMINA.bwa.CHS.low_coverage.20101123.bam)
* Medium (~15GB)
    - [Low coverage, phase1, WGS HG00119](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/data/HG00119/alignment/HG00119.mapped.ILLUMINA.bwa.GBR.low_coverage.20101123.bam)
    - [Low coverage, phase1, WGS NA19700](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/data/NA19700/alignment/NA19700.mapped.ILLUMINA.bwa.ASW.low_coverage.20101123.bam)
    - [Low coverage, phase1, WGS HG00513](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/data/HG00513/alignment/HG00513.unmapped.ILLUMINA.bwa.CHS.low_coverage.20101123.bam)
* Large (~60GB)
    - [High coverage, HGSV, WGS HG00513](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/data/CHS/HG00513/high_cov_alignment/HG00513.alt_bwamem_GRCh38DH.20150715.CHS.high_coverage.cram)
