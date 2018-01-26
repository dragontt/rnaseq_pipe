# Human RNA-seq pipeline

### Package installation
Follow instructions from https://github.com/griffithlab/rnaseq_tutorial/wiki/Installation. Note: use OSX version instead of Linux, if working on a Mac.

### Configure tools
Required for each run. This configures the paths for executable tools.
```
export PATH=tools/hisat2-2.0.0-beta/:$PATH
export PATH=tools/samtools-1.3.1/:$PATH
export PATH=tools/stringtie-1.3.3b.OSX_x86_64/:$PATH
```

### Build genome index
This builds the genome index by `hisat2-build` tool.
```
hisat2-build GRCh38/chr22_ERCC92.fa GRCh38/hisat2_index/chr22_ERCC92
```

### Make job script to run alignment
This automatically makes the job scripts to process all samples based on the sample summary file.
```
python tools/build_job_script.py -s metadata/sample_summary.txt -i GRCh38/hisat2_index/chr22_ERCC92 > job_scripts/alignment.makefile
```

### Run jobs
It uses `make` mechanism to process samples with parallelization if applicable. Define number of jobs running in parallel with `-j`. The configuration of tool paths may require adjustment if working with task management system, e.g. dispatching distributed jobs.
```
make -j 4 -f job_scripts/alignment.makefile all > log/alignment.out 2> log/alignment.err
```