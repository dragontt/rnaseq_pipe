# Human RNA-seq pipeline

### Dependency installation
Follow instructions in [RNA-seq tutorial](https://github.com/griffithlab/rnaseq_tutorial/wiki/Installation). Note: use OSX version instead of Linux, if working on a Mac.

Additional tools used in this pipeline are [HTSeq](http://htseq.readthedocs.io/en/master/count.html#usage) for calculating raw read counts, [VarScan](http://varscan.sourceforge.net/using-varscan.html#v2.3_mpileup2snp) for tumor-specific variant calling, and [VEP](http://useast.ensembl.org/info/docs/tools/vep/script/vep_options.html) for variant annotation.

### Configure tools
Update `$HOME/.profile` to configure the paths for executable tools. For example:
```
export PATH=<path/to/rnaseq_pipe>/tools/<tool/bin>/:$PATH
```

### Make directories
```
mkdir -p {metadata,log,job_scripts,report,sequence,alignment/hisat2,expression/htseq,var_calling/varscan_vep}
```

### Build genome index
This builds the genome index by `hisat2-build` tool. You may also swap this for other aligners, and use their genome indexing tools.
```
hisat2-build GRCh38/chr22_ERCC92.fa GRCh38/hisat2_index/chr22_ERCC92
```

### Create job script for batch processing
This step automatically creates the job scripts (in makefile fashion) that is used for processing a batch of samples annotated in the sample summary file.
```
python tools/build_job_script.py -s metadata/sample_summary.txt -i GRCh38/hisat2_index/chr22_ERCC92 -r GRCh38/genes_chr22_ERCC92.gtf > job_scripts/<process_batch>.makefile
```

### Run jobs
This step uses the `make` mechanism as instructed by makefile (from above) to process samples in parallelized fashion. Define number of jobs running in parallel with `-j`. The configuration of tool paths may require adjustment if working with task management system, e.g. dispatching distributed jobs.
```
make -j 4 -f job_scripts/<process_batch>.makefile all > log/<process_batch>.out 2> log/<process_batch>.err
```