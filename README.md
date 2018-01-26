# Human RNA-seq pipeline

### Package installation
Follow instructions from https://github.com/griffithlab/rnaseq_tutorial/wiki/Installation
*** Use OSX version instead of Linux if needed

### Configure tools
```
export PATH=tools/hisat2-2.0.0-beta/:$PATH
export PATH=tools/samtools-1.3.1/:$PATH
export PATH=tools/stringtie-1.3.3b.OSX_x86_64/:$PATH
```

### Build genome index
```
./tools/hisat2-2.0.0-beta/hisat2-build GRCh38/chr22_ERCC92.fa GRCh38/hisat2_index/chr22_ERCC92
```

### Make job script to run alignment
```
python tools/build_job_script.py -s metadata/sample_summary.txt -i GRCh38/hisat2_index/chr22_ERCC92 > job_scripts/alignment.makefile
```