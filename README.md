# scChrBin
Single-cell transcriptome gene expression values are binned according to chromatin position

# Use
```python
# python main.py chr_length_path gene_gtf_genename_path width counts_path
python main.py chr_length.txt gene_gtf_genename.txt 10 counts.txt
# chr_length_path: chr_length.txt (format:txt)
# gene_gtf_genename: gene transfer format(GTF) (format:txt)
# width: bin width (example:10)
# counts_path: single cell RNA-Seq expression (format:txt)
```

# Requirements
- ray=1.11.0
- pandas
- numpy

