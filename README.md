# scChrBin
Single-cell transcriptome gene expression values are binned according to chromatin position

# Use

Our tool bins each chromatin into an equal number of bins based on the length of the chromatin, using the gene annotation files (GTF file) and the lengths of each chromatin.

```python
# example
python main.py --chr_length_file chromelength.txt --gene_gtf_file Homo_sapiens.GRCh37.87.gtf  --width 10 --input_file rpkm.txt
# get output file: rpkm_scChrBin_bin_counts.txt

# help
python main.py --help
Usage: main.py [OPTIONS]

# Options:
#   --chr_length_file PATH  Path to the chromosome length file.  [required]
#   --gene_gtf_file PATH    Path to the gene GTF gene name file.  [required]
#   --width INTEGER         gene bin width.  [required]
#   --input_file PATH       Path for the input file.  [required]
#   --help                  Show this message and exit.
```

Users can simplify the operation through our main.sh to obtain data with different numbers of bins. Users can modify the maximum bin number, starting bin number, and step size in the main.sh script to customize the binning data.

```
sh main.sh
```

For optimal bin number, we recommend finding the optimal number of bins by making predictions on the data using a machine learning classifier. By utilizing the [Feature-scML](https://github.com/liameihao/Feature-scML) tool, it is possible to automate machine learning for binned data, and then proceed with further calculations using the optimal number of bins. 

# Requirements

- ray=1.11.0
- pandas
- numpy
- click

# Reference

Feature-scML：**Liang P**, Wang H, Liang Y, et al. Feature-scML: An Open-source Python Package for the Feature Importance Visualization of Single-Cell Omics with Machine Learning[J]. Current Bioinformatics, 2022, 17(7): 578-585 (https://doi.org/10.2174/1574893617666220608123804).
