import pandas as pd
import numpy as np
import re
import ray
import sys
import os


def gene_gtf_info(gene_length_path, gene_gtf_genename_path):
    gene_length = pd.read_csv(gene_length_path,sep="\t")
    gene_gtf = pd.read_csv(gene_gtf_genename_path, sep="\t", comment="#", header=None, low_memory=False)
    pattern = re.compile('gene_name "(.*?)";')
    gene_gtf_gene = gene_gtf[gene_gtf.iloc[:, 2] == "gene"]
    gene_gtf_gene.head()
    gene_gtf_genename = gene_gtf_gene.iloc[:, [0, 3, 4, 6, 8]].copy()
    for i in range(gene_gtf_genename.shape[0]):
        result = pattern.findall(gene_gtf_genename.iloc[i, 4])
        gene_gtf_genename.iloc[i, 4] = result[0]
    gene_gtf_genename.columns = ["chr", "start", "end", "strand","gene_name"]
    return gene_gtf_genename, gene_length

def bin_split(width, gene_gtf_genename, gene_length):
    gene_length_copy = gene_length.copy()
    bins_dict = {}
    for i in range(gene_length_copy.shape[0]):
        bins = []
        length = gene_length_copy.iloc[i, 1]
        for j in range(0, length-length%width+1, length//width):
            bins.append(j)
        bins.pop()
        bins.append(gene_length_copy.iloc[i, 1])
        bins_dict[gene_length_copy.iloc[i, 0]] = bins

    gene_bin = dict()
    for i in bins_dict.keys():
        bins = bins_dict[i]
        simple_chr = gene_gtf_genename[gene_gtf_genename.chr == i].copy()
        simple_chr["bins"] = pd.cut(simple_chr.start, bins=bins, labels=range(len(bins) - 1))
        gene_bin[i] = simple_chr
    return gene_bin, bins_dict

@ray.remote
def cal_bin_counts(one_chr_name, one_chr, counts, width):
    q = np.zeros([1,counts.shape[1]+1])
    for j in range(width):
        one_bins = one_chr[one_chr.bins == j]
        chr_count = counts[counts.gene_name.isin(one_bins.gene_name)]
        res = np.append(np.array([one_chr_name, j]), np.sum(chr_count.iloc[:,1:].values, axis=0))
        q = np.append(q, np.expand_dims(res, axis=0), axis=0)
    return pd.DataFrame(q[1:,], columns=["chr","bins"]+list(counts.columns[1:]))

def bin_sum(counts, gene_bin, width):
    ray.init(num_cpus=30)
    ray.put(gene_bin)
    ray.put(counts)
    ray.put(width)
    result_id = [cal_bin_counts.remote(i,gene_bin[i],counts,width) for i in gene_bin.keys()]
    results = ray.get(result_id)
    ray.shutdown()
    return pd.concat(results, axis=0)

def filename_rename(file_path,suffix):
    base_name, extension = os.path.splitext(file_path)
    return base_name+suffix+extension

def main(gene_length_path, gene_gtf_genename_path, width, counts_path):
    gene_gtf_genename, gene_length = gene_gtf_info(gene_length_path, gene_gtf_genename_path)
    gene_length.to_csv(filename_rename(gene_length_path,"_scChrBin_gene_length"), index=None)
    gene_gtf_genename.to_csv(filename_rename(gene_gtf_genename_path,"_scChrBin_gene_gtf_genename"), index=None)
    #
    gene_bin, bins_dict = bin_split(width, gene_gtf_genename, gene_length)
    data_counts = pd.read_csv(counts_path, sep="\t")
    data_counts.columns = ["gene_name"]+list(data_counts.columns[1:])
    results = bin_sum(data_counts, gene_bin, width)
    results.to_csv(filename_rename(counts_path,"_scChrBin_bin_counts"), index=None)

if __name__ == '__main__':
    gene_length_path = str(sys.argv[1])
    gene_gtf_genename_path = str(sys.argv[2])
    width = int(sys.argv[3])
    counts_path = str(sys.argv[4])
    # run
    print("#################")
    print("start")
    main(gene_length_path, gene_gtf_genename_path, width, counts_path)
    print("end")
    print("#################")

