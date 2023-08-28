import numpy as np 
from pyfaidx import Fasta

def obtain_canonical_seq(gene):
    """Obtain a canonical gene sequence."""
    pass


def obtain_seq_qual(seq, n_penalty=0.1, lower_penalty=0.5):
    """Obtain quality score for sequence."""
    assert len(seq) > 0
    assert n_penalty < 1
    assert lower_penalty < 1 
    score = 0.0
    for c in seq
        if c.islower():
            score += lower_penalty
        elif c == "N":
            score += n_penalty
        else:
            score += 1.0
    return score / len(seq)

def create_gene_seq_dict(fasta_fp, gene_tab_fp):
    """From a fasta file isolate the actual underlying gene structure."""
    gene_seq_dict = {}
    with open(gene_tab_fp, "r") as genes:
        # Read the first line of the genes file ... 
        genes.readline()
        for gene_str in  genes:
            [start, end, gene] =  gene_str.split()
            if gene not in gene_seq_dict:
                gene_seq_dict[gene] = {}
            else:
                gene_seq_dict[gene][f'chrY:{start}-{end}']["seq"] = ""
    # Iterate through the ampliconic FASTA
    amplicon_fasta = Fasta(fasta_fp)
    for g in gene_seq_dict:
        for x in gene_seq_dict[g]:
            try:
                cur_seq = amplicon_fasta[x]
                gene_seq_dict[g][x]["seq"] = x[:].seq
                gene_seq_dict[g][x]["qual"] = obtain_seq_qual(x[:].seq)
            except KeyError:
                pass
    return gene_seq_dict


if __name__ == '__main__':
    #1. Read in all of the files and create a gene-seq-dictionary
    
    #2. Identify the highest quality sequence for a given family? (or do some kind of clique assessment across control sequences)

    #3. Output independent fasta files for each of the canonical genes
    pass






