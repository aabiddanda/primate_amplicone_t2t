import numpy as np 
import pandas as pd 

def get_fasta_length(fasta_fp):
    seq_len = 0
    with open(fasta_fp, "r") as f:
        for i, line in enumerate(f):
            if i > 0:
                seq_len += len(line)
    return seq_len 

def rescale_positions(fasta_fps, dat_fps, out_fp):
    assert len(fasta_fps) == len(dat_fps)
    pos = 0
    cumulative_lines = []
    for cur_fasta_fp, cur_dat_fp in zip(fasta_fps, dat_fps):
        n = get_fasta_length(cur_fasta_fp)
        with open(cur_dat_fp, "r") as dat_fp:
            for i,line in enumerate(dat_fp): 
                line = line.split()
                if i > 0:
                    line[0] = str(int(line[0]) + pos)
                    line[1] = str(int(line[1]) + pos)
                cumulative_lines.append(' '.join(line))
        pos += n
    
    with open(out_fp, "w") as out:
        for line in cumulative_lines:
            out.write(line+"\n")
    
if __name__ == '__main__':
    rescale_positions(snakemake.input["fastas"], snakemake.input["trf_dats"], snakemake.output["meta_trf_dat"])








