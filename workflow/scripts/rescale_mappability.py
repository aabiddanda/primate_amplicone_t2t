"""Rescale mappability scores for defined regions based on the underlying kmer analysis."""
import polars as pl

def convert_raw_df_to_full(raw_df):
    """Convert the raw data frame to a full dataframe."""
    positions = []
    gc_percentage = []
    mappability = []
    for s,e,m in zip(raw_df['START'].to_numpy(), raw_df['END'].to_numpy(), raw_df['Mappability'].to_numpy()):
        for i in range(s,e):
            positions.append(i)
            gc_percentage.append(50.0)
            mappability.append(m)
    full_df = pl.DataFrame({'Position': positions, 'GCpercentage': gc_percentage, 'Mappability': mappability, 'InformativeSites': 1})
    return full_df

def set_mappability(start, end, anno_df, raw_mappability_df):
    assert start <= end
    df = anno_df.filter((pl.col('Position') <= e) & (pl.col('Position') >= s))
    raw_df = raw_mappability_df.filter((pl.col('START') <= e) & (pl.col('END') >= s))
    if anno_df.shape[0] == 0:
        fix_df = convert_raw_df_to_full(raw_df)
    else:
        fix_df = df.with_columns(pl.lit(1.0).alias('InformativeSites'))
    return fix_df

if __name__ == '__main__':
    anno_df = pl.read_csv(snakemake.input['chrY_rescaled_anno'], separator='\t')
    gene_df = pl.read_csv(snakemake.input['gene_def'], separator='\t')
    mappability_df = pl.read_csv(snakemake.input['mappability_bed'], separator='\t')
    collection_dfs =  []
    for s,e,t,l in zip(gene_df['START'].to_numpy(), gene_df['END'].to_numpy(), gene_df['TYPE'].to_numpy()):
        print(s,e,t)
        fix_df = set_mappability(s,e, anno_df, mappability_df)
        collection_dfs.append(fix_df)
    # tying together all of these results     
    rev_anno_df = pl.concat([anno_df] + collection_dfs).unique(subset='Position', keep='last').sort('Position')
    rev_anno_df.write_csv(snakemake.output['chrY_rescaled_anno'], separator='\t')
        