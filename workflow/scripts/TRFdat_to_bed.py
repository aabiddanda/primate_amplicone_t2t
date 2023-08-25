#!/usr/bin/env python
# Downloaded from: https://github.com/hdashnow/TandemRepeatFinder_scripts/blob/master/TRFdat_to_bed.py on 8/22/23

def main():
    # Parse command line arguments
    datfile = snakemake.input["tandem_repeat_dat"]
    bedfile = snakemake.output["tandem_repeat_bed"]
    chrom = "chrY"
    with open(bedfile, 'w') as bed:
        with open(datfile, 'r') as dat:
            for line in dat:
                splitline = line.split()
                if line.startswith("Sequence:"):
                    pass
                    # chrom = line.split()[1]
                else:
                    # Catch index errors when line is blank
                    try:
                        # Check if in header sequence (all non-header lines start with an int: start pos)
                        try:
                            int(splitline[0])
                        except ValueError:
                            continue
                        start = splitline[0]
                        end = splitline[1]
                        motif = splitline[13]
                        bed.write('\t'.join([chrom,start,end,motif]) + '\n')
                    except IndexError:
                        pass

if __name__ == '__main__':
    main()
