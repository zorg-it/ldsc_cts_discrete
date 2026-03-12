#!/usr/bin/env python
from __future__ import print_function
import pandas as pd
import numpy as np
import argparse
from pybedtools import BedTool
import gzip
from os import path
import resource

### !!! the script is fitted in order to use more bed files, but not for gene set ###
def read_in_chunks(file_object, chunk_size=1073741824):
    """Lazy function (generator) to read a file piece by piece.
    Default chunk size: 1GB."""
    while True:
        data = file_object.read(chunk_size)
        if not data:
            break
        yield data

def gene_set_to_bed(args):
    print('making gene set bed file')
    GeneSet = pd.read_csv(args.gene_set_file, header = None, names = ['GENE'])
    all_genes = pd.read_csv(args.gene_coord_file, delim_whitespace = True)
    df = pd.merge(GeneSet, all_genes, on = 'GENE', how = 'inner')
    df['START'] = np.maximum(1, df['START'] - args.windowsize)
    df['END'] = df['END'] + args.windowsize
    iter_df = [[(str(x1).lstrip('chr')), x2 - 1, x3] for (x1,x2,x3) in np.array(df[['CHR', 'START', 'END']])]
    return BedTool(iter_df).sort().merge()

def make_annot_files(args, bed_for_annot):
    print('making annot file')
    #loading bim
    df_bim = pd.read_csv(args.bimfile,
            delim_whitespace=True, usecols = [0,1,2,3], names = ['CHR','SNP','CM','BP'])
    iter_bim = [[str(x1), x2 - 1, x2] for (x1, x2) in np.array(df_bim[['CHR', 'BP']])]
    bimbed = BedTool(iter_bim)

    #intersecting with bed and extracting relevant info
    annotbed = bimbed.intersect(bed_for_annot, loj=True)
    bp = [x.start + 1 for x in annotbed]
    if args.nomerge:
        annotnames= [x[6] for x in annotbed]
    else:
        annotnames= ['.' if x[5]=='-1' else 'merge' for x in annotbed]
    thisannot= [1 for x in annotbed]
    df_int = pd.DataFrame({'BP': bp, 'ANNOT': thisannot, 'ANNOTNAME': annotnames})

    #merge with bim info
    df_annot = pd.merge(df_bim, df_int, on='BP',how='left')

    #casting to pivot table
    dfannot= df_annot.pivot_table(index=['CHR','SNP','CM','BP'], columns='ANNOTNAME', values='ANNOT').reset_index()
    annotlevels=list(set(annotnames))
    missing_annots=[annot for annot in annotlevels if annot not in dfannot.columns]
    dfannot = dfannot.reindex(['CHR','SNP','CM','BP']+sorted(list(dfannot.columns[4:])+missing_annots), axis=1)
    #cleaning up the pivot table
    if '.' in dfannot.columns:
        print(str(dfannot.loc[:,'.'].sum())+' out of '+str(dfannot.shape[0])+' SNPs without annotations')
        dfannot.drop('.',1,inplace=True)
    dfannot.fillna(0, inplace=True)
    dfannot.sort_values(by=['BP'], inplace=True)
    dfannot.drop_duplicates(subset=['BP'], inplace=True)

    if args.add_base:
        dfannot['base']= 1
    dfannot.iloc[:,4:]= dfannot.iloc[:,4:].astype(int, inplace=True)

    #printing
    if args.annot_file.endswith('.gz'):
        with gzip.open(args.annot_file, 'wb') as f:
            dfannot.to_csv(f, sep = "\t", index = False)
    else:
        dfannot.to_csv(args.annot_file, sep="\t", index=False)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--gene-set-file', type=str, help='a file of gene names, one line per gene.')
    parser.add_argument('--gene-coord-file', type=str, default='ENSG_coord.txt', help='a file with columns GENE, CHR, START, and END, where START and END are base pair coordinates of TSS and TES. This file can contain more genes than are in the gene set. We provide ENSG_coord.txt as a default.')
    parser.add_argument('--windowsize', type=int, help='how many base pairs to add around the transcribed region to make the annotation?')
    parser.add_argument('--named-bed-file', type=str, help='the UCSC bed file with named regions for each of your annotations. Only one named bed can be provided, file is assumed to be sorted with \'sort -k1,1 -k2,2n\'')
    parser.add_argument('--bed-file', type=str, nargs='+', help='the UCSC bed file with the regions that make up your annotation, you can add more bed files')
    parser.add_argument('--nomerge', action='store_true', default=False, help='don\'t merge the bed file; make an annot file with multiple columns containing values for each input bedfile.')
    parser.add_argument('--chunks', action='store_true', default=False, help='read input named-bed-file in chunks of 1GB')
    parser.add_argument('--mem', type=int, default=None, help='maximum memory to be used for this process in mb. when reaching the limit the process is killed, defaults to no limit')
    parser.add_argument('--bimfile', type=str, help='plink bim file for the dataset you will use to compute LD scores.')
    parser.add_argument('--annot-file', type=str, help='the name of the annot file to output.')
    parser.add_argument('--add-base', action='store_true', default=False, help='add a base column annotation for every snp')

    args = parser.parse_args()    
    if args.mem is not None:
        resource.setrlimit(resource.RLIMIT_AS, (args.mem*1049000, args.mem*1049000))

    if args.gene_set_file is not None:
        bed_for_annot = gene_set_to_bed(args.gene_set_file)
    elif args.named_bed_file is not None:
        if args.chunks:
            bed_for_annot=BedTool([])
            with open(args.named_bed_file) as f:
                for piece in f.readlines(1073741824):
                    bed_for_annot = bed_for_annot.cat(BedTool(piece),postmerge=False)
        else:
            bed_for_annot=BedTool(args.named_bed_file)
    else:
        bed_for_annot=BedTool([])
        for f in args.bed_file:
            bed_iter=[l.strip().split('\t') + [path.basename(f)] for l in open(f).readlines()]
            bed_for_annot=bed_for_annot.cat(BedTool(bed_iter),postmerge=False)
        bed_for_annot=bed_for_annot.sort()
    if not args.nomerge:
            bed_for_annot = bed_for_annot.merge()
    make_annot_files(args, bed_for_annot)
