#!/usr/bin/env python

from typing import Dict, List, Set, Union, Tuple, Iterable
from collections import Counter, defaultdict

import click
import pysam
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from pysam.libcalignmentfile import AlignmentFile



def midchar(s: str) -> str:
    '''Extract middle character

    >>> midchar('ABCDE')
    'C'
    >>> midchar('GGCTCCATCCTTAGACTTGGTCGGTAAAATCTA')
    'T'
    >>> midchar('GGCTCCATCCTTAGACCTGGTCGGTAAAATCTA')
    'C'
    '''
    return s[int(len(s)/2)]


def _parse_scheme_marker(header: str, sequence: str) -> Tuple[int, Dict[str, str]]:
    '''Parse marker info from fasta header and sequence
    
    >>> _parse_scheme_marker('negative1234-1.1.2', 'AGTGA')
    (1234, {'-': 'T', 'subtype': '1.1.2'})
    >>> _parse_scheme_marker('1234-1.1.2', 'AGCGA')
    (1234, {'+': 'C', 'subtype': '1.1.2'})
    '''
    out = {}
    idx, subtype = header.split('-')
    i = int(idx.replace('negative', ''))
    snv = midchar(sequence)
    if idx.startswith('negative'):
        out['-'] = snv
    else:
        out['+'] = snv
    out['subtype'] = subtype
    return i, out


def parse_scheme(scheme_path) -> Dict[int, Dict[str, str]]:
    out = {}
    with open(scheme_path) as f:
        for h, s in SeqIO.FastaIO.SimpleFastaParser(f):
            i, d = _parse_scheme_marker(h, s)
            p = out.get(i, {})
            out[i] = {**p, **d}
    return out


def find_markers(bam: AlignmentFile, positions: Iterable[int]) -> Dict[int, Dict[str, int]]:
    marker_pos_counts = {}
    ref_name = bam.get_reference_name(0)
    for idx in sorted(positions):
        for pcol in bam.pileup(ref_name, start=idx-1, stop=idx):
            if pcol.pos != idx-1:
                continue
            counter = Counter((pread.alignment.query_sequence[pread.query_position] for pread in pcol.pileups if not pread.is_del and not pread.is_refskip))
            marker_pos_counts[pcol.pos] = dict(counter)
    return marker_pos_counts


def snippy_marker_results(snippy_consensus_subs_fa: SeqRecord,
                          marker_pos_counts: Dict[int, Dict[str, int]]) -> pd.DataFrame:
    '''Snippy results for SNV markers

    Extract final Snippy call at each marker position from the Snippy 
    .consensus.subs.fa file. For each marker position, report count for each
    possible base. Output results in a Pandas DataFrame
    '''
    records = []
    for pos, marker_count in marker_pos_counts.items():
        records.append(
            dict(position=pos,
                 snippy_call=snippy_consensus_subs_fa.seq[pos],
                 snippy_A=marker_count.get('A', 0),
                 snippy_C=marker_count.get('C', 0),
                 snippy_G=marker_count.get('G', 0),
                 snippy_T=marker_count.get('T', 0),
                 snippy_depth=sum(marker_count.values())))
    return pd.DataFrame(records)



def bh_base_counts(df_bh: pd.DataFrame) -> Dict[int, Dict[str, int]]:
    '''BioHansel results to dict of position to dict of base counts'''
    bh = defaultdict(dict)
    for i, row in df_bh.iterrows():
        pos = row.refposition
        s = midchar(row.seq)
        bh[pos][s] = row.freq
    return dict(bh)



def biohansel_marker_results(bh_results: Dict[int, Dict[str, int]]) -> pd.DataFrame:
    '''BioHansel results into base calls table'''
    records = []
    for pos, marker_count in bh_results.items():
        records.append(
            dict(position=pos-1,
                 biohansel_A=marker_count.get('A', 0),
                 biohansel_C=marker_count.get('C', 0),
                 biohansel_G=marker_count.get('G', 0),
                 biohansel_T=marker_count.get('T', 0),
                 biohansel_depth=sum(marker_count.values())))
    return pd.DataFrame(records)

@click.command()
@click.option('-b', '--biohansel-results', required=True, help='BioHansel results')
@click.option('-B', '--bam-file', required=True, help='Snippy BAM file')
@click.option('-c', '--snippy-consensus-subs-fa', required=True, help='Snippy consensus FASTA with subsitutions only')
@click.option('-s', '--scheme-fasta', required=True, help='BioHansel scheme FASTA')
@click.option('-o', '--output-csv', required=True, help='Output CSV table path')
def main(biohansel_results, bam_file, snippy_consensus_subs_fa, scheme_fasta, output_csv):
    bam = pysam.AlignmentFile(bam_file)
    scheme = parse_scheme(scheme_fasta)
    seq_rec = SeqIO.read(snippy_consensus_subs_fa, 'fasta')
    snippy_counts = find_markers(bam, scheme.keys())
    df_snippy = snippy_marker_results(seq_rec, snippy_counts)
    df_bh_detailed_results = pd.read_csv(biohansel_results, sep='\t')
    bh_counts = bh_base_counts(df_bh_detailed_results)
    df_bh = biohansel_marker_results(bh_counts)
    df_out = pd.merge(df_bh, df_snippy, how='outer', on='position')
    df_out.to_csv(output_csv)



if __name__ == '__main__':
    main()
