import pandas as pd
import pyranges as pr


def _bed2ranges(bed: pd.DataFrame) -> pr.PyRanges:
    names_map = {
        'chr': 'Chromosome',
        'start': 'Start',
        'end': 'End',
        'name': 'Name',
        'score': 'Score',
        'strand': 'Strand'
    }

    return pr.PyRanges(bed.rename(columns=names_map))


def bed_intersect(bed1: pd.DataFrame, bed2: pd.DataFrame) -> pd.DataFrame:
    result = _bed2ranges(bed1).join(
        _bed2ranges(bed2),
        strandedness='same',
        report_overlap=True,
        suffix='_b'
    ).df

    names_map = {
        'Chromosome': 'chr',
        'Start': 'start1',
        'End': 'end1',
        'Name': 'name1',
        'Score': 'score1',
        'Strand': 'strand1',
        'Start_b': 'start2',
        'End_b': 'end2',
        'Name_b': 'name2',
        'Score_b': 'score2',
        'Strand_b': 'strand2'
    }
    result = result.rename(columns=names_map)

    return result
