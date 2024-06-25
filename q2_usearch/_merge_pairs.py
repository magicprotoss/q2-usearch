# ----------------------------------------------------------------------------
# Modified from q2-vsearch
#
# Copyright (c) 2016-2024, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os.path
import yaml
from typing import List
import gzip
import tempfile
import shutil
import subprocess
import pandas as pd
import numpy as np
from q2_types.per_sample_sequences import (
    SingleLanePerSampleSingleEndFastqDirFmt,
    SingleLanePerSamplePairedEndFastqDirFmt,
    FastqManifestFormat, YamlFormat)

def run_command(cmd, verbose=True):
    print("Running external command line application. This may print "
          "messages to stdout and/or stderr.")
    print("The command being run is below. This command cannot "
          "be manually re-run as it will depend on temporary files that "
          "no longer exist.")
    print("\nCommand:", end=' ')
    print(" ".join(cmd), end='\n\n')
    subprocess.run(cmd, check=True)


_mp_defaults = {
    'truncqual': 2, # -fastq_trunctail (vsearch default is None)
    'minlen': 64, # -fastq_minlen (vsearch is 1)
    # 'maxns': None, # don't exist in usearch
    'allowmergestagger': False, # -fastq_nostagger (vsearch defalue is false, ussearch default is to trim  overhangs (non-biological sequence))
    'minovlen': 16, # -fastq_minovlen (vsearch 10)
    'maxdiffs': 5, # -fastq_maxdiffs (vsearch 10)
    'minmergelen': None, # -fastq_minmergelen
    'maxmergelen': None, # fastq_maxmergelen
    # 'maxee': None, # do not exist in usearch
    'threads': "auto", # vsearch default is 1
    'percent_identity': 90,
    # 'preset': None
}


def merge_pairs(
    demultiplexed_seqs: SingleLanePerSamplePairedEndFastqDirFmt,
    truncqual: int = _mp_defaults['truncqual'],
    minlen: int = _mp_defaults['minlen'],
    allowmergestagger: bool = _mp_defaults['allowmergestagger'],
    minovlen: int = _mp_defaults['minovlen'],
    maxdiffs: int = _mp_defaults['maxdiffs'],
    percent_identity: int = _mp_defaults['percent_identity'],
    minmergelen: int = _mp_defaults['minmergelen'],
    maxmergelen: int = _mp_defaults['maxmergelen'],
    preset: int = _mp_defaults['preset'],
    threads: str = _mp_defaults['threads'],
) -> (
    SingleLanePerSampleSingleEndFastqDirFmt,
    SingleLanePerSamplePairedEndFastqDirFmt
):
    _, merged, unmerged = _merge_pairs_cli(
        demultiplexed_seqs, truncqual, minlen, allowmergestagger,
        minovlen, maxdiffs, percent_identity, minmergelen, maxmergelen, threads
    )

    # if preset == 'double_reigon_short_overlap':
    #     pass
    # elif preset == 'double_reigon_long_overlap':
    #     pass
    # elif preset == 'signle_reigon_long_overlap':
    #     allowmergestagger = True
    #     maxdiffs = 10
    #     percent_identity = 80
        
    return merged, unmerged


def _merge_pairs_cli(
    demultiplexed_seqs: SingleLanePerSamplePairedEndFastqDirFmt,
    truncqual: int = _mp_defaults['truncqual'],
    minlen: int = _mp_defaults['minlen'],
    allowmergestagger: bool = _mp_defaults['allowmergestagger'],
    minovlen: int = _mp_defaults['minovlen'],
    maxdiffs: int = _mp_defaults['maxdiffs'],
    percent_identity: int = _mp_defaults['percent_identity'],
    minmergelen: int = _mp_defaults['minmergelen'],
    maxmergelen: int = _mp_defaults['maxmergelen'],
    threads = _mp_defaults['threads'],
) -> (
    List[str],
    SingleLanePerSampleSingleEndFastqDirFmt,
    SingleLanePerSamplePairedEndFastqDirFmt
):
    # this function exists only to simplify unit testing

    # create formats
    merged = SingleLanePerSampleSingleEndFastqDirFmt()
    unmerged = SingleLanePerSamplePairedEndFastqDirFmt()

    # create manifests
    merged_manifest = FastqManifestFormat()
    merged_manifest_fh = merged_manifest.open()
    unmerged_manifest = FastqManifestFormat()
    unmerged_manifest_fh = unmerged_manifest.open()

    # write manifest headers
    _write_manifest_header(merged_manifest_fh, add_warning=True)
    _write_manifest_header(unmerged_manifest_fh)

    # generate input reads iterable
    manifest = pd.read_csv(
        os.path.join(
            str(demultiplexed_seqs), demultiplexed_seqs.manifest.pathspec
        ),
        header=0,
        comment='#'
    )

    manifest.filename = manifest.filename.apply(
        lambda x: os.path.join(str(demultiplexed_seqs), x)
    )
    # how usearch 11 handles old illumina 1.3 dataset?
    phred_offset = yaml.load(
        open(os.path.join(
            str(demultiplexed_seqs), demultiplexed_seqs.metadata.pathspec
        )),
        Loader=yaml.SafeLoader
    )['phred-offset']

    id_to_fps = manifest.pivot(
        index='sample-id', columns='direction', values='filename'
    )

    # create a temp folder to avoid corrupting Artifact filepath
    with tempfile.TemporaryDirectory() as temp_dir:
        os.mkdir(os.path.join(temp_dir, 'input'))
        os.mkdir(os.path.join(temp_dir, 'merged'))
        os.mkdir(os.path.join(temp_dir, 'unmerged'))

        for i, (sample_id, (gzipped_fwd_fp, gzipped_rev_fp)) in enumerate(id_to_fps.iterrows()):
            # The barcode id and lane number are not relevant for either format.
            # We might ultimately want to use a dir format other than these which
            # doesn't care about this information.
            # The read number (direction) is only relevant for the unmerged reads.
    
            # prep input fps
            fwd_fp, rev_fp = _unzip_seqs_for_usearch(sample_id, gzipped_fwd_fp, gzipped_rev_fp, os.path.join(temp_dir, "input"))
                
            # prep output fps
            gz_merged_path, fq_merged_path = _get_output_paths(
                merged, sample_id, i, 1, os.path.join(temp_dir, "merged")
            )
            gz_unmerged_fwd_path, fq_unmerged_fwd_path = _get_output_paths(
                unmerged, sample_id, i, 1, os.path.join(temp_dir, "unmerged")
            )
            gz_unmerged_rev_path, fq_unmerged_rev_path = _get_output_paths(
                unmerged, sample_id, i, 2, os.path.join(temp_dir, "unmerged")
            )
    
            # build command
            cmd = [
                'usearch',
                '-fastq_mergepairs', fwd_fp,
                '-reverse', rev_fp,
                '-fastqout', fq_merged_path,
                '-fastqout_notmerged_fwd', fq_unmerged_fwd_path,
                '-fastqout_notmerged_rev', fq_unmerged_rev_path,
                '-fastq_minlen', str(minlen),
                '-fastq_minovlen', str(minovlen),
                '-fastq_maxdiffs', str(maxdiffs)
            ]
            
            if percent_identity is not None:
                cmd += ['-fastq_pctid', str(percent_identity)]
            if truncqual is not None:
                cmd += ['-fastq_trunctail', str(truncqual)]
           # no such option in usearch
           # if maxns is not None:
           #     cmd += ['--fastq_maxns', str(maxns)]
            if minmergelen is not None:
                cmd += ['-fastq_minmergelen', str(minmergelen)]
            if maxmergelen is not None:
                cmd += ['-fastq_maxmergelen', str(maxmergelen)]
            # if maxee is not None:
            #     cmd += ['--fastq_maxee', str(maxee)]
            if threads != "auto":
                cmd += ['-threads', str(threads)]
            if not allowmergestagger:
                cmd.append('-fastq_nostagger')
    
            run_command(cmd)
            
            # remove input files
            os.remove(fwd_fp)
            os.remove(rev_fp)
            
            # zip all output files
            with open(fq_merged_path, 'rb') as f_in:
                with gzip.open(gz_merged_path, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
                    
            os.remove(fq_merged_path)
            
            with open(fq_unmerged_fwd_path, 'rb') as f_in:
                with gzip.open(gz_unmerged_fwd_path, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            
            os.remove(fq_unmerged_fwd_path)
            
            with open(fq_unmerged_rev_path, 'rb') as f_in:
                with gzip.open(gz_unmerged_rev_path, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
                    
            os.remove(fq_unmerged_rev_path)
            
            # if it works why change it TT
    
            merged_manifest_fh.write(
                '%s,%s,%s\n' % (sample_id, gz_merged_path.name, 'forward')
            )
            unmerged_manifest_fh.write(
                '%s,%s,%s\n' % (sample_id, gz_unmerged_fwd_path.name, 'forward')
            )
            unmerged_manifest_fh.write(
                '%s,%s,%s\n' % (sample_id, gz_unmerged_rev_path.name, 'reverse')
            )

    merged_manifest_fh.close()
    unmerged_manifest_fh.close()
    merged.manifest.write_data(merged_manifest, FastqManifestFormat)
    unmerged.manifest.write_data(unmerged_manifest, FastqManifestFormat)

    metadata = YamlFormat()
    metadata.path.write_text(yaml.dump({'phred-offset': phred_offset}))
    merged.metadata.write_data(metadata, YamlFormat)
    unmerged.metadata.write_data(metadata, YamlFormat)

    return cmd, merged, unmerged


def _get_output_paths(format_, sample_id, barcode_id, direction, temp_folder):
    path = format_.sequences.path_maker(
        sample_id=sample_id,
        barcode_id=barcode_id,
        lane_number=1,
        read_number=direction
    )
    # since it's only a fn it doesn't matter here
    tmp_path = os.path.join(temp_folder, sample_id + '_R' + str(direction) + '.fastq')
    return path, tmp_path

def _write_manifest_header(manifest_fh, add_warning=False):
    manifest_fh.write('sample-id,filename,direction\n')
    if add_warning:
        manifest_fh.write('')
        manifest_fh.write('# direction is not meaningful for joined reads\n')
        
def _unzip_seqs_for_usearch(sample_id, gzipped_fwd_fq_fp, gzipped_rev_fq_fp, temp_dir):
    uzipped_fwd_fq_fp = os.path.join(temp_dir, sample_id + '_R1.fastq')
    uzipped_rev_fq_fp = os.path.join(temp_dir, sample_id + '_R2.fastq')
    with gzip.open(gzipped_fwd_fq_fp, 'rb') as f_in:
        with open(uzipped_fwd_fq_fp, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    with gzip.open(gzipped_rev_fq_fp, 'rb') as f_in:
        with open(uzipped_rev_fq_fp, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    return uzipped_fwd_fq_fp, uzipped_rev_fq_fp
