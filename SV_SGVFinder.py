from glob import glob
from os.path import join, splitext, basename
import sys
sys.path.append('/path/to/software/SGVFinder/src/')
import os
import pandas as pd
import numpy as np
from ICRA import single_file
from SV_SGVFinder import get_sample_map, work_on_collection
INPUT_FOLDER = ''
OUTPUT_FOLDER = ''

fastqs_1 = glob(join(INPUT_FOLDER, '*_1.fastq'))
fastqs = [f for f in glob(join(INPUT_FOLDER, '*.fastq')) if f not in fastqs_1 and f.replace('_1.fastq', '_2.fastq') not in fastqs_1]

for f in fastqs:
        single_file(INPUT_FOLDER+f+'_paired_1.fastq',INPUT_FOLDER+f+'_paired_2.fastq',OUTPUT_FOLDER,8, True, 1e-6, 100, 10, 100, 100, 60, 1e5, 2e7, 'genomes', False)

os.mkdir(join(OUTPUT_FOLDER, 'browser'))
samp_to_map = {splitext(basename(f))[0]: get_sample_map(f, .01, 100, 10) \
               for f in glob(join(OUTPUT_FOLDER, '*.jsdel'))}

vsgv, dsgv = work_on_collection(samp_to_map, 10, 0.25, 0.95, 0.25, 0.125, 0.02, 0.95, 'betaprime', 0.01, 10, 85, join(OUTPUT_FOLDER, 'browser'))
vsgv.to_pickle(join(OUTPUT_FOLDER, 'vsgv.df'))
dsgv.to_pickle(join(OUTPUT_FOLDER, 'dsgv.df'))

