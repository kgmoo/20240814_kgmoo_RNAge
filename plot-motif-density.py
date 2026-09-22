#!/usr/bin/env python
# coding: utf-8

import os
import sys
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import glob
import pybedtools as bt
from scipy.stats import norm

plt.rcParams.update({'font.size': 20})

MOTIF_BED, DESEQ2, PREFIX = "/scratch/scjp_root/scjp0/porchard/aging-motif-density/data/motifs/NR3C1_known18.bed", "/scratch/scjp_root/scjp0/porchard/aging-motif-density/data/age-association/M.Type_2a.tsv", "M.Type_2a.NR3C1_known18"

TSS = '/nfs/turbo/umms-scjp/shared_data/reference/human/hg38/topmed/tss/hg38.gene-id.tss.bed'

PEAK_CENTER_FLANK = 1000
BINS = 200
PEAK_CENTER_FLANK = int(PEAK_CENTER_FLANK) # use peak center +/- this number of bps

MOTIF = os.path.basename(MOTIF_BED).replace('.bed', '')



deseq2 = pd.read_csv(DESEQ2, sep='\t')
deseq2 = deseq2[['estimate', 'se', 'z', 'p_GLMMNB']].rename(columns={'estimate': 'beta', 'se': 'beta_se', 'z': 'beta_z', 'p_GLMMNB': 'beta_p'})
deseq2[['chrom', 'start', 'end']] = deseq2.index.to_series().str.split(':', expand=True)
deseq2.start = deseq2.start.astype(int)
deseq2.end = deseq2.end.astype(int)
deseq2 = deseq2[['chrom', 'start', 'end', 'beta', 'beta_se', 'beta_z', 'beta_p']].sort_values('beta_z')

peaks = deseq2[['chrom', 'start', 'end']]
fcs = dict(deseq2.beta_z)
peak_rank = {peak: rank for rank, peak in enumerate(deseq2.index)}

# get peak center and flanking
peak_centers = peaks.copy()
peak_centers['center'] = peak_centers[['start', 'end']].mean(axis=1).map(round)
new_peak_to_old_peak = {'{}:{}:{}'.format(chrom, center-PEAK_CENTER_FLANK, center+PEAK_CENTER_FLANK): '{}:{}:{}'.format(chrom, start, end) for chrom, start, end, center in zip(peak_centers.chrom, peak_centers.start, peak_centers.end, peak_centers.center)}
peak_centers.start = peak_centers.center - PEAK_CENTER_FLANK
peak_centers.end = peak_centers.center + PEAK_CENTER_FLANK


motifs = pd.read_csv(MOTIF_BED, sep='\t', header=None, names=['chrom', 'start', 'end'], usecols=[0, 1, 2])
MOTIF_LENGTH = (motifs.end - motifs.start).unique()[0]


# assign peaks as TSS distal or proximal
peaks_bt = bt.BedTool().from_dataframe(peaks[['chrom', 'start', 'end']]).sort()
closest = peaks_bt.closest(bt.BedTool(TSS).sort(), d=True).to_dataframe()
closest['peak'] = closest.chrom + ':' + closest.start.astype(str) + ':' + closest.end.astype(str)
peak_to_distal_proximal = {peak: 'TSS distal' if abs(distance) > 5000 else 'TSS proximal' for peak, distance in zip(closest.peak, closest.blockCount)}


motifs_bt = bt.BedTool().from_dataframe(motifs).sort()
peaks_bt = bt.BedTool().from_dataframe(peak_centers[['chrom', 'start', 'end']]).sort()
intersect = motifs_bt.intersect(peaks_bt, wa=True, wb=True).to_dataframe()
intersect['relative_motif_start'] = intersect['start'] - intersect['score']
intersect['relative_motif_end'] = intersect['relative_motif_start'] + MOTIF_LENGTH + 1
intersect['peak'] = intersect.name.astype(str) + ':' + intersect.score.astype(str) + ':' + intersect.strand.astype(str)
intersect.relative_motif_start = intersect.relative_motif_start.map(lambda x: 0 if x < 0 else x)
intersect.relative_motif_end = intersect.relative_motif_end.map(lambda x: PEAK_CENTER_FLANK*2 if x > PEAK_CENTER_FLANK*2 else x)

# convert to a matrix
# each row is a peak
# each col is a basepair
# each value indicates whether a motif overlaps that bp of that peak
motif_overlap = np.zeros((len(peak_rank), PEAK_CENTER_FLANK*2))
for peak, relative_motif_start, relative_motif_end in zip(intersect.peak, intersect.relative_motif_start, intersect.relative_motif_end):
    motif_overlap[peak_rank[new_peak_to_old_peak[peak]],relative_motif_start:(relative_motif_end+1)] = 1


PEAKS_PER_BIN = round(len(motif_overlap) / BINS)



df = pd.DataFrame(motif_overlap)
dfbin = df.groupby(lambda x: divmod(x, PEAKS_PER_BIN)[0]).sum()
correlation = dfbin[PEAK_CENTER_FLANK].reset_index().corr().iloc[0,1]
pd.DataFrame([[MOTIF, correlation]]).to_csv(f'/nfs/turbo/umms-scjp/kgmoo/analysis/20250428_kgmoo_ATACge/results/peterPlots/{PREFIX}.correlation.txt', index=False, header=False, sep='\t')

df = pd.DataFrame(motif_overlap, index=deseq2.index.to_list()) # TODO: peak_rank must be sorted by value
df['distal_proximal'] = df.index.map(peak_to_distal_proximal)
df['fc'] = df.index.map(fcs)

distal = df[df.distal_proximal=='TSS distal'].drop(columns=['distal_proximal'])
proximal = df[df.distal_proximal=='TSS proximal'].drop(columns=['distal_proximal'])



fig, axs = plt.subplots(ncols=2, nrows=2, figsize=(6, 10), gridspec_kw={'width_ratios': [1, 2], 'height_ratios': [len(distal)/len(proximal), 1], 'wspace':0.1})
fig.suptitle(f'{MOTIF} motif density')

# Distal first
df = distal.copy()
df['bin'] = [divmod(x, PEAKS_PER_BIN)[0] for x in range(len(df))]
dfbin = df.drop(columns=['fc']).groupby('bin').sum()
bin_fc = df.groupby('bin').fc.median().sort_index(ascending=False)

ax = axs[0,0]
ax.imshow(dfbin.values, aspect='auto')
ax.set_xticks([])
bottom, top = ax.get_ylim()
ax.set_yticks([bottom, top])
ax.set_yticklabels([f'', f''])
print(f'peaks per bin: {PEAKS_PER_BIN}')
ax.set_ylabel(f'{len(dfbin)} bins')
ax.set_title('TSS distal')

ax = axs[0,1]
ax.plot(bin_fc.to_list(), bin_fc.index)
ax.axvline(0, linestyle='dashed', color='black')
ax.set_xlabel('')
ax.set_ylim(bottom=axs[0,0].get_ylim()[0], top=axs[0,0].get_ylim()[1])
ax.set_xlim(left=-1*bin_fc.abs().max(), right=bin_fc.abs().max())
ax.set_yticks([])

# Proximal second
df = proximal.copy()
df['bin'] = [divmod(x, PEAKS_PER_BIN)[0] for x in range(len(df))]
dfbin = df.drop(columns=['fc']).groupby('bin').sum()
bin_fc = df.groupby('bin').fc.median().sort_index(ascending=False)

ax = axs[1,0]
ax.imshow(dfbin.values, aspect='auto')
ax.set_xticks([0, PEAK_CENTER_FLANK, PEAK_CENTER_FLANK*2])
ax.set_xticklabels([-1*PEAK_CENTER_FLANK, 0, PEAK_CENTER_FLANK])
bottom, top = ax.get_ylim()
ax.set_yticks([bottom, top])
ax.set_yticklabels([f'', f''])
ax.set_ylabel(f'{len(dfbin)} bins')
ax.set_title('TSS proximal')

ax = axs[1,1]
ax.plot(bin_fc.to_list(), bin_fc.index)
ax.axvline(0, linestyle='dashed', color='black')
ax.set_xlabel(f'Z-score for OCR association with age')
ax.set_ylim(bottom=axs[1,0].get_ylim()[0], top=axs[1,0].get_ylim()[1])
ax.set_xlim(left=-1*bin_fc.abs().max(), right=bin_fc.abs().max())
ax.set_yticks([])

fig.savefig(f'/nfs/turbo/umms-scjp/kgmoo/analysis/20250428_kgmoo_ATACge/results/peterPlots/{PREFIX}.png', bbox_inches='tight')
