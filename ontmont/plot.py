from collections import defaultdict

import pysam
import numpy as np
import pandas as pd
import matplotlib
from matplotlib.path import Path
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from matplotlib.patches import ConnectionPatch
from scipy import interpolate

from .datatypes import Breakpoint
from .collect import get_svtype, normalize_sv_table


##################
##  Plot genes  ##
##################

def get_edgecolor(svtype):
    edgecolor_map = {
        'TRA': 'grey',
        'DEL': 'tab:blue',
        'DUP': 'red',
        'INS': 'green',
        'INV': 'purple',
    }
    return edgecolor_map[svtype]

def is_overlap(interval1, interval2):
    start1, end1 = interval1
    start2, end2 = interval2
    start1, end1, start2, end2 = float(start1), float(end1), float(start2), float(end2)
    overlap = (start2 < end1) and (start1 < end2)
    return overlap

def convert_tupleproxy_to_pyranges(exon):
    assert len(exon) == 9, f'exon (length:{len(exon)}): {tuple(exon)}'
    exon_field = exon[:8]
    info = exon[8].rstrip(';')
    annots = info.split('; ')
    annot_map = {v.split(' ')[0]: v.split(' ')[1][1:-1] for v in annots} # [1:-1] hack to quick-remove ""

    annot_cols = ['gene_id', 'gene_version', 'gene_name', 'gene_source',
        'gene_biotype', 'transcript_id', 'transcript_version',
        'transcript_name', 'transcript_source', 'transcript_biotype', 'tag',
        'transcript_support_level', 'exon_number', 'exon_id', 'exon_version',
        'protein_id', 'protein_version', 'ccds_id']

    annot_field = [annot_map[k] if k in annot_map else np.nan for k in annot_cols]
    field = exon_field + annot_field
    return field

def parse_gtf_region(gtf, region):
    pyranges_cols = [
        'Chromosome', 'Source', 'Feature', 'Start', 'End', 'Score', 'Strand',
        'Frame', 'gene_id', 'gene_version', 'gene_name', 'gene_source',
        'gene_biotype', 'transcript_id', 'transcript_version',
        'transcript_name', 'transcript_source', 'transcript_biotype', 'tag',
        'transcript_support_level', 'exon_number', 'exon_id', 'exon_version',
        'protein_id', 'protein_version', 'ccds_id'
    ]
    int_cols = ['Start', 'End',]
    df = pd.DataFrame(columns=pyranges_cols)
    for exon in gtf.fetch(*region, parser=pysam.asTuple()):
        field = convert_tupleproxy_to_pyranges(exon)
        df.loc[df.shape[0]] = field
    df[int_cols] = df[int_cols].astype(int)
    return df

def get_gene_repr_exons(gtf, gene, lenient=False):
    transcript = get_repr_transcript_id(gtf, gene, lenient=lenient)
    exons = get_transcript_exons(gtf, transcript)
    return exons

def get_transcript_exons(gtf, transcript_id):
    exons = gtf[
        (gtf['transcript_id']==transcript_id) &
        (gtf['Feature'].isin(['exon', 'CDS']))
    ]
    return exons

def get_repr_transcript_id(gtf, gene_name, lenient=False):
    if not lenient:
        transcripts = (
            gtf
            .query(f"gene_name == '{gene_name}'")
            .query("Source == 'protein_coding'")
            # .query("Feature == 'transcript'")
            # .query("transcript_biotype == 'protein_coding'")
            # .query("transcript_support_level == '1'")
        ).copy()
        transcripts = transcripts[
            ~transcripts['gene_name'].str.startswith('OR')
        ]
    else:
        transcripts = (
            gtf
            .query(f"gene_name == '{gene_name}'")
            # .query("Feature == 'transcript'")
        ).copy()
    transcript_id = None
    if transcripts.shape[0] > 0:
        transcripts['length'] = transcripts['End'] - transcripts['Start']
        transcript = transcripts.sort_values(by=['length'], ascending=False).iloc[0]
        transcript_id = transcript['transcript_id']
    return transcript_id

def _get_gene_intervals(gene_exons):
    # gather gene start-end intervals
    intervals = []
    for gene, exons in gene_exons.items():
        start_min = exons['Start'].min()
        end_max = exons['End'].max()
        if not np.isnan(start_min) and not np.isnan(end_max):
            interval = (start_min, end_max, gene)
            intervals.append(interval)
    intervals.sort(key=lambda x: (x[0], x[1]))
    return intervals

def _get_gene_offsets(intervals):
    # set gene offset according to overlapping intervals
    gene2offset = defaultdict(int)
    n_intervals = len(intervals)
    margin = 5e7
    for i in range(n_intervals):
        for j in range(i+1, n_intervals):
            start_i, end_i, gene_i = intervals[i]
            start_j, end_j, gene_j = intervals[j]
            if is_overlap((start_i-margin, end_i+margin), (start_j-margin, end_j+margin)):
                if gene2offset[gene_i] == gene2offset[gene_j]:
                    gene2offset[gene_j] += 1
    return gene2offset

def add_gene_annotations(ax, gtf, chromosome, start, end, genes=None, already_plotted=set(), gene_font_size=8):
    chromosome = chromosome.replace('chr', '')
    # ax.set_xlim((start, end))
    # ax.spines[['top', 'bottom', 'right', 'left']].set_visible(False)
    # ax.set_xticks([]); ax.set_yticks([])

    # setup which genes to plot
    region = (chromosome, start, end)
    gene_df = parse_gtf_region(gtf, region)
    if genes is None: 
        genes = gene_df['gene_name'].unique() # if no genes put in, use all genes in region
    gene_exons = {gene: get_gene_repr_exons(gene_df, gene, lenient=False)
                  for gene in genes}
    intervals = _get_gene_intervals(gene_exons)
    gene2offset = _get_gene_offsets(intervals)
    max_offset = 1
    if len(gene2offset):
        max_offset = max(gene2offset.values())
    
    # apply offset and plot exon rectangles
    for gene, exons in gene_exons.items():
        if gene in already_plotted:
            continue
        already_plotted.add(gene)
        if exons.shape[0] == 0:
            continue # filter out if none
        strands = exons['Strand'].unique()
        assert len(strands) == 1, strands
        strand = strands[0]
        transcript_s, transcript_e = exons['Start'].min(), exons['End'].max()
        offset = gene2offset[gene]
        for i, (_, row) in enumerate(exons.iterrows()):
            add_annot = (i==0)
            _add_transcript_rectangles(ax, row, transcript_s, transcript_e, offset, strand, add_annot=add_annot, gene_font_size=gene_font_size) # add rectangles and a line

    return already_plotted, max_offset

def _add_quiver(ax, x_quiver_starts, y_quiver_starts, x_quiver_directs, y_quiver_directs, box_color='#090A76'):
    ax.quiver(
        x_quiver_starts, y_quiver_starts, x_quiver_directs, y_quiver_directs,
        color=box_color,
        width=0.005, 
        headwidth=30, 
        headaxislength=5,
        headlength=17, 
        pivot='tip',
        scale=100,
    )

def _add_transcript_rectangles(ax, row, transcript_s, transcript_e, y_offset, strand, add_annot=False, gene_font_size=8, zorder=-1):
    box_color = '#090A76'
    strand_color = 'blue'
    x_direct, y_direct = -1, 0
    if strand == '+': 
        x_direct = 1
        strand_color = 'red'
    x_min, x_max = ax.get_xlim()
    y_offset_factor = 3
    y_offset *= y_offset_factor
    exon_s, exon_e = row['Start'], row['End']
    feature = row['Feature']
    rect_y, rect_height = 1+y_offset, 1
    transcript_line_y = rect_y + 0.5 * rect_height
    text_y = rect_y - rect_height + 3
    if feature == 'CDS':
        rect_y = 0.5 + y_offset
        rect_height = 2
    exon_len = exon_e - exon_s
    rect_x = exon_s #max(exon_s, start)
    
    # for each exon / CDS
    rectangle = matplotlib.patches.Rectangle((rect_x, rect_y), exon_len, rect_height,
            linewidth=1, edgecolor=strand_color, color=box_color, zorder=zorder)
    ax.add_patch(rectangle)
    if add_annot: 
        quiver_interval = int((x_max - x_min) / 50)
        if quiver_interval > 0:
            x_quiver_starts = np.arange(transcript_s, transcript_e, quiver_interval)[1:-1]
            n_quivers = len(x_quiver_starts)
            y_quiver_starts = [transcript_line_y] * n_quivers
            x_quiver_directs = [x_direct] * n_quivers
            y_quiver_directs = [y_direct] * n_quivers
            _add_quiver(ax, x_quiver_starts, y_quiver_starts, x_quiver_directs, y_quiver_directs)
            ax.plot([transcript_s+50, transcript_e-100], 
                    [transcript_line_y, transcript_line_y], 
                    linewidth=1, color=box_color, zorder=0)
        text_x = (transcript_s + transcript_e) / 2
        if text_x < x_min:
            text_x = x_min
        if text_x > x_max:
            text_x = x_max
        ax.text(x=text_x, y=text_y, s=row['gene_name'], color=box_color, horizontalalignment='center', fontdict={'size':gene_font_size}, zorder=zorder)


#######################
##  Plot quasijabba  ##
#######################
def plot_quasijabba(input_svs, cn, cn_clones=['Pseudobulk'], 
                    figsize=(12,4), svlineh=0.8, width_ratios=None, alpha=0.3, 
                    loc1=(0.83, 0.75), default_rad=0.25, suptitle='', read_length=None, gtf=None, gene_margin=50000,
                    bin_size=50000, linewidth=0, show_support_legend=True, gene_font_size=8, sv_x_length_margin=0.2,
                    intra_y_offset=1, linetick_divfactor=3, ylim=(0, 10)):
    # alpha = 0.3
    read_chroms = get_sv_table_chroms(input_svs)
    chrom_xlims = get_chrom_xlims(input_svs, read_chroms, sv_x_length_margin)
    chrom_cns = get_cn_blocks_for_read_xlim(cn, read_chroms, chrom_xlims, bin_size=bin_size, clones=cn_clones)
    cn_colors = {2:'#CCCCCC', 0: 'tab:blue', 1:'lightblue', 3:'#FA8072', 4:'#FF2400', 5:'#C21807', 6:'#800000'}
    cn_colors.update({i:'#5E1914' for i in range(7, 50)})
    
    svcolors = {}
    with matplotlib.rc_context({'font.family':'Arial'}):
        # fig, axes = format_pseudojabba_axes(figsize, read_chroms, width_ratios)
        # fig, mode_axes = make_axes_for_breakpoints(read_chroms, figsize=figsize, rowsizes=(2, 6), modes=['vaf', 'sv'], width_ratios=width_ratios)
        fig, mode_axes = make_axes_for_breakpoints(read_chroms, figsize=figsize, rowsizes=(5, 2), modes=['sv', 'gene'], width_ratios=width_ratios, hspace=0)
        if suptitle: 
            fig.suptitle(suptitle, y=1.07)
        if read_length:
            read_length_info = f'(read length: {int(round(read_length / 1000, 0))} kb)'
            fig.text(0.5, 1, read_length_info, fontsize=10, ha='center')
        # chrom_axes = get_chrom_axes(axes, read_chroms, chrom_xlims)
        axes = mode_axes['sv']
        gene_axes = mode_axes['gene']

        for chrom in read_chroms:
            ax = axes[chrom]
            chrom_cn = chrom_cns[chrom]
            plot_jabba_cn(chrom_cn, ax, cn_colors, ylim=ylim)
    
        svcolors = plot_svs_to_cn_segments(input_svs, axes, chrom_cns, svlineh, alpha, default_rad=default_rad, sv_x_length_margin=sv_x_length_margin,
                                           _linewidth=linewidth, intra_y_offset=intra_y_offset, linetick_divfactor=linetick_divfactor)

        if 'support' in input_svs.columns:
            add_sv_legend_to_axes(axes, svcolors, input_svs['support'].unique(), loc1=loc1, loc2=(.83, .5), alignment='left', size=8, show_support_legend=show_support_legend)
        else:
            add_sv_legend_to_axes(axes, svcolors, [], loc1=loc1, loc2=(.83, .5), alignment='left', size=8, show_support_legend=show_support_legend)
        # print(svcolors)

        if type(gtf) == pysam.libctabix.TabixFile:
            for chrom in read_chroms:
                ax = gene_axes[chrom]
                plot_gene_annotations(ax, gtf, chrom, input_svs, margin=gene_margin, gene_font_size=gene_font_size)

        # plt.tight_layout()
        fix_coordinate_xticklabels(gene_axes)
        return fig

def plot_gene_annotations(ax, gtf, chrom, input_svs, margin=50000, gene_font_size=8):
    max_offset = -1
    genes_already_plotted=set()
    for rix, row in input_svs.iterrows():
        chrom1 = row['chromosome_1']
        chrom2 = row['chromosome_2']
        pos1 = row['position_1']
        pos2 = row['position_2']
        start1, end1 = pos1 - margin, pos1 + margin
        start2, end2 = pos2 - margin, pos2 + margin
        if chrom1 == chrom:
            genes_already_plotted, _max_offset = add_gene_annotations(ax, gtf, chrom1, start1, end1, already_plotted=genes_already_plotted, gene_font_size=gene_font_size)
            if max_offset < _max_offset:
                max_offset = _max_offset
        if chrom2 == chrom:
            genes_already_plotted, _max_offset = add_gene_annotations(ax, gtf, chrom2, start2, end2, already_plotted=genes_already_plotted, gene_font_size=gene_font_size)
            if max_offset < _max_offset:
                max_offset = _max_offset
    ax.set_ylim((0, 8 * max_offset))

def _extract_cn_matching_sv(chrom_cn, pos, default_ploidy=2):
    sv_cn = chrom_cn[(chrom_cn['start'] <= pos) & (chrom_cn['end'] > pos)]
    if sv_cn.shape[0] > 0:
        sv_cn = sv_cn.squeeze()['state']
    else:
        sv_cn = default_ploidy
    return sv_cn

def plot_svs_to_cn_segments(input_svs, axes, chrom_cns, svlineh, alpha, _linewidth=0, intra_y_offset=1, linetick_divfactor=3, default_rad=0.25, sv_x_length_margin=0.2, sv_x_offset=100000):
    svcolors = {}
    tra_color = 'black'
    chrom_sv_cnt = defaultdict(int)
    chrom_xlims = pd.DataFrame(columns=['chrom', 'pos'])
    for rix, row in input_svs.iterrows():
        chrom1, pos1, ori1, chrom2, pos2, ori2 = row.iloc[:6]
        chrom_cn1 = chrom_cns[chrom1]
        chrom_cn2 = chrom_cns[chrom2]
        chrom_xlims.loc[chrom_xlims.shape[0]] = [chrom1, pos1]
        chrom_xlims.loc[chrom_xlims.shape[0]] = [chrom2, pos2]
        sv_cn1 = _extract_cn_matching_sv(chrom_cn1, pos1)
        sv_cn2 = _extract_cn_matching_sv(chrom_cn2, pos2)
        svtype = row['type']
        svcolor = get_edgecolor(svtype)
        svcolors[svtype] = svcolor
        support = 1
        in_source = 0
        if 'in_source' in row.index: in_source = int(row['in_source'])
        if 'support' in row.index: support = row['support']
        if _linewidth == 0: linewidth = 0.5+np.log2(support)
        linewidth = _linewidth + 1 + in_source
        # print(linewidth, svtype)
        
        if chrom1 != chrom2: # TRA
            yoffset1 = svlineh/2 if ori1 == '+' else -svlineh/2
            yoffset2 = svlineh/2 if ori2 == '+' else -svlineh/2
            rad = -default_rad if ori1 == '+' else default_rad
            assert chrom1 in axes and chrom2 in axes
            ax1, ax2 = axes[chrom1], axes[chrom2]
            ax1.plot([pos1, pos1], [sv_cn1-svlineh/2, sv_cn1+svlineh/2], color=tra_color, linestyle='--', linewidth=0.5, solid_capstyle='butt')
            ax2.plot([pos2, pos2], [sv_cn2-svlineh/2, sv_cn2+svlineh/2], color=tra_color, linestyle='--', linewidth=0.5, solid_capstyle='butt')
            sv_patch = ConnectionPatch(xyA=[pos1, sv_cn1+yoffset1], xyB=[pos2, sv_cn2+yoffset2], 
                                       coordsA='data', coordsB='data',
                                       axesA=ax1, axesB=ax2, arrowstyle='-', capstyle='butt', alpha=alpha,
                                       color=tra_color, lw=linewidth, connectionstyle=f"arc3,rad={rad}")
            fig = ax1.get_figure()
            fig.add_artist(sv_patch)

        else: # DEL DUP INV
            ax = axes[chrom1]
            x_cnt_offset = chrom_sv_cnt[svtype+chrom1]
            y_cnt_offset = chrom_sv_cnt[svtype+chrom1] * 0.3
            chrom_sv_cnt[svtype+chrom1] += 1
            ax.plot([pos1, pos1], [sv_cn1-svlineh/2, sv_cn1+svlineh/2], color=svcolor, linestyle='--', linewidth=0.5, solid_capstyle='butt')
            ax.plot([pos2, pos2], [sv_cn2-svlineh/2, sv_cn2+svlineh/2], color=svcolor, linestyle='--', linewidth=0.5, solid_capstyle='butt')
            linetick = (pos2-pos1) / linetick_divfactor * (1 * (x_cnt_offset))
            secondx = pos1 + linetick if ori1 == '+' else pos1 - linetick
            thirdx = pos2 + linetick if ori2 == '+' else pos2 - linetick
            secondy = sv_cn1 + y_cnt_offset
            thirdy = sv_cn2 + y_cnt_offset
            if sv_cn1 == sv_cn2:
                secondy = sv_cn1 + y_cnt_offset + intra_y_offset
                thirdy = sv_cn2 + y_cnt_offset + intra_y_offset
                
            pp = matplotlib.patches.PathPatch(
                Path([(pos1, sv_cn1), (secondx, secondy), (thirdx, thirdy), (pos2, sv_cn2)],
                     [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4]),
                fc="none", transform=ax.transData, color=svcolor, alpha=alpha, lw=linewidth)
            ax.add_patch(pp)
    for chrom, chromdf in chrom_xlims.groupby('chrom'):
        pos_min = chromdf['pos'].min()
        pos_max = chromdf['pos'].max()
        length = pos_max - pos_min
        ax = axes[chrom]
        xmin = max(0, pos_min - length * sv_x_length_margin)
        xmax = pos_max + length * sv_x_length_margin
        x_length = xmax-xmin
        if x_length < 2 * sv_x_offset:
            xmin = max(0, pos_min - sv_x_offset)
            xmax = pos_max + sv_x_offset
        xlim = (xmin, xmax)
        ax.set_xlim(xlim)
    return svcolors

def get_cn_blocks_for_read_xlim(cn, read_chroms, chrom_xlims, bin_size=50000, clones=['Pseudobulk']):
    chrom_cns = {}
    for chrom in read_chroms:
        sv_min, sv_max = chrom_xlims[chrom]
        x_range = sv_max - sv_min
        if x_range < 5*bin_size:
            cn_min, cn_max = sv_min - bin_size, sv_max + bin_size
        else:
            cn_min = sv_min
            cn_max = sv_max
        chrom_xlims[chrom] = (cn_min, cn_max)
        chrom_cn = cn[
            (cn['chr'] == chrom) &
            (cn['start'] < cn_max) &
            (cn['end'] > cn_min)
        ]
        # chrom_cn = get_blocks_from_cn(chrom_cn)
        if len(clones):
            chrom_cn = chrom_cn[chrom_cn['clone_id'].isin(clones)]
        chrom_cns[chrom] = chrom_cn#.reset_index(drop=True)
    return chrom_cns

def get_chrom_xlims(read_svs, read_chroms, sv_x_length_margin):
    chrom_xlims = {}
    for chrom in read_chroms:
        brk1s = read_svs[['chromosome_1', 'position_1']]
        brk1s = brk1s[brk1s['chromosome_1']==chrom]
        pos1s = brk1s['position_1'].tolist()
        brk2s = read_svs[['chromosome_2', 'position_2']]
        brk2s = brk2s[brk2s['chromosome_2']==chrom]
        pos2s = brk2s['position_2'].tolist()
        positions = pos1s + pos2s
        min_pos = np.min(positions)
        max_pos = np.max(positions)
        length = max_pos - min_pos
        xmin = max(0, min_pos - length * sv_x_length_margin)
        xmax = min_pos + length * sv_x_length_margin
        chrom_xlims[chrom] = (xmin, xmax)
    return chrom_xlims

def get_sv_table_chroms(sv):
    chroms = sorted(set(sv['chromosome_1'].unique().tolist() + sv['chromosome_2'].unique().tolist()), 
                    key=lambda x: Breakpoint.chroms.index(x))
    return chroms

def plot_jabba_cn(chrom_cn, ax, cn_colors, ylim=(0, 6)):
    for rix, row in chrom_cn.iterrows():
        start = row['start']
        end = row['end']
        width = end-start
        height = 0.5
        cn_state = row['state']
        cn_color = cn_colors[cn_state]
        rectangle = matplotlib.patches.Rectangle(
            xy=(start, cn_state-height/2), width=width, height=height,
            capstyle='butt', facecolor=cn_color, edgecolor='black', linewidth=0.5,
        )
        ax.add_patch(rectangle)
        ax.set_ylim(ylim)
        # ax.plot([start, end], [cn_state, cn_state], linewidth=10, color='grey', solid_capstyle='butt')


######################
##  Plot VAF-CN-SV  ##
######################
def plot_vaf_cn_sv_for_reads(plot_data, sv, cn, clone_ids,
                             fig_row_sizes=(15, 40, 50), cn_window=200000, cn_metric='state'):
    colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:pink', 'tab:purple', 'tab:olive', 'tab:brown']
    brk_chroms = get_chromosomes_from_bundle(plot_data)
    clone_colors = dict(zip(clone_ids, colors[:len(clone_ids)]))
    
    fig, axes = make_axes_for_breakpoints(brk_chroms, 
        figsize=(4*len(brk_chroms), 4), rowsizes=fig_row_sizes)
    
    plot_segments_on_axes(plot_data, axes['sv'])
    plot_transition_on_axes(plot_data, axes['sv'])
    plot_vaf_on_axes(plot_data, axes['vaf'], clone_colors, sv)
    plot_cn_on_axes(cn, axes['cn'], clone_colors, bin_size=cn_window, metric=cn_metric)
    plot_sv_lines_on_cn(plot_data, axes['cn'])
    fix_coordinate_xticklabels(axes['sv'])
    return fig

def get_chromosomes_from_bundle(bundle):
    brk_chroms = set()
    for brks in bundle:
        for brk in brks:
            _chrom = brk.chrom
            brk_chroms.add(_chrom)
    brk_chroms = list(sorted(brk_chroms, key=lambda x: Breakpoint.chroms.index(x)))
    return brk_chroms

def make_axes_for_breakpoints(brk_chroms, figsize=(10, 4), hspace=1.5, rowsizes=(1, 4, 5, 1), modes=['vaf', 'cn', 'sv', 'gene'], width_ratios=None):
    assert len(rowsizes) == len(modes), (rowsizes, modes)
    fig = plt.figure(figsize=figsize)
    y_intervals = [(sum(rowsizes[:i]), sum(rowsizes[:i+1])) for i in range(len(rowsizes))]
    gs = {}
    gs['main'] = gridspec.GridSpec(1, len(brk_chroms), figure=fig, width_ratios=width_ratios)
        # left=0.1, right=0.9, top=5, bottom=1, 
        # wspace=10, hspace=20)
    axes = {mode:{} for mode in modes}
    mode_titles = {'vaf':'DLP\nVAF', 'cn':'DLP\nCN', 'sv':'ONT\nSVs', '_gap1':'', 'gene':'Genes'}
    mode_labelpads = {'vaf':15, 'cn':15, 'sv':15, 'gene':50}
    if len(modes) == 2 and 'cn' not in modes and 'vaf' not in modes: #
        mode_titles['sv'] = 'DLP CN\nw/ ONT SVs'
        mode_labelpads['sv'] = 30
    with matplotlib.rc_context({'font.family':'Arial'}):
        for ix, chrom in enumerate(brk_chroms): # x
            gs[chrom] = gridspec.GridSpecFromSubplotSpec(sum(rowsizes), 1, subplot_spec=gs['main'][ix])

            mode_axes = []
            for jx, mode in enumerate(modes):
                labelpad = mode_labelpads[mode]
                mode_title = mode_titles[mode]
                ax = fig.add_subplot(gs[chrom][y_intervals[jx][0]:y_intervals[jx][1], :], facecolor='white')
                if jx == 0: # mode axis index
                    ax.set_title('chr'+chrom.replace('chr', ''))
                if ix == 0: # chrom axis index
                    ax.set_ylabel(mode_title, rotation=0, labelpad=labelpad, va='center')
                if mode == 'vaf':
                    ax.set_ylim((0, 1))
                if mode == 'sv':
                    ax.set_zorder(0)
                    if ix != 0: # if not first chrom index
                        ax.set_yticks([])
                        ax.set_yticklabels([])
                        ax.spines[['left']].set_visible(False)
                    elif len(modes) > 1 and set(modes) != {'sv', 'gene'}:
                        ax.set_yticks([])
                        ax.set_yticklabels([])
                        ax.spines[['left']].set_visible(False)
                if mode == 'gene':
                    ax.set_yticks([])
                    ax.set_yticklabels([])
                    ax.spines[['left']].set_visible(False)
    
                mode_axes.append(ax)
                
            for jx, mode in enumerate(modes): # hide xticks, share x for all but last m-ax
                ax = mode_axes[jx]
                last_mode_ax = mode_axes[-1]
                if jx != len(modes)-1:
                    ax.spines[['bottom']].set_visible(False)
                    ax.set_xticks([])
                    ax.tick_params('x', labelbottom=False, color=[1,1,1,1])
                    plt.setp(ax.get_xticklabels(), visible=False)
                    ax.sharex(last_mode_ax)

                ax.spines[['right', 'top']].set_visible(False)
                ax.autoscale(enable=False, axis='y')
                axes[mode][chrom] = ax
        fig.subplots_adjust(hspace=hspace)

        return fig, axes

def fix_coordinate_xticklabels(sv_axes):
    for chrom, ax in sv_axes.items():
        xticklabels = ax.get_xticklabels()
        new_xticklabels = []
        for xticklabel in xticklabels:
            text = str(int(xticklabel.get_position()[0]))
            xticklabel.set_text(text)
            xticklabel.set_c('black')
            new_xticklabels.append(xticklabel)
        ax.set_xticklabels(new_xticklabels, rotation=90)
        # for each transition, check DLP overlap, and plot vaf as bars

def plot_sv_lines_on_cn(plot_data, cn_axes):
    sv_df = get_unique_sv_with_support(plot_data)
    for rix, row in sv_df.iterrows():
        chrom1, pos1, ori1 = row['chromosome_1'], row['position_1'], row['strand_1']
        chrom2, pos2, ori2 = row['chromosome_2'], row['position_2'], row['strand_2']
        svtype = row['type']
        svcolor = get_edgecolor(svtype)
        support = row['support']
        ax1, ax2 = cn_axes[chrom1], cn_axes[chrom2]
        ax1.plot([pos1, pos1], [-1, 50], color=svcolor, linestyle='--', linewidth=0.5)
        ax2.plot([pos2, pos2], [-1, 50], color=svcolor, linestyle='--', linewidth=0.5)

def fix_yticks_with_integers(cn_axes):
    for chrom, ax in cn_axes.items():
        ymin, ymax = ax.get_ylim()
        ymax_round = int(round(ymax, 0))
        if ymax_round < 8:
            step = 1
        else:
            step = 2
        ax.set_yticks(np.arange(0, ymax_round+1, step))
        ax.set_yticklabels(np.arange(0, ymax_round+1, step))

def plot_cn_on_axes(cn, cn_axes, clone_colors, bin_size=50000, metric='copy', max_cn=6):
    min_xsize = 2 * bin_size
    clone_ids = list(sorted(clone_colors.keys()))
    _max_cn = max_cn
    axes_fixed = False
    for cix, clone_id in enumerate(clone_ids):
        clone_color = clone_colors[clone_id]
        ccn = cn[cn['clone_id']==clone_id]
        for chrom, ax in cn_axes.items():
            max_cn = _max_cn
            xmin, xmax = ax.get_xlim()
            xsize = xmax - xmin
            if cix == 0:
                if xsize < min_xsize:
                    diff = min_xsize - xsize
                    _xmin = xmin - diff / 2
                    _xmax = xmax + diff / 2
                else:
                    _xmin = xmin - xsize*0.15
                    _xmax = xmax + xsize*0.15
                ax.set_xlim((_xmin, _xmax))
                # print(chrom, _xmin, _xmax)
            else:
                _xmin, _xmax = xmin, xmax
            _cn = ccn[
                (ccn['chr']==chrom) & 
                (ccn['start'] < _xmax + bin_size) & 
                (ccn['end'] > _xmin - bin_size)
            ]
            if _cn.shape[0] > 0:
                new_max_cn = int(_cn[metric].max()) + 1
                if new_max_cn > max_cn:
                    max_cn = new_max_cn
            ax.set_ylim((0, max_cn))
            for rix, row in _cn.iterrows():
                start = row['start']
                end = row['end']
                state = row[metric]
                ax.plot([start, end], [state, state], color=clone_color, alpha=0.5)
    fix_yticks_with_integers(cn_axes)

def extract_vaf_from_annotated_sv_table(tra, sv, margin=10):
    brk1, brk2 = tra.brk1, tra.brk2
    vaf_cols = [c for c in sv.columns if c.startswith('vaf_')]
    clone_vafs = {}
    if brk1 > brk2:
        brk1, brk2 = brk2, brk1
    flt = sv[
        (sv['chromosome_1']==brk1.chrom) &
        (sv['chromosome_2']==brk2.chrom) &
        (sv['strand_1']==brk1.ori) &
        (sv['strand_2']==brk2.ori) &
        (np.abs(sv['position_1']-brk1.pos) < margin) &
        (np.abs(sv['position_2']-brk2.pos) < margin)
    ]
    if flt.shape[0] > 0:
        if flt.shape[0] > 1:
            flt = flt.iloc[:1]
        flt = flt.squeeze()
        for vaf_col in vaf_cols:
            vaf = flt[vaf_col]
            clone_id = vaf_col.replace('vaf_', '')
            clone_vafs[clone_id] = vaf
    return clone_vafs

def plot_vaf_on_axes(plot_data, vaf_axes, clone_colors, flt_ont_sv):
    with matplotlib.rc_context({'font.family':'Arial'}):
        for brks in plot_data:
            for tra in brks.tras:
                svtype = get_svtype(tra)
                svcolor = get_edgecolor(svtype)
                brk1, brk2 = tra.brk1, tra.brk2
                chrom1, pos1, ori1 = brk1.chrom, brk1.pos, brk1.ori
                chrom2, pos2, ori2 = brk2.chrom, brk2.pos, brk2.ori
                clone_vafs = extract_vaf_from_annotated_sv_table(tra, flt_ont_sv)
                vaf_y1 = 0
                vaf_y0 = 0
                for ix, (clone_id, clone_vaf) in enumerate(clone_vafs.items()):
                    clone_color = clone_colors[clone_id]
                    vaf_y1 += clone_vaf
                    zorder = -ix
                    vaf_axes[chrom1].plot([pos1, pos1], [vaf_y0, vaf_y1], color=clone_color, linewidth=2)
                    vaf_axes[chrom2].plot([pos2, pos2], [vaf_y0, vaf_y1], color=clone_color, linewidth=2)
                    vaf_y0 = vaf_y1
        brk_chroms = sorted(vaf_axes.keys(), key=lambda x: Breakpoint.chroms.index(x))
        end_chrom = brk_chroms[-1]
        ncols = 1 if len(clone_colors.keys()) <= 4 else 2
        alignment = 'left' if ncols == 1 else 'center'
        # print(ncols, clone_colors.keys(), clone_colors)
        handles = [plt.plot([], [],
                   color=clone_colors[label], marker="s", ms=4, ls="")[0] for label in clone_colors]
        # print(handles)
        legend = vaf_axes[end_chrom].legend(handles=handles, labels=clone_colors.keys(), title="clones", handletextpad=0.1, alignment=alignment, columnspacing=0.2,
                                            frameon=False, loc='upper left', bbox_to_anchor=(0.82-(ncols-1)*0.1,1.5), prop={'size': 8}, ncols=ncols)
        # print(legend)
        vaf_axes[end_chrom].add_artist(legend)
        # sns.move_legend(ax, loc='upper left', bbox_to_anchor=(1, 1))

def add_sv_legend_to_axes(sv_axes, svcolors, supports=None, alpha=0.5, loc1=(0.82, 0.5), loc2=(0.82, 0.0), size=8, alignment='left', show_support_legend=True):
    if 'TRA' in svcolors:
        svcolors['TRA'] = 'black'
    with matplotlib.rc_context({'font.family':'Arial'}):
        brk_chroms = list(sorted(sv_axes.keys(), key=lambda x: Breakpoint.chroms.index(x)))
        end_chrom = brk_chroms[-1]
        ax = sv_axes[end_chrom]
        
        handle1 = [plt.plot([], [], alpha=alpha,
            color=svcolors[label], marker="s", ms=4, ls="")[0] for label in svcolors]
        legend1 = ax.legend(handles=handle1, labels=svcolors.keys(), title="types", loc=loc1, prop={'size':size}, frameon=False, handletextpad=0.1, alignment=alignment)
        ax.add_artist(legend1);

        if type(supports) == type(None):
            supports = [int(x) for x in sorted(supports)]
            min_support, max_support = min(supports), max(supports)
            sizes = {min_support:min_support, max_support:max_support}
            if len(supports) >= 5:
                med_support = int(np.median(supports))
                sizes = {min_support:min_support, med_support:med_support, max_support:max_support}
            if show_support_legend:
                handle2 = [plt.plot([], [], alpha=alpha,
                           color="black", marker="s", ms=np.sqrt(sizes[s]/2), ls="")[0] for s in sizes]
                legend2 = ax.legend(handles=handle2, labels=sizes.keys(), title="supports", loc=loc2, prop={'size':size}, frameon=False, handletextpad=0.1, alignment=alignment)
                ax.add_artist(legend2);

def make_spline_coordinates(pos1, pos2, sv_y=1, margin_div=5, y_offset=0.05):
    dist = (pos2-pos1)
    margin_x = dist / margin_div
    sv_xs = [
        pos1, 
        pos1-margin_x, 
        # pos1+dist*2/4, 
        pos2+margin_x, 
        pos2,
    ]
    sv_ys = [
        sv_y+sv_y/2, 
        sv_y+sv_y/2+y_offset, 
        # sv_y+sv_y/2+0.15, 
        sv_y+sv_y/2+y_offset, 
        sv_y+sv_y/2, 
    ]
    try:
        tck, u = interpolate.splprep([sv_xs, sv_ys], s=0, k=2)
        unew = np.arange(0, 1.01, 0.01)
        out = interpolate.splev(unew, tck)
        new_xs = out[0]
        new_ys = out[1]
    except ValueError:
        new_xs = sv_xs
        new_ys = sv_ys
    return (new_xs, new_ys)

def get_unique_sv_with_support(bundle):
    sv_cols = ['chromosome_1', 'position_1', 'strand_1', 'chromosome_2', 'position_2', 'strand_2', 'type']
    sv_df = pd.DataFrame(columns=sv_cols)
    for brks in bundle:
        for tra in brks.tras:
            svtype = get_svtype(tra)
            brk1, brk2 = tra.brk1, tra.brk2
            chrom1, pos1, ori1 = brk1.chrom, brk1.pos, brk1.ori
            chrom2, pos2, ori2 = brk2.chrom, brk2.pos, brk2.ori
            field = [chrom1, pos1, ori1, chrom2, pos2, ori2, svtype]
            sv_df.loc[sv_df.shape[0]] = field
    sv_df = normalize_sv_table(sv_df)
    sv_df = pd.DataFrame(sv_df.value_counts()).reset_index()
    sv_df.rename(columns={'count':'support'}, inplace=True)
    return sv_df

def plot_transition_on_axes(plot_data, sv_axes, alpha=0.5):
    tra_color = 'black'
    sv_df = get_unique_sv_with_support(plot_data)
    chrompair_y = defaultdict(float)
    svcolors = {}
    supports = set()
    for rix, row in sv_df.iterrows():
        chrom1, pos1, ori1 = row['chromosome_1'], row['position_1'], row['strand_1']
        chrom2, pos2, ori2 = row['chromosome_2'], row['position_2'], row['strand_2']
        svtype = row['type']
        svcolor = get_edgecolor(svtype)
        support = row['support']
        svcolors[svtype] = svcolor
        supports.add(support)
        chrompair_y[(chrom1, chrom2)] += 1
        sv_y = 1 
        sv_rad = chrompair_y[(chrom1, chrom2)] + 0.2
        if chrom1 != chrom2:
            ax1, ax2 = sv_axes[chrom1], sv_axes[chrom2]
            ax1.plot([pos1, pos1], [sv_y-sv_y/2, sv_y+sv_y/2], color=tra_color, linestyle='--', linewidth=0.5, solid_capstyle='butt')
            ax2.plot([pos2, pos2], [sv_y-sv_y/2, sv_y+sv_y/2], color=tra_color, linestyle='--', linewidth=0.5, solid_capstyle='butt')
            sv_patch = ConnectionPatch(xyA=[pos1, sv_y+sv_y/5], xyB=[pos2, sv_y+sv_y/5], 
                                       coordsA='data', coordsB='data',
                                       axesA=ax1, axesB=ax2, arrowstyle='-', capstyle='butt', alpha=alpha,
                                       color=tra_color, lw=0.5+np.log2(support), connectionstyle=f"arc3,rad={0.1+np.log10(sv_rad)}")
            fig = ax2.get_figure()
            fig.add_artist(sv_patch)
        else:
            if pos1 > pos2:
                pos1, pos2 = pos2, pos1
                ori1, ori2 = ori2, ori1
            ax = sv_axes[chrom1]
            # ymax = ax1.get_ylim()[-1]
            ax.plot([pos1, pos1], [sv_y-sv_y/2, sv_y+sv_y/2], color=svcolor, linestyle='--', linewidth=0.5)
            ax.plot([pos2, pos2], [sv_y-sv_y/2, sv_y+sv_y/2], color=svcolor, linestyle='--', linewidth=0.5)
            new_xs, new_ys = make_spline_coordinates(pos1, pos2, y_offset=chrompair_y[(chrom1, chrom2)] / 150, margin_div=10*np.log10(pos2-pos1))
            ax.plot(new_xs, new_ys, color=svcolor, linewidth=0.5+np.log2(support), solid_capstyle='butt', alpha=alpha)
    add_sv_legend_to_axes(sv_axes, svcolors, supports)

def plot_segments_on_axes(plot_data, sv_axes):
    segment_y = 1
    offset_y = 2
    segment_linewidth = 5
    segment_color = 'tab:grey'
    _breakpoint_margin = 5000
    chrom_xlims = {}
    for brks in plot_data:
        if len(brks.segs) == 0:
            brks.get_segments()
        for seg in brks.segs: # plot segments
            brk1, brk2 = seg.brk1, seg.brk2
            chrom1, pos1 = brk1.chrom, brk1.pos
            chrom2, pos2 = brk2.chrom, brk2.pos
            assert chrom1 == chrom2, (chrom1, chrom2)
            ax_seg = sv_axes[chrom1]
            ax_seg.plot([pos1, pos2], [segment_y, segment_y], linewidth=segment_linewidth, color=segment_color, solid_capstyle='butt')
            ax_seg.tick_params(axis='x', labelrotation=90)
            ax_seg.set_ylim((segment_y-offset_y, segment_y+offset_y))
            chrom_xlims[chrom1] = ax_seg.get_xlim()
        for brk in [brks[0], brks[-1]]: # plot dangles for unpaired breakpoints
            chrom, pos, ori = brk.chrom, brk.pos, brk.ori
            pos1 = pos
            if chrom in chrom_xlims:
                (xmin, xmax) = chrom_xlims[chrom]
                breakpoint_margin = (xmax-xmin) / 20
            else:
                breakpoint_margin = _breakpoint_margin
            if ori == '+':
                pos2 = pos1 - breakpoint_margin
                pos3 = pos1 + breakpoint_margin
            else:
                pos2 = pos1 + breakpoint_margin
                pos3 = pos1 - breakpoint_margin
            ax_seg = sv_axes[chrom]
            ax_seg.plot([pos1, pos2], [segment_y, segment_y], linewidth=segment_linewidth, color=segment_color, solid_capstyle='butt')
            ax_seg.plot([pos1, pos3], [segment_y, segment_y], linewidth=segment_linewidth, alpha=0)
            ax_seg.tick_params(axis='x', labelrotation=90)
            ax_seg.set_ylim((segment_y-offset_y, segment_y+offset_y))
            with matplotlib.rc_context({'path.sketch': (1, 5, 1)}):
                ax_seg.plot([pos2, pos2], [segment_y - segment_y/4, segment_y + segment_y/4], linewidth=segment_linewidth/3, color='black')
