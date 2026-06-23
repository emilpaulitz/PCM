from pycirclize import Circos  # THIS CHANGES MARGINS!!!
import numpy as np
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage, leaves_list, to_tree
from scipy.spatial.distance import pdist
import matplotlib as mpl
from matplotlib.patches import Patch
from random import shuffle
import pandas as pd

# # TODO make giving a tree accessible from wrapper function Function definition for circular plot

def categorize(x, categories):
    if x == 0:
        return float('Nan')
    for i, c in enumerate(categories):
        if x >= c:
            return i

def format_data(df, categories):
    # categorize each cell into how often the corresponding reaction occurs
    data = (df.mean() * df).map(lambda x: categorize(x, categories))

    # count occurrences of each category and summarize in a df
    return pd.concat({i: (data == i).sum(axis=1) for i in range(len(categories))}, axis=1)

# from https://stackoverflow.com/questions/28222179/save-dendrogram-to-newick-format/31878514#31878514
def get_newick(node, parent_dist, leaf_names, newick='') -> str:
    """
    Convert sciply.cluster.hierarchy.to_tree()-output to Newick format.

    :param node: output of sciply.cluster.hierarchy.to_tree()
    :param parent_dist: output of sciply.cluster.hierarchy.to_tree().dist
    :param leaf_names: list of leaf names
    :param newick: leave empty, this variable is used in recursion.
    :returns: tree in Newick format
    """
    if node.is_leaf():
        return "%s:%.2f%s" % (leaf_names[node.id], parent_dist - node.dist, newick)
    else:
        if len(newick) > 0:
            newick = "):%.2f%s" % (parent_dist - node.dist, newick)
        else:
            newick = ");"
        newick = get_newick(node.get_left(), node.dist, leaf_names, newick=newick)
        newick = get_newick(node.get_right(), node.dist, leaf_names, newick=",%s" % (newick))
        newick = "(%s" % (newick)
        return newick

def original_value(x):
    return np.exp2(x) - 0.5
    #return x  # for difference data

def scaled_value(x):
    return np.log2(x+0.5)

def scale_order_data(data, scale=True, dists=None, nwk_file=None, given_order = None, index_transform=None):

    if given_order is not None and not all([acc in data.index for acc in given_order]):
        print('Error: Given order of accessions is invalid: not all accessions are in the data index')
        return None, None

    # transform the data to make large values appear dark and reduce the dominance of accessory genes
    if scale:
        data = data.map(scaled_value)

    # order data according to a clustering, s.th. the distance between neighbors is minimal
    if nwk_file is None:
        if dists is None:
            print(f'Error: {__name__} needs either dists or nwk_file')
        Z = linkage(dists, optimal_ordering = True)
        tree = to_tree(Z)
        nwk = get_newick(tree, tree.dist, data.index)
        if given_order is None:
            return nwk, data.iloc[leaves_list(Z), ]
        return nwk, data.loc[given_order, ]
    
    else:
        # read tree file
        nwk = ''
        with open(nwk_file, 'r') as f:
            for line in f:
                nwk += line

        # we only want the order of species
        nnwk = ''
        in_score = False
        for c in nwk:
            if c in '();':
                continue
            if c == ':':
                in_score = True
                continue
            if c == ',':
                in_score = False
            if in_score and (c.isnumeric() or c == '.'):
                continue
            nnwk += c

        # sort dataframe according to this order
        if index_transform is None:
            return data.loc[nnwk.split(','),]
        else:
            return data.loc[[index_transform[n] for n in nnwk.split(',')], ]

# Draw circular figure
def draw_circular_plot(data, filename, categories, round_legend_labels_ndigits = None,
                       obj ='Gene', label_formatter=lambda x: x):
    colormap = "inferno_r"

    # Define sectors
    num_genomes = data.shape[0]
    sectors = {"Heatmaps": num_genomes, 'Labels': 0}

    # Initialize Circos plot
    inner_space = 0.35
    circos = Circos(sectors, space=0, start=0.25 * 360, end = 1 * 360)

    # Add heatmap tracks
    vmin = data.min().min()
    vmax = data.max().max()

    # these determine the inner/outer border of tracks
    tmp = [1.] + categories
    ranges = [(tmp[-(i + 1)], tmp[-(i + 2)]) for i in range(len(categories))]

    # normalize ranges
    ranges = [((a - 1.) * (1 - inner_space) + 1, (b - 1.) * (1 - inner_space) + 1) for a,b in ranges]
    ncategories = len(ranges)
    
    for sector in circos.sectors:

        if sector.name == 'Heatmaps':
            heatmap_tracks = {}
            for i, (l, u) in enumerate(ranges):

                # add track
                heatmap_track = sector.add_track((l * 100, u * 100))
                heatmap_track.axis(fc="none")
                # reverse order of data just for looks
                heatmap_track.heatmap(data[i].to_numpy()[::-1], cmap=colormap, vmin=vmin, vmax=vmax)

                # draw ticks and labels
                if i == ncategories - 1:
                    heatmap_track.xticks(
                        [i + 0.5 for i in range(num_genomes)], [label_formatter(label) for label in data.index[::-1]], # remove 'Cyc'
                        label_orientation="vertical", show_bottom_line=True, label_size=10, line_kws={'ec': 'grey'})
                    
                # for labels
                heatmap_tracks[i / ncategories] = heatmap_track
            
        if sector.name == 'Labels':
            padding = ' ' * 2
            percs = [100, 80, 60, 40, 20, 0]
            for i, track in heatmap_tracks.items():
                for perc in percs:
                    if perc == round(i*100):
                        # the i ranges from 1 to 0 but the other way around than we need it
                        circos.text(f'{padding}{100 - perc}%', r = track.r_lim[0], deg = 0, horizontalalignment='left')

            circos.text(f'{padding}{obj} in 0%\n{padding}of accessions', r = track.r_lim[1], deg = 0, horizontalalignment='left')


    # colorbar does not work becuase we can not get a logarithmic scale :(
    #circos.colorbar(vmin=original_value(vmin), vmax=original_value(vmax), cmap=colormap, bounds=(0.75, 0.55, 0.1, 0.4))

    # sample colormap
    nlabels = 8
    cmap_sampled = mpl.colormaps[colormap].resampled(nlabels)

    # define my own range function because otherwise the rounding error will become too large
    def float_range(vmin, vmax, step):
        i = vmin
        while i <= vmax:
            yield i
            i += step

    step = (vmax - vmin) / (nlabels - 1)
    handles = [Patch(color=cmap_sampled((i-vmin)/(vmax-vmin)), linewidth=10,
                    #label=int(original_value(i))) for i in float_range(vmin, vmax, step)]
                    label=round(original_value(i), round_legend_labels_ndigits)) for i in float_range(vmin, vmax + 1, step)]

    fig = circos.plotfig()  # needed for some reason
    _ = circos.ax.legend(handles=handles, bbox_to_anchor=(0.75, 0.55), loc="lower left", fontsize=12)

    # Save and display the plot
    fig.savefig(filename, dpi = 400)

# Draw circular figure
def draw_circular_plot_with_tree(data, filename, categories, nwk, round_legend_labels_ndigits = None,
                                 obj = 'Gene', accs = None, fontsize=12, fontsize_legend=12,
                                 legend_bbox = (0.75, 0.55), label_formatter = lambda x: x, 
                                 tree_kws = {}, heatmap_kws = {}, axis_kws = None, species_accs = 'accessions'):

    colormap = "inferno_r"

    # Define sectors
    num_genomes = data.shape[0]

    # Initialize Circos plot
    inner_space = 0.35
    circos, tv = Circos.initialize_from_tree(nwk, r_lim = (0, inner_space * 100),
                                                start=0.25 * 360, end=1*360,
                                                leaf_label_rmargin=(1-inner_space + 0.05) * 100,
                                                ignore_branch_length=True, 
                                                label_formatter=lambda s: '', 
                                                line_kws = tree_kws)
    sector = tv.track.parent_sector

    # Add heatmap tracks
    vmin = data.min().min()
    vmax = data.max().max()

    # these determine the inner/outer border of tracks
    tmp = [1.] + categories 
    ranges = [(tmp[-(i + 1)], tmp[-(i + 2)]) for i in range(len(categories))]

    # normalize ranges
    ranges = [((a - 1.) * (1 - inner_space) + 1, (b - 1.) * (1 - inner_space) + 1) for a,b in ranges]
    ncategories = len(ranges)

    heatmap_tracks = {}
    for i, (l, u) in enumerate(ranges):

        # add track
        heatmap_track = sector.add_track((l * 100, u * 100))
        if axis_kws is None:
            axis_kws = {}
        if 'fc' not in axis_kws:
            axis_kws['fc'] = 'none'
        heatmap_track.axis(**axis_kws)
        # reverse order of data just for looks
        heatmap_track.heatmap(data[i].to_numpy()[::-1], cmap=colormap,
                              vmin=vmin, vmax=vmax, rect_kws = heatmap_kws)

        # draw ticks and labels
        if i == ncategories - 1 :
            xticks = list()
            xticks_names = list()
            for j, label in zip(range(num_genomes), data.index[::-1]):
                if accs is None or label in accs:
                    xticks.append(j + 0.5)
                    xticks_names.append(label_formatter(label))

            heatmap_track.xticks(
                xticks, xticks_names,
                label_orientation="vertical", show_bottom_line=True, label_size=fontsize, line_kws={'ec': 'grey'})#, text_kws={'fontsize':fontsize})
            
        # for labels
        heatmap_tracks[i / ncategories] = heatmap_track

    # add percentages labels
    padding = ' ' * 2
    percs = [100, 80, 60, 40, 20, 0]
    for i, track in heatmap_tracks.items():
        for perc in percs:
            if perc == round(i*100):
                # the i ranges from 1 to 0 but the other way around than we need it
                circos.text(f'{padding}{100 - perc}%', r = track.r_lim[0], deg = 0, horizontalalignment='left', fontsize=fontsize_legend)
    circos.text(f'{padding}{obj} in 0%\n{padding}of {species_accs}', r = track.r_lim[1], deg = 0, horizontalalignment='left', fontsize=fontsize_legend)

    # sample colormap
    nlabels = 8
    cmap_sampled = mpl.colormaps[colormap].resampled(nlabels)

    # define my own range function because otherwise the rounding error will become too large
    def float_range(vmin, vmax, step):
        i = vmin
        while i < vmax:
            yield i
            i += step

    step = (vmax - vmin) / (nlabels - 1)
    handles = [Patch(color=cmap_sampled((i-vmin)/(vmax-vmin)), linewidth=10,
                    #label=int(original_value(i))) for i in float_range(vmin, vmax, step)]
                    label=round(original_value(i), round_legend_labels_ndigits)) for i in float_range(vmin, vmax + 0.001, step)]

    fig = circos.plotfig()  # needed for some reason
    _ = circos.ax.legend(handles=handles, bbox_to_anchor=legend_bbox, loc="lower left", fontsize=fontsize_legend)

    # Save and display the plot
    if filename:
        fig.savefig(filename, dpi = 400, transparent=True)

def plot_circular(pav, filename, ncategories, obj='Gene', accs = None, given_order = None,
                  fontsize = 12,  fontsize_legend=12, label_formatter=lambda x: x, 
                  tree_kws = {}, heatmap_kws = {}, axis_kws = None, species_accs = 'accessions'):

    # giving uneven categories has the problem that more percentages will be summed and thus the color will be misleading
    categories = [x/ncategories for x in range(ncategories)][::-1]
    # genome nwk: 'SpeciesTree_rooted.txt' # for genomes
    data = format_data(pav, categories = categories)
    nwk, data = scale_order_data(data, dists = pdist(pav), given_order = given_order)

    #draw_circular_plot(data, filename, nwk = nwk, categories = categories, obj=obj)
    # accs determines which accs will be labeled
    draw_circular_plot_with_tree(data, filename, nwk = nwk, categories = categories, obj=obj,
                                 accs = accs, fontsize=fontsize, fontsize_legend=fontsize_legend, 
                                 label_formatter = label_formatter, tree_kws = tree_kws,
                                 heatmap_kws = heatmap_kws, axis_kws = axis_kws, species_accs = species_accs)
    plt.show()

    return data
