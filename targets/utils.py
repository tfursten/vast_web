import os
import pandas as pd
import numpy as np
from scipy.spatial.distance import hamming

import plotly.express as px
import plotly.graph_objects as go
from plotly.colors import hex_to_rgb, label_rgb
from grapetree import grapetree
from io import StringIO

from Bio import Phylo
import igraph

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
# from ete3 import (
#     Tree, TreeStyle,
#     NodeStyle,  TextFace,
#     RectFace, CircleFace,
#     PieChartFace)

# os.environ['QT_QPA_PLATFORM']='offscreen'


def replace_missing_with_unique(arr1, arr2, val):
    """
    Recategorizes missing data as all unique values
    input:
        arr1: Encoded array for alleles for genome 1
        arr2: Encoded array for alleles for genome 2
        val: Value representing the missing category
    output:
        arr1 and arr2 with missing category number replaced 
        with unique values.
    """
    max_value = np.concatenate([arr1, arr2]).max() + 1
    new_arr_1 = []
    new_arr_2 = []
    for i in arr1:
        if i == val:
            i = max_value
            max_value += 1
        new_arr_1.append(i)
    for i in arr2:
        if i == val:
            i = max_value
            max_value += 1
        new_arr_2.append(i)
    return np.array(new_arr_1), np.array(new_arr_2)
            
    

def get_distance_matrix(allele_table):
    """
    Calculate the distance between each genome.
    Hamming distance is proportion of matches.
    All missing data is treated like a mismatch.

    Input:
        allele_table: DataFrame with genome ids as index and Loci as columns
    Output:
        2D distance matrix with genome ids as columns and index 
    """
    n_genomes = allele_table.shape[0]
    dist = np.zeros((n_genomes, n_genomes))
    for i in range(n_genomes - 1):
        for j in range(i + 1, n_genomes):
            genome1_alleles = allele_table.iloc[i]
            genome2_alleles = allele_table.iloc[j]

            categories = np.unique(
                np.concatenate([genome1_alleles.values, genome2_alleles.values]))
            genome1_alleles_enc = np.searchsorted(categories, genome1_alleles)
            genome2_alleles_enc = np.searchsorted(categories, genome2_alleles)
            # to avoid genomes with missing data from being more similar to eachother
            # add unique values for each missing so they are treated as differences

            blank_idx = np.where(categories=='')[0]
            if len(blank_idx):
                genome1_alleles_enc, genome2_alleles_enc = replace_missing_with_unique(
                    genome1_alleles_enc, genome2_alleles_enc, blank_idx[0])
            d = hamming(
                        genome1_alleles_enc,
                        genome2_alleles_enc)
            dist[i][j] = d
            dist[j][i] = d

    return pd.DataFrame(dist, index=allele_table.index, columns=allele_table.index)


def get_grapetree_distance_matrix(alleles):
    # fill in missing values using graptree format convert to string
    file_buffer = StringIO()
    alleles.replace("", "-").reset_index().to_csv(file_buffer, sep="\t", index=False)
    file_str = file_buffer.getvalue()

    dist = pd.read_csv(StringIO(grapetree.backend(
        profile=file_str,
        method="distance",
        matrix_type="symmetric")),
        sep="\\s+", skiprows=1, header=None, index_col=0)
    dist.columns = dist.index
    return dist
    

def get_grapetree_newick(alleles):
    """
    Get newick tree from grapetree
    """
    file_buffer = StringIO()
    alleles.replace("", "-").reset_index().to_csv(file_buffer, sep="\t", index=False)
    file_str = file_buffer.getvalue()
    return StringIO(grapetree.backend(
        profile=file_str,
        method="MSTreeV2")).getvalue()
 
        
def get_weighted_likelihood_allele(alleles, dists):
    """
    Get most likely allele weighted by distance
    """
    dat = pd.DataFrame.from_dict({'Alleles': alleles, 'Dists': dists})
    dat['Inverse'] = (1 - dat['Dists'])**3
    dat = dat[dat['Alleles'] != '']
    counts = dat.groupby('Alleles')['Inverse'].sum().idxmax()
    return str(counts)

def impute_missing_alleles(allele_table):
    """
    Impute missing alleles, returns a dataframe with rows 
    indicating [genome, locus, imputed allele]
    """
    print(allele_table.head())
    print("Getting Distance Matrix")
    dist = get_grapetree_distance_matrix(allele_table)
    print("Done Getting Matrix")
    imput = []
    for genome, row in allele_table.iterrows():
        # list loci with missing data for this genome
        missing = row[row == ""]
        if len(missing):
            for locus in missing.index:
                # impute missing allele weighted by nearest neighbors
                # update allele table with imputed value
                imput.append([genome, locus, get_weighted_likelihood_allele(
                    allele_table[locus], dist.loc[genome])])
                
    return pd.DataFrame(imput, columns=['Genome', 'Locus', 'Allele'])



# def add_imputed_vals_to_table(alleles_df, imputed_df):
#     """
#     Add imputed values to alleles table
#     """
#     for allele in alleles_df:
#         col = alleles_df[allele]
#         empty = col[col == ''].index
#         for genome in empty:
#             print(genome, allele)
#             alleles_df.loc[genome, allele] = str(
#                 imputed_df.loc[(genome, allele)].Allele)
#     return alleles_df


def alleles_to_numerical(alleles_df):
    """
    Convert alleles to numerical values
    """
    for column in alleles_df.columns:
        alleles_df[column] = alleles_df[column].astype('category').cat.codes
    return alleles_df


def format_data_for_optimization(alleles_df):
    """
    Convert alleles to numerical category values.
    """
    return alleles_to_numerical(alleles_df)
        


def get_starting_resolution(genomes):
    """
    Starting resolution is all genomes in same group == all zeros
    """
    return pd.DataFrame(np.zeros(len(genomes)), index=genomes, columns=['Start'])


def get_resolution(alleles_df):
    """
    Given a matrix of alleles for different genomes, return an 
    array that clusters genomes with the same genotype
    for the provided SNPs.
    Given: 
            [LOCI]
    [GENOMES]  0  0  0 
               1  0  1  
               1  1  2
               1  0  1
    Return:
    [0, 1, 2, 1]
    """
    return np.unique(
        alleles_df.values.astype('U'), axis=0, return_inverse=True)[-1]
    

def calculate_scores(opt_patterns):
    scores = []
    for row in opt_patterns:
        _, counts = np.unique(row, return_counts=True)
        scores.append(sum([c**2 - c for c in counts]))
    scores = np.array(scores)
    if np.max(scores) == 0:
        return np.zeros(len(scores))
    else:
        return (scores-np.min(scores))/(np.max(scores)-np.min(scores))

def calculate_gini(opt_patterns, metadata):
    meta_index = np.where(metadata != "nan")[0]
    meta_values = np.unique(metadata.iloc[meta_index], return_inverse=True)[1]
    scores = []
    for i, row in enumerate(opt_patterns):
        sizes = []
        gini_impurity = []
        row = row[meta_index]
        for group in np.unique(row):
            group = np.where(row==group)[0]
            group_sz = len(group)
            sizes.append(group_sz)
            gini_impurity.append(
                1 - sum([(n/group_sz)**2 for n in
                         np.unique(meta_values[group],
                                   return_counts=True)[1]]))        
        sizes = [s/len(row) for s in sizes]
        scores.append(sum([g*s for g, s in zip(gini_impurity, sizes)]))
    if np.max(scores) == 0:
        return np.zeros(len(scores))
    else:
        return (scores-np.min(scores))/(np.max(scores)-np.min(scores))

    

def optimization_loop(
    patterns, starting_pattern, n_targets, randomize, metadata=None):
    """
    Run optimization loop and return a list of chosen targets
    """
    # Add constant value to starting pattern to easily add to pattern matrix
    # Constant needs to be order of magnitude higher than any value in pattern matrix
    const = 100 ** len(str(patterns.shape[0]))
    result_pattern = starting_pattern
    targets = []
    target_names = []
    iterations = min(n_targets, patterns.shape[1]) # reduce iterations to remaining targets if less than n_targets
    patterns_arr = patterns.values.T
    for i in range(iterations):
        # Add const to current pattern
        cur_pattern = np.multiply(result_pattern, const)
        # Add current pattern to all available patterns
        opt_pattern = np.add(patterns_arr, cur_pattern)
        
        # Calculate entropy score for each pattern combined with current pattern
        # Use gini to calculate scores based on metadata categories otherwise use diversity
        if metadata is None:
            scores = calculate_scores(opt_pattern)
        else:
            # pick scores based mostly by gini impurity but slightly modified by resoltuion
            scores = calculate_gini(opt_pattern, metadata) + calculate_scores(opt_pattern) ** 3

        # Find minimum score, remove any targets that have already been selected in prev. iters.
        min_scores = np.setdiff1d(np.where(scores == scores.min())[0], targets) # index of best scores (lowest)
        # Find index of pattern with min score
        # If randomize, return a random min score else pick first
        min_score_target = np.random.choice(min_scores) if randomize else min_scores[0]
        targets.append(min_score_target) # append index of best score
        target_names.append(patterns.columns.values[min_score_target])
        # update result pattern for next iteration
        result_pattern = np.unique(opt_pattern[min_score_target], return_inverse=True)[-1]
        
    return target_names



def draw_par_cats(
        resolution,
        metadata_cat=None,
        font_size=12, palette='inferno'):
    """
    Draw parallel categories plot from resolution and metadata.
    Returns HTML of figure from plotly plot
    """
    
    resolution = resolution.reset_index()
    resolution['count'] = resolution.groupby('Resolution')['Resolution'].transform('count')
    if metadata_cat:
        resolution[f'{metadata_cat}_numerical'] = resolution[metadata_cat].astype('category').cat.codes
        resolution = resolution.sort_values([metadata_cat, 'count', 'Resolution'], ascending=False)
        color = resolution[f'{metadata_cat}_numerical'] # don't do any resorting after color is set
        resolution = resolution[[metadata_cat, 'Resolution', 'Genome']]
    else:
        resolution = resolution.sort_values(['count', 'Resolution'], ascending=False)
        resolution['genome_numerical'] = resolution['Genome'].astype('category').cat.codes
        color = resolution['genome_numerical'] # Don't do any sorting after color is set
        resolution = resolution[['Resolution', 'Genome']]


    dims = []

    for c in resolution:
        dims.append(
            go.parcats.Dimension(
                values = resolution[c],
                label = 'Resolution with selected loci' if c=="Resolution" else c
            ))

    fig = go.Figure(data = [go.Parcats(dimensions=dims,
            line={'color': color, 'colorscale': palette},
            bundlecolors=False,
            hoveron='category', hoverinfo='count',
            labelfont={'size': font_size, 'family': 'Arial', 'color': 'black'},
            tickfont={'size': font_size, 'family': 'Arial', 'color': 'black'},
            arrangement='freeform')])

    fig.update_layout(
        height=max(800, (resolution.shape[0] * 10))
    )
    
    html = fig.to_html()
    return html

def drop_nans(df, thresh_percent, axis=0):
    """
    Drop rows or columns of a dataframe that have Nans over a threashold
    """
    # Calculate the minimum number of non-NaN values required based on the threshold percentage
    threshold_count = int(df.shape[(not bool(axis))] * thresh_percent)
    # Drop rows with more than the specified threshold count of NaN values
    res = df.dropna(thresh=threshold_count, axis=axis)
    return res


# def draw_tree(tree, metadata, style, font_size, palette, legend=True):

#     # Set up colormaps
#     metadata.index = metadata.index.map(str)

#     meta_cat = metadata.columns.values[0]
#     meta_data_code_map = {
#         v: k for k, v in enumerate(
#             metadata[meta_cat].value_counts().index)}
#     metadata['code'] = metadata[meta_cat].apply(lambda x: meta_data_code_map.get(x))
#     # Choose a built-in colormap 
#     cmap = plt.get_cmap(palette)
#     print(metadata)
#     # Normalize values to be between 0 and 1
#     norm = plt.Normalize(0, metadata['code'].nunique())

#     # Create a ScalarMappable object to map values to colors
#     sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
#     meta_colors_rgba = {
#         k: v for k, v in zip(
#             metadata['code'].unique(),
#             sm.to_rgba(metadata['code'].unique()))}

#     meta_colors_hex = {k: mcolors.to_hex(v)
#                        for k, v in meta_colors_rgba.items()}

#     # Get the color-mapped values
#     metadata['color'] = metadata['code'].map(meta_colors_hex)
    
#     tree = Tree(tree, format=1)
#     tree.ladderize()
#     for node in tree.traverse():
#         nstyle = NodeStyle()
#         if node.is_leaf():
#             # Get metadata
#             genome_id = node.name

#     #         Add sample name face
#             name_face = TextFace(
#                 genome_id, fgcolor="black",
#                 fsize=font_size, penwidth=6, ftype='Arial')
#             name_face.margin_left = 5
#             node.add_face(name_face, column=2)

#             # Rectangle sizes
#             rect_height = 10
#             rect_width = 20
#             outline = '#f0f0f0'
            
#             # Metadata Face
#             meta_col = metadata.loc[genome_id]['color']
            
#             meta_rect = RectFace(rect_width, rect_height, outline, meta_col)
#             node.add_face(meta_rect, column=1)

#         # Add pie chart to each node
#         if node.is_leaf():
#             continue
#     #         sample_leaf_nodes = [node.name]
#         else:
#     #     if nodeid in pick_nodes:
#             leaf_names = node.get_leaf_names()
#             genome_leaf_nodes = list(leaf_names)
#             meta_pie_data = np.array([
#                 metadata.loc[gen]['code'] for gen in genome_leaf_nodes])
#             meta_pie_data.sort()
#             if len(meta_pie_data):
#                 percent = [(1/len(meta_pie_data)) * 100 for _ in meta_pie_data]
#                 colors = [meta_colors_hex.get(val) for val in meta_pie_data]
#                 line_color = None if len(percent) > 15 else "white"
#                 pie_size = 20
#                 pie = PieChartFace(
#                             percent, pie_size, pie_size,
#                             colors=colors, line_color=line_color)
#                 node.add_face(pie, column=1)
#         # Last set node style
#         nstyle['hz_line_width'] = 1
#         nstyle['vt_line_width'] = 1
#         nstyle['fgcolor'] = "black"
#         nstyle['size'] = 1

#         node.set_style(nstyle)
        
#     ts = TreeStyle()
#     ts.mode = 'c' if style == "Circle" else 'r'
#     ts.margin_bottom = 15
#     ts.margin_left = 30
#     ts.margin_top = 10
#     ts.margin_right = 10
#     ts.min_leaf_separation = 1
#     ts.allow_face_overlap = False
#     ts.show_leaf_name = False
#     # ts.optimal_scale_level = "full"

    
#     if legend:
#         ts.legend_position=3

#         # Draw Legend
#         dot_size = 6
#         t_margin_right = 10
#         t_margin_left = 8
#         t_margin_top = 2

#         #     meta_data_code_map (name to code) meta_colors_hex (code to color)
#         for name, code in meta_data_code_map.items():
#             color = meta_colors_hex.get(code)
#             ts.legend.add_face(CircleFace(dot_size, color), column=2)
#             textface = TextFace(name, fsize=font_size, tight_text=True, ftype='Arial')
#             textface.margin_top = t_margin_top
#             textface.margin_right = t_margin_right
#             textface.margin_left = t_margin_left
#             ts.legend.add_face(textface, column=1)
        
#     return tree, ts




def draw_tree(tree, metadata, style, font_size, palette, legend=True):
    tree = Phylo.read(StringIO(tree), "newick")
    G = convert_to_igraph(tree)

    # Get edges and vertices
    verticies = [v.index for v in G.vs]
    edges = [e.tuple for e in G.es]
    node_colors = [
        'Black' for v in G.vs]

    labels = [v['name'] for v in G.vs]
    display_labels = ['' if 'Inner' in label else label for label in labels]
    node_size = [1 if d_label == '' else 12 for d_label in display_labels]

    layt = G.layout('kk')
    N = len(layt)# N is equal to len(V)=len(G.vs)
    



def convert_to_igraph(tree):
    #Convert a Biopython Tree object to an igraph Graph.
    def add_edge(graph, node1, node2):
        graph.add_edge(node1.name, node2.name)
           

    def build_subgraph(graph, top):
        """Traverse  the Tree, and retrieve  graph edges and nodes."""
        for clade in top:
            graph.add_vertex(name=clade.root.name)
            add_edge(graph, top.root, clade.root)
            build_subgraph(graph, clade)

    if tree.rooted:
        G = igraph.Graph(directed=True)
    else:
        G = igraph.Graph()
    G.add_vertex(name=str(tree.root))
    build_subgraph(G, tree.root)
    return G 




# def draw_par_cats(resolution, metadata, alleles,
#                    metadata_cat=None, font_size=12, palette='inferno'):
#     alleles = alleles.astype(str)
#     resolution = resolution.astype(str)
#     # resolution = resolution.T
#     alleles.index = alleles.index.astype(str)
#     resolution.index = resolution.index.astype(str)
#     metadata.index = metadata.index.astype(str)
#     print("IN PARCAT")
#     print(metadata, alleles, resolution)

#     for genome in resolution.index:
#         for locus in resolution.columns:
#             resolution.loc[genome, locus] = f"{alleles.loc[genome, locus]}-{resolution.loc[genome, locus]}"
    
#     if metadata_cat:
#         metadata[f'{metadata_cat}_numerical'] = metadata[metadata_cat].astype('category').cat.codes
#         resolution = metadata[[metadata_cat, f'{metadata_cat}_numerical' ]].join(resolution)    
#         resolution = resolution.sort_values(list(resolution.columns.values))
#         color = resolution[f'{metadata_cat}_numerical']
#         resolution = resolution.drop(f'{metadata_cat}_numerical', axis=1)
#     else:
#         resolution['Genomes'] = np.arange(0, resolution.shape[0])
#         color = resolution['Genomes']
#         resolution = resolution.drop('Genomes', axis=1)

#     dims = []

#     for c in resolution:
#         dims.append(
#             go.parcats.Dimension(
#                 values = resolution[c],
#                 label = c
#             ))

#     dims.append(
#         go.parcats.Dimension(
#             values = resolution.index,
#             label = "Genomes",
#         )
#     )

#     fig = go.Figure(data = [go.Parcats(dimensions=dims,
#             line={'color': color, 'colorscale': palette},
#             bundlecolors=True,         
#             hoveron='category', hoverinfo='count',
#             labelfont={'size': font_size, 'family': 'Arial', 'color': 'black'},
#             tickfont={'size': font_size, 'family': 'Arial', 'color': 'black'},
#             arrangement='freeform')]) 

#     fig.update_layout(
#         height=max(800, (resolution.shape[0] * 10))
#     )
# #     fig.update_layout(
# #         margin=dict(l=200, r=200, t=50, b=20),
# #     )
#     html = fig.to_html()
#     return html



