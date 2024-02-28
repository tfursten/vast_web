
import pandas as pd
import numpy as np
from scipy.spatial.distance import hamming

import plotly.express as px
import plotly.graph_objects as go
from plotly.colors import hex_to_rgb, label_rgb


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
    dist = get_distance_matrix(allele_table)
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



def add_imputed_vals_to_table(alleles_df, imputed_df):
    """
    Add imputed values to alleles table
    """
    for allele in alleles_df:
        col = alleles_df[allele]
        empty = col[col == ''].index
        for genome in empty:
            print(genome, allele)
            alleles_df.loc[genome, allele] = str(
                imputed_df.loc[(genome, allele)].Allele)
    return alleles_df


def alleles_to_numerical(alleles_df):
    """
    Convert alleles to numerical values
    """
    for column in alleles_df.columns:
        alleles_df[column] = alleles_df[column].astype('category').cat.codes
    return alleles_df


def format_data_for_optimization(alleles_df, imputed_df):
    """
    Format allele table by adding back in imputed values, 
    and converting alleles to numerical values.
    """
    return alleles_to_numerical(
        add_imputed_vals_to_table(alleles_df, imputed_df))


def get_starting_resolution(genomes):
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
    return np.array(scores)


def calculate_gini(opt_patterns, metadata):

    meta_index = np.where(metadata != "")[0]
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
    return np.array(scores)
    

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
        print(cur_pattern)
        # Add current pattern to all available patterns
        opt_pattern = np.add(patterns_arr, cur_pattern)
        print(opt_pattern)
        
        # Calculate entropy score for each pattern combined with current pattern
        # Use gini to calculate scores based on metadata categories otherwise use diversity
        if metadata is None:
            scores = calculate_scores(opt_pattern)
        else:
            scores = calculate_gini(opt_pattern, metadata)
        # Find minimum score, remove any targets that have already been selected in prev. iters.
        min_scores = np.setdiff1d(np.where(scores == scores.min())[0], targets)
        # Find index of pattern with min score
        # If randomize, return a random min score else pick first
        min_score_target = np.random.choice(min_scores) if randomize else min_scores[0]
        targets.append(min_score_target) # append index of best score
        target_names.append(patterns.columns.values[min_score_target])
        # update result pattern for next iteration
        result_pattern = np.unique(opt_pattern[min_score_target], return_inverse=True)[-1]
        
    return target_names





def draw_par_cats(resolution, metadata, alleles,
                   metadata_cat=None, font_size=12, palette='inferno'):
    alleles = alleles.astype(str)
    resolution = resolution.astype(str)
    resolution = resolution.T
    alleles.index = alleles.index.astype(str)
    resolution.index = resolution.index.astype(str)
    metadata.index = metadata.index.astype(str)
    # print(resolution)
    # print(alleles)
    # print(metadata)
    # print(metadata_cat)
    for genome in resolution.index:
        for locus in resolution.columns:
            resolution.loc[genome, locus] = f"{alleles.loc[genome, locus]}-{resolution.loc[genome, locus]}"
    
    print(resolution)
    if metadata_cat:
        metadata[f'{metadata_cat}_numerical'] = metadata[metadata_cat].astype('category').cat.codes
        resolution = metadata[[metadata_cat, f'{metadata_cat}_numerical' ]].join(resolution)    
        print("HERE", resolution)    
        resolution = resolution.sort_values(list(resolution.columns.values))
        color = resolution[f'{metadata_cat}_numerical']
        resolution = resolution.drop(f'{metadata_cat}_numerical', axis=1)
    else:
        resolution['Genomes'] = np.arange(0, resolution.shape[0])
        color = resolution['Genomes']
        resolution = resolution.drop('Genomes', axis=1)

    dims = []

    for c in resolution:
        dims.append(
            go.parcats.Dimension(
                values = resolution[c],
                label = c
            ))

    dims.append(
        go.parcats.Dimension(
            values = resolution.index,
            label = "Genomes",
        )
    )

    fig = go.Figure(data = [go.Parcats(dimensions=dims,
            line={'color': color, 'colorscale': palette},
            bundlecolors=True,         
            hoveron='category', hoverinfo='count',
            labelfont={'size': font_size, 'family': 'Arial', 'color': 'black'},
            tickfont={'size': font_size, 'family': 'Arial', 'color': 'black'},
            arrangement='freeform')]) 

    fig.update_layout(
        height=max(800, (resolution.shape[0] * 10))
    )
#     fig.update_layout(
#         margin=dict(l=200, r=200, t=50, b=20),
#     )
    html = fig.to_html()
    return html


