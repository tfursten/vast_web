
import pandas as pd
import numpy as np
from scipy.spatial.distance import hamming




from scipy.spatial.distance import hamming
import numpy as np

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
    return counts

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
