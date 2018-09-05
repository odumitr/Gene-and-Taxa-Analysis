import pandas as pd
import matplotlib.pyplot as plt
import math
import statsmodels.formula.api as smf

def load_data():
	
    inFile = open("../../data/Scores_Sizes.txt")

    scores = []
    gene_profile_sizes = []
    taxon_profile_sizes = []

    for line in inFile:
        if "Score" not in line:
            data = line.strip().split("\t")
            score = float(data[6])
            scores.append(score)
            gene_profile_sizes.append(int(data[1]))
            taxon_profile_sizes.append(int(data[4]))
    inFile.close()
    
    return scores, gene_profile_sizes, taxon_profile_sizes


def generate_avg_scores(scores, taxon_profile_sizes, gene_profile_sizes):

    if (len(scores) == 0 or len(taxon_profile_sizes) == 0 or len(gene_profile_sizes) == 0):
        print('missing values for scores/taxon profile sizes/gene profile sizes')
        return
    
    gene_taxon_scores = pd.DataFrame(
    {
        'gene': gene_profile_sizes,
        'taxon': taxon_profile_sizes,
        'scores': scores
    })    

    gene_taxon_scores = gene_taxon_scores.groupby(['taxon','gene'], as_index=False)['scores'].mean() 
    gene_taxon_scores.rename(columns={'scores':'avg_scores'}, inplace=True)    
    print(gene_taxon_scores.head())
    return gene_taxon_scores

def ols_multiple_regression(input_data):

    results = smf.ols(formula='avg_scores ~ taxon + gene', data=input_data).fit()
    print(results.params)
    print(results.summary())

    return results

def compute_residuals(regression_results, gene_taxon_avg_scores):

    

def main():
	
    scores, gene_profile_sizes, taxon_profile_sizes = load_data()
    gene_taxon_avg_scores = generate_avg_scores(scores, taxon_profile_sizes, gene_profile_sizes)
    regression_results = ols_multiple_regression(gene_taxon_avg_scores)    
    compute_residuals(regression_results, gene_taxon_avg_scores)

main()