def main():
	getScores()
	scores, rawTaxonProfileSizes, rawGeneProfileSizes = loadProfiles()

	# singleRegression(scores, rawGeneProfileSizes)
	# plotRawResiduals(scores, rawTaxonProfileSizes, rawGeneProfileSizes)

	# results = sampleBins(scores, rawTaxonProfileSizes, rawGeneProfileSizes)

	plotGrid(scores, rawTaxonProfileSizes, rawGeneProfileSizes)


def getScores():
	resultDir = "../results/"
	if not os.path.exists(resultDir):
		os.makedirs(resultDir)
	size = loadProfileSizes()


def loadProfileSizes():
	profileSize = dict()
	inFile = open("../data/ProfileSizes.txt")

	for line in inFile:
		entity, size = line.strip().split("\t")
		entity = entity.replace("#profile","")
		entity = entity.replace("http://purl.obolibrary.org/obo/","")
		profileSize[entity] = int(size)

	inFile.close()
	return profileSize


def loadProfiles():
	inFile = open("../data/Scores_Sizes.txt")

	scores = []
	geneProfileSizes = []
	taxonProfileSizes = []

	rawscores = []
	rawGeneProfileSizes = []
	rawTaxonProfileSizes = []

	for line in inFile:
		if "Score" not in line:
			data = line.strip().split("\t")
			score = float(data[6])
			scores.append(score)
			rawGeneProfileSizes.append(int(data[1]))
			rawTaxonProfileSizes.append(int(data[4]))
	inFile.close()
	return scores, rawTaxonProfileSizes, rawGeneProfileSizes


def singleRegression(scores,rawGeneProfileSizes):
	residuals = []
	sizes = [rawGeneProfileSizes]
	ones = np.ones(len(sizes[0]))
	X = sm.add_constant(np.column_stack((sizes[0], ones)))
	for ele in sizes[1:]:
		X = sm.add_constant(np.column_stack((ele, X)))
	results = sm.OLS(scores, X).fit()
	print
	print
	print "Single Regression Results"
	print results.summary()

	predicted = results.fittedvalues
	for i in range(0,len(scores)):
		actualscore = scores[i]
		predictedscore = predicted[i]
		residual = actualscore - predictedscore
		residuals.append(residual)

	plt.scatter(np.array(predicted), np.array(residuals))
	plt.xlabel('Predicted')
	plt.ylabel('Raw Residual')
	plt.savefig('../ACMManuscript/Figures/SingleRegression_PredictedvsResiduals.png')


def multiRegression(scores, taxonProfileSizes, geneProfileSizes):
	sizes = [taxonProfileSizes, geneProfileSizes]
	ones = np.ones(len(sizes[0]))
	X = sm.add_constant(np.column_stack((sizes[0], ones)))
	for ele in sizes[1:]:
		X = sm.add_constant(np.column_stack((ele, X)))
	results = sm.OLS(scores, X).fit()
	return results


def plotRawResiduals(scores, rawTaxonProfileSizes, rawGeneProfileSizes):
	geneAxis = []
	taxonAxis = []
	residuals = []
	count = 0
	rawResults = multiRegression(scores, rawTaxonProfileSizes, rawGeneProfileSizes)

	print
	print
	print "Raw Regression Results"
	print rawResults.summary()
	predicted = rawResults.fittedvalues
	for i in range(0, len(scores)):
		predictedScore = predicted[i]
		actualScore = scores[i]
		residual = actualScore - predictedScore
		residuals.append(residual)

	plt.scatter(np.array(predicted), np.array(residuals))
	plt.xlabel('Predicted Scores')
	plt.ylabel('Raw Residual')
	plt.savefig('../ACMManuscript/Figures/PredictedvsResiduals.png')
	# sys.exit()

	plt.scatter(np.array(rawGeneProfileSizes), np.array(residuals))
	plt.xlabel('Gene Profile Size')
	plt.ylabel('Raw Residual')
	plt.savefig('../ACMManuscript/Figures/RawResidualPlot_GeneSizes.png')

	plt.scatter(np.array(rawTaxonProfileSizes), np.array(residuals))
	plt.xlabel('Taxon Profile Size')
	plt.ylabel('Raw Residual')
	plt.savefig('../ACMManuscript/Figures/RawResidualPlot_TaxonSizes.png')


def sampleBins(scores, rawTaxonProfileSizes, rawGeneProfileSizes):
	logTaxonSize = [math.log(x) for x in rawTaxonProfileSizes]
	logGeneSize = [math.log(x) for x in rawGeneProfileSizes]
	taxonCount = dict()
	geneCount = dict()
	gridCount = dict()
	sample = dict()

	for i in range(0, len(logTaxonSize)):
		taxonRange = getRange(logTaxonSize[i], 1)
		geneRange = getRange(logGeneSize[i], 2)

		if taxonRange not in sample:
			sample[taxonRange] = dict()
		if geneRange not in sample[taxonRange]:
			sample[taxonRange][geneRange] = []
		sample[taxonRange][geneRange].append((logTaxonSize[i], logGeneSize[i], scores[i]))

		if taxonRange not in gridCount:
			gridCount[taxonRange] = dict()
		if geneRange not in gridCount[taxonRange]:
			gridCount[taxonRange][geneRange] = 0
		gridCount[taxonRange][geneRange] += 1

		n = 10000
		randomDataSet = selectRandom(sample, n)
		results = randomRegression(randomDataSet, "Random_ResidualPlot_TaxonSizes.png", "Random_ResidualPlot_GeneSizes.png")
		return results


def getRange(size, flag):
	scoreRange = 0
	if size >= 0 and size < 1:
		scoreRange = 1
	elif size >= 1 and size < 2:
		scoreRange = 2
	elif size >= 2 and size < 3:
		scoreRange = 3
	elif size >= 3 and size < 4:
		scoreRange = 4
	elif size >= 4 and size < 5:
		scoreRange = 5
	elif size >= 5 and size < 6:
		if flag == 2:
			scoreRange = 5
		else:
			scoreRange = 6
	elif size >= 6:
		scoreRange = 7

	return scoreRange


def selectRandom(sample, n):
	randomDataSet = []
	for taxonRange in sample:
		for geneRange in sample[taxonRange]:
			randomDataSet = randomDataSet + random.sample(sample[taxonRange][geneRange], n)
	return randomDataSet


def randomRegression(randomDataSet, taxonPlot, genePlot):
	taxonProfileSizes = []
	geneProfileSizes = []
	figureDir = "../ACMManuscript/Figures/"
	scores = []

	for i in range(0, len(randomDataSet)):
		taxonProfileSizes.append(randomDataSet[i][0])
		geneProfileSizes.append(randomDataSet[i][1])
		scores.append(randomDataSet[i][2])

	results = multiRegression(scores, taxonProfileSizes, geneProfileSizes)
	return results

	residuals = []
	predicted = results.fittedvalues
	for i in range(0, len(scores)):
		actualScore = scores[i]
		predictedScore = predicted[i]
		residual = actualScore - predictedScore
		residuals.append(residual)

	sizes = [geneProfileSizes]
	ones = np.ones(len(sizes[0]))
	X = sm.add_constant(np.column_stack((sizes[0], ones)))
	for ele in sizes[1:]:
		x = sm.add_constant(np.column_stack((ele, X)))
	results = sm.OLS(scores, X).fit()
	print "Random Gene Size to Residuals", results.summary()

	sizes = [taxonProfileSizes]
	ones = np.ones(len(sizes[0]))
	X = sm.add_constant(np.column_stack((sizes[0], ones)))
	for ele in sizes[1:]:
		X = sm.add_constant(np.column_stack((ele, X)))
	results = sm.OLS(scores, X).fit()
	print "Random Taxon Size to Residuals", results.summary()

	plt.scatter(np.array(geneProfileSizes), np.array(residuals), alpha=0.4)
	plt.xlabel('Gene Profile Size')
	plt.ylabel('Raw Residual')
	plt.savefig(figureDir + genePlot)

	plt.scatter(np.array(taxonProfileSizes), np.array(residuals), alpha=0.4)
	plt.xlabel('Taxon Profile Size')
	plt.ylabel('Raw Residual')
	plt.savefig(figureDir + taxonPlot)


def plotGrid(scores, rawTaxonProfileSizes, rawGeneProfileSizes):
	figureDir = "../ACMManuscript/Figures/"

	plt.subplot(221)
	plt.scatter(np.array(rawTaxonProfileSizes), np.array(scores), s=12)
	plt.xlabel('Taxon Size')
	plt.ylabel('Score')
	plt.grid(True)

	plt.subplot(222)
	plt.scatter(np.array(rawGeneProfileSizes), np.array(scores), s=12)
	plt.xlabel('Gene Size')
	plt.ylabel('Score')
	plt.grid(True)

	plt.savefig(figureDir + 'GridPlot.png')


if __name__=='__main__':
	from scipy.stats import norm
	import matplotlib.mlab as mlab
	import sys
	import random
	import os
	from statsmodels.stats.outliers_influence import OLSInfluence
	import math
	import matplotlib.pyplot as plt
	import numpy as np
	import statsmodels.api as sm
	import statsmodels.stats.api as sms
	from scipy import stats
	from scipy.stats import boxcox
	from statsmodels.stats import stattools
	main()
