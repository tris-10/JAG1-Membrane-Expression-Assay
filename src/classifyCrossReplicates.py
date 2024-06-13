# Name:         classifyCrossReplicates.py
# Purpose:      Read in and process perVariantScore tables across replicates, plot, and generate classification thresholds
#
# Author:       Tristan J. Hayeck, Christopher J. Sottolano, and Melissa A. Gilbert
#
# Created:      2024
# Copyright:    (c) Tristan J. Hayeck 2024
# Licence:     <your licence>
# -------------------------------------------------------------------------------

import argparse
import datetime
import numpy as np
import pandas as pd
import scipy.stats as stats
import labelsAndConstants as LABELS
import plotter


def createVariantCol(inputTable, appendChr=True):
    """
    Append "chr" if not present in header
    :param inputTable: Variant control tables
    :param appendChr: (bool) Append chr prefix to variant IDs
    :return: inputTable
    """
    if appendChr:
        inputTable[LABELS.VARIANT_COL] = "chr"+inputTable[LABELS.CHROM_COL].astype(str) + "_" + inputTable[LABELS.POS_COL].astype(str) +\
                                   "_" + inputTable[LABELS.REF_COL] + "_" + inputTable[LABELS.ALT_COL]
    else:
        inputTable[LABELS.VARIANT_COL] = inputTable[LABELS.CHROM_COL].astype(str) + "_" + inputTable[LABELS.POS_COL].astype(str) +\
                                   "_" + inputTable[LABELS.REF_COL] + "_" + inputTable[LABELS.ALT_COL]


def combineReplicateFiles(replicateTableFile, numSorts, DEBUG):
    """
    Combine all replicate perVariantScore files into combineReplicates table
    :param replicateTableFile: Generic name for the different replicates, expecting to put the sort/replicate
    number in the string twice, ie replicateTableFile % curSort
    :return: combineReplicates
    """
    combineReplicates = None
    identifierCols = [LABELS.VARIANT_COL, LABELS.CHROM_COL, LABELS.POS_COL, LABELS.REF_COL, LABELS.ALT_COL,
                      LABELS.AA_COL, LABELS.CONSEQ_COL, LABELS.AA_POSITION_COL, LABELS.REF_CODON_COL, LABELS.NEW_CODON_COL,
                      LABELS.CADD_PHRED_COL, LABELS.CADD_RAW_COL]
    for curSort in range(1, numSorts+1):
        if DEBUG:
            print(curSort)
        repTable = pd.read_csv(replicateTableFile % curSort, sep="\t")
        repTable = repTable[identifierCols + [LABELS.VARIANT_SCORE_COL]]
        repTable.rename(columns={LABELS.VARIANT_SCORE_COL: LABELS.VARIANT_SCORE_COL+f"_{curSort}"}, inplace=True)
        if combineReplicates is None:
            combineReplicates = repTable
        else:
            combineReplicates = pd.merge(combineReplicates, repTable, on=identifierCols, how='outer')
    return combineReplicates


def filterMinVariants(combineReplicates, replicateScoreCols, logFile, numSorts, minNumReplicates, DEBUG):
    """
    Find and remove variants with less than minimum number of replicates
    :param combineReplicates: Variant table containing scores for each replicate
    :param replicateScoreCols: Replicate variant score column names
    :param logFile: Log containing filtered/missing variants
    :param numSorts: Number of replicate perVariantScore (sort) files
    :param minNumReplicates: Minimum number of replicates per variant allowed above threshold filter
    :param DEBUG: Toggle debug
    :return: combineReplicatesCleaned, logFile
    """
    combineReplicates[LABELS.NULL_COUNT_COL] = (combineReplicates[replicateScoreCols]).apply(lambda row: row.isnull().sum(), axis=1)
    setToRemove = combineReplicates[combineReplicates[LABELS.NULL_COUNT_COL] > numSorts - minNumReplicates]
    combineReplicatesCleaned = combineReplicates[combineReplicates[LABELS.NULL_COUNT_COL] <= numSorts - minNumReplicates]
    logOut = f"Warning: removing {setToRemove.shape[0]} variant(s) with only 1 replicate"
    print(logOut)
    logFile = pd.concat([logFile, pd.DataFrame([logOut])])
    if DEBUG:
        print(setToRemove[LABELS.VARIANT_COL])
    return combineReplicatesCleaned, logFile


def filterNonsenseVariants(combineReplicatesCleaned, logFile, DEBUG):
    """
    Filter and log all nonsense variants
    :param combineReplicatesCleaned: Filtered variant table containing scores for each replicate
    :param logFile: Log containing filtered/missing variants
    :param DEBUG: Toggle debug
    :return: combineReplicatesCleaned, logFile
    """
    # Remove nonsense variants
    nonsense_var = combineReplicatesCleaned[combineReplicatesCleaned[LABELS.CONSEQ_COL].str.contains(LABELS.NONSENSE_VAR)]
    combineReplicatesCleaned = combineReplicatesCleaned[~combineReplicatesCleaned[LABELS.CONSEQ_COL].str.contains(LABELS.NONSENSE_VAR)]
    logOut = f"Warning: removing {nonsense_var.shape[0]} nonsense variants"
    print(logOut)
    logFile = pd.concat([logFile, pd.DataFrame([logOut])])
    if DEBUG:
        print(nonsense_var[LABELS.VARIANT_COL])
    return combineReplicatesCleaned, logFile


def filterVariantsByDatabase(pathTable, benignTable, combineReplicatesCleaned, logFile, dbBenign, dbPath):
    """
    Filter variants by database
    :param pathTruthFile: Mutant variant control file
    :param benignTruthFile: Wild-type variant control file
    :param combineReplicatesCleaned: Filtered variant table containing scores for each replicate
    :param logFile: Log containing filtered/missing variants
    :param dbBenign: Wild-type variants databases used in analysis (Default: "" [all variants])
    :param dbPath: Mutant variants databases used in analysis (Default: "" [all variants])
    :return: combineReplicatesCleaned, pathTable, benignTable, logFile
    """
    if (dbBenign != "") | (dbPath != ""):
        logOut = f"Wild-type variant database(s): {dbBenign}, mutant variant database(s): {dbPath}"
        print(logOut)
        logFile = pd.concat([logFile, pd.DataFrame([logOut])])
        dbPath = dbPath.split(',')
        dbBenign = dbBenign.split(',')
        pathMask = pathTable[LABELS.SOURCE_COL] == dbPath[0]
        benignMask = benignTable[LABELS.DATABASE_COL] == dbBenign[0]
        for db in dbPath:
            pathMask = pathMask.mask(pathTable[LABELS.SOURCE_COL] == db, True)
        for db in dbBenign:
            benignMask = benignMask.mask(benignTable[LABELS.DATABASE_COL] == db, True)
        pathTable = pathTable[pathMask]
        benignTable = benignTable[benignMask]
    return combineReplicatesCleaned, pathTable, benignTable, logFile


def calcStats(combineReplicatesCleaned, replicateScoreCols):
    """
    Calculate mean, standard dev, and confident intervals using variant scores
    :param combineReplicatesCleaned: Filtered variant table containing scores for each replicate
    :param replicateScoreCols: Replicate variant score column names
    :return: combineReplicatesCleaned
    """
    combineReplicatesCleaned[LABELS.MEAN_ACROSS_REPS_COL] = combineReplicatesCleaned[replicateScoreCols].mean(axis=1)
    combineReplicatesCleaned[LABELS.STD_ACROSS_REPS_COL] = combineReplicatesCleaned[replicateScoreCols].std(axis=1)

    # Calculate confidence intervals, ie mean+/-z*std/sqrt(#valid replicates at variant)
    combineReplicatesCleaned[LABELS.UPPER_CONF_INT_COL] = combineReplicatesCleaned[LABELS.MEAN_ACROSS_REPS_COL] + \
                                                        stats.norm.ppf(0.95) * \
                                                        combineReplicatesCleaned[LABELS.STD_ACROSS_REPS_COL] / \
                                                        combineReplicatesCleaned[replicateScoreCols].count(axis=1).apply(np.sqrt)
    combineReplicatesCleaned[LABELS.LOWER_CONF_INT_COL] = combineReplicatesCleaned[LABELS.MEAN_ACROSS_REPS_COL] - \
                                                        stats.norm.ppf(0.95) * \
                                                        combineReplicatesCleaned[LABELS.STD_ACROSS_REPS_COL] / \
                                                        combineReplicatesCleaned[replicateScoreCols].count(axis=1).apply(np.sqrt)
    return combineReplicatesCleaned


def mergeAndFilterTables(pathTruthFile, benignTruthFile, replicateTableFile, outputPrefix, numSorts, minNumReplicates,
                         dbBenign, dbPath, DEBUG):
    """
    Filter variants and merge tables
    :param pathTruthFile: Mutant variant control file
    :param benignTruthFile: Wild-type variant control file
    :param outputPrefix: Output directory
    :param numSorts: Number of replicate perVariantScore (sort) files
    :param minNumReplicates: Minimum number of replicates per variant allowed above threshold filter
    :param dbBenign: Wild-type variants databases used in analysis (Default: "" [all variants])
    :param dbPath: Mutant variants databases used in analysis (Default: "" [all variants])
    :param DEBUG: Toggle debug
    :return: combineReplicatesCleaned, logFile, pathSet, benignSet, librarySet, numBenign, numPath
    """
    # Import truth sets
    pathTable = pd.read_csv(pathTruthFile, sep="\t")
    benignTable = pd.read_csv(benignTruthFile, sep="\t")
    createVariantCol(benignTable)
    createVariantCol(pathTable)
    logFile = pd.DataFrame()

    # Log input data in log file
    inpLogs = [f"Benign variant control file: {benignTruthFile}",
               f"Mutant variant control file: {pathTruthFile}",
               f"Output file prefix: {outputPrefix}"]
    for line in inpLogs:
        logFile = pd.concat([logFile, pd.DataFrame([line])])

    # Import sort files for different replicates and append them to table
    combineReplicates = combineReplicateFiles(replicateTableFile, numSorts, DEBUG)
    replicateScoreCols = combineReplicates.columns[pd.Series(combineReplicates.columns).str.startswith(LABELS.VARIANT_SCORE_COL)]

    # Filter variants
    combineReplicatesCleaned, logFile = filterMinVariants(combineReplicates, replicateScoreCols, logFile, numSorts, minNumReplicates, DEBUG)
    combineReplicatesCleaned, logFile = filterNonsenseVariants(combineReplicatesCleaned, logFile, DEBUG)
    combineReplicatesCleaned, pathTable, benignTable, logFile = filterVariantsByDatabase(pathTable, benignTable, combineReplicatesCleaned, logFile, dbBenign, dbPath)

    # Calculate mean, standard dev, and confident intervals
    combineReplicatesCleaned = calcStats(combineReplicatesCleaned, replicateScoreCols)

    # Add labels from both benign and mutant truth sets
    combineReplicatesCleaned[LABELS.CLASSIFICATION_CAT_COL] = LABELS.LIBRARY_CAT
    pathVar = pathTable[LABELS.VARIANT_COL].unique().tolist()
    combineReplicatesCleaned.loc[combineReplicatesCleaned[LABELS.VARIANT_COL].isin(pathVar), LABELS.CLASSIFICATION_CAT_COL] = LABELS.MUT_CAT
    benignVar = benignTable[LABELS.VARIANT_COL].unique().tolist()
    combineReplicatesCleaned.loc[combineReplicatesCleaned[LABELS.VARIANT_COL].isin(benignVar), LABELS.CLASSIFICATION_CAT_COL] = LABELS.WT_CAT

    # Check variant counts after filtering
    numBenign = (combineReplicatesCleaned[LABELS.CLASSIFICATION_CAT_COL] == LABELS.WT_CAT).sum()
    numPath = (combineReplicatesCleaned[LABELS.CLASSIFICATION_CAT_COL] == LABELS.MUT_CAT).sum()
    logOut = f"Pathogenic file variants: {len(pathVar)}, pathogenic variants after filtering: {numPath}"
    print(logOut)
    logFile = pd.concat([logFile, pd.DataFrame([logOut])])
    logOut = f"Benign file variants: {len(benignVar)}, benign variants after filtering: {numBenign}"
    print(logOut)
    logFile = pd.concat([logFile, pd.DataFrame([logOut])])

    # Generate tables for pathogenic, benign, and library sets
    pathSet = combineReplicatesCleaned[combineReplicatesCleaned[LABELS.CLASSIFICATION_CAT_COL] == LABELS.MUT_CAT]
    benignSet = combineReplicatesCleaned[combineReplicatesCleaned[LABELS.CLASSIFICATION_CAT_COL] == LABELS.WT_CAT]
    librarySet = combineReplicatesCleaned[combineReplicatesCleaned[LABELS.CLASSIFICATION_CAT_COL] == LABELS.LIBRARY_CAT]
    return combineReplicatesCleaned, logFile, pathSet, benignSet, librarySet, numBenign, numPath


def calcThresholds(controlSet, confInt, varString, logFile):
    """
    Calculate thresholds based on variant control set means
    :param controlSet: Variant control table
    :param confInt: (constant) Upper or lower confidence interval label
    :param varString: (string) "Mutant" or "Benign"
    :param logFile: Log containing filtered/missing variants
    :return: controlMean, controlCI, controlLowerPerc, controlUpperPerc, logFile
    """
    controlMean = controlSet[LABELS.MEAN_ACROSS_REPS_COL].mean()
    controlCI = controlSet[confInt].mean()
    controlLowerPerc = np.percentile(controlSet[LABELS.MEAN_ACROSS_REPS_COL].values, 5)
    controlUpperPerc = np.percentile(controlSet[LABELS.MEAN_ACROSS_REPS_COL].values, 95)
    logOut = f"{varString} means (95% percentile threshold) {controlMean} ({controlUpperPerc})"
    print(logOut)
    logFile = pd.concat([logFile, pd.DataFrame([logOut])])
    return controlMean, controlCI, controlLowerPerc, controlUpperPerc, logFile


def calcOddsPath(numPath, numBenign, numTruePredPath, numPredPath, numTruePathPredBen, logFile):
    """
    Calculate odds path
    :param numPath: Number of mutant control variants after filtering
    :param numBenign: Number of wild-type control variants after filtering
    :param numTruePredPath: Number of correctly predicted pathogenic variants
    :param numPredPath: Number of predicted pathogenic variants
    :param numTruePathPredBen: Number correctly predicted benign variants
    :param logFile: Log containing filtered/missing variants
    :return: logFile
    """
    # P1 = # true pathogenic / (true pathogenic + true benign)
    p1 = numPath / (numPath + numBenign)
    # if numTruePredPath == numPredPath, add correction factor (+1)
    if numPredPath == numTruePredPath:
        p2 = numTruePredPath / (numPredPath + 1)
    else:
        p2 = numTruePredPath / numPredPath
    p2norm = numTruePathPredBen / numBenign
    # Calculate odds path
    oddsPath = (p2 * (1 - p1)) / ((1 - p2) * p1)
    oddsPathnorm = (p2norm * (1 - p1)) / ((1 - p2norm) * p1)
    logOut = f"OddsPath = {oddsPath}"
    print(logOut)
    logFile = pd.concat([logFile, pd.DataFrame([logOut])])
    logOut = f"OddsPathnorm = {oddsPathnorm}"
    print(logOut)
    logFile = pd.concat([logFile, pd.DataFrame([logOut])])
    return logFile


def classifyVariants(combineReplicatesCleaned, lowerPercBenign, upperPercPath, includeLikely, numPath, numBenign, logFile):
    """
    Classify variants by confidence interval
    :param combineReplicatesCleaned: Filtered variant table containing scores for each replicate
    :param lowerPercBenign: Lower threshold (5th percentile of benign variants)
    :param upperPercPath: Upper threshold (95th percentile of mutant variants)
    :param includeLikely: Include likely pathogenic variants in OddsPath score
    :param numPath: Number of mutant control variants after filtering
    :param numBenign: Number of wild-type control variants after filtering
    :param logFile: Log containing filtered/missing variants
    :return: combineReplicatesCleaned, logFile
    """
    combineReplicatesCleaned[LABELS.CLASSIFICATION_PERC_COL] = np.nan
    varClass = {LABELS.AMBIG_CAT: combineReplicatesCleaned[LABELS.MEAN_ACROSS_REPS_COL].between(lowerPercBenign, upperPercPath, inclusive="both"),
                LABELS.POSS_PATH_CAT: (combineReplicatesCleaned[LABELS.MEAN_ACROSS_REPS_COL] < lowerPercBenign) &
                                      (combineReplicatesCleaned[LABELS.UPPER_CONF_INT_COL] >= lowerPercBenign),
                LABELS.POSS_BENIGN_CAT: (combineReplicatesCleaned[LABELS.MEAN_ACROSS_REPS_COL] > upperPercPath) &
                                      (combineReplicatesCleaned[LABELS.LOWER_CONF_INT_COL] <= upperPercPath),
                LABELS.PATH_CAT: combineReplicatesCleaned[LABELS.UPPER_CONF_INT_COL] < lowerPercBenign,
                LABELS.BENIGN_CAT: combineReplicatesCleaned[LABELS.LOWER_CONF_INT_COL] > upperPercPath}

    for cat in varClass:
        combineReplicatesCleaned[LABELS.CLASSIFICATION_PERC_COL].mask(varClass[cat], cat, inplace=True)
    if includeLikely:
        # Include pathogenic and likely pathogenic
        numTruePredPath = combineReplicatesCleaned[(combineReplicatesCleaned[LABELS.CLASSIFICATION_CAT_COL] == LABELS.MUT_CAT) &
                                                   ((combineReplicatesCleaned[LABELS.CLASSIFICATION_PERC_COL] == LABELS.PATH_CAT) |
                                                    (combineReplicatesCleaned[LABELS.CLASSIFICATION_PERC_COL] == LABELS.POSS_PATH_CAT))].shape[0]
        numPredPath = combineReplicatesCleaned[(combineReplicatesCleaned[LABELS.CLASSIFICATION_CAT_COL] != LABELS.LIBRARY_CAT) &
                                               ((combineReplicatesCleaned[LABELS.CLASSIFICATION_PERC_COL] == LABELS.PATH_CAT) |
                                                (combineReplicatesCleaned[LABELS.CLASSIFICATION_PERC_COL] == LABELS.POSS_PATH_CAT))].shape[0]
    else:
        # Include pathogenic only
        numTruePredPath = combineReplicatesCleaned[(combineReplicatesCleaned[LABELS.CLASSIFICATION_CAT_COL] == LABELS.MUT_CAT) &
                                                   (combineReplicatesCleaned[LABELS.CLASSIFICATION_PERC_COL] == LABELS.PATH_CAT)].shape[0]
        numPredPath = combineReplicatesCleaned[(combineReplicatesCleaned[LABELS.CLASSIFICATION_CAT_COL] != LABELS.LIBRARY_CAT) &
                                               (combineReplicatesCleaned[LABELS.CLASSIFICATION_PERC_COL] == LABELS.PATH_CAT)].shape[0]

    numTruePathPredBen = combineReplicatesCleaned[(combineReplicatesCleaned[LABELS.CLASSIFICATION_CAT_COL] == LABELS.MUT_CAT) &
                                                  (combineReplicatesCleaned[LABELS.CLASSIFICATION_PERC_COL] == LABELS.BENIGN_CAT)].shape[0]
    logFile = calcOddsPath(numPath, numBenign, numTruePredPath, numPredPath, numTruePathPredBen, logFile)
    return combineReplicatesCleaned, logFile


def run(pathTruthFile, benignTruthFile, replicateTableFile, outputPrefix, numSorts, minNumReplicates, dbBenign, dbPath,
        labels, colors, N_bin, fixedYAxis, yMaxOffset, includeLikely, plotAllHist, DEBUG):
    """
    Read in and process perVariantScore tables across replicates, plot, and generate classification thresholds
    :param pathTruthFile: Mutant variant control file
    :param benignTruthFile: Wild-type variant control file
    :param replicateTableFile: Replicate perVariantScore file (substitute replicate number with %s)
    :param outputPrefix: Output directory
    :param numSorts: Number of replicate perVariantScore (sort) files
    :param minNumReplicates: Minimum number of replicates per variant allowed above threshold filter
    :param dbBenign: Wild-type variants databases used in analysis (Default: "" [all variants])
    :param dbPath: Mutant variants databases used in analysis (Default: "" [all variants])
    :param labels: Plot legend labels (Separate by commas)
    :param colors: Plot colors (Separate by commas)
    :param N_bin: Number of bins used for plotting histograms
    :param fixedYAxis: (bool) If True, sets y-axis as fixed the maximum histogram value
    :param yMaxOffset: Scalar multiplier of y-axis
    :param includeLikely: Include likely pathogenic variants in OddsPath score
    :param plotAllHist: Plot all variants in single histogram (requires only one color/labels)
    :param DEBUG: Toggle debug
    :return: (None)
    """
    # Merge and filter variants
    combineReplicatesCleaned, logFile, pathSet, \
    benignSet, librarySet, numBenign, numPath = mergeAndFilterTables(pathTruthFile, benignTruthFile, replicateTableFile, outputPrefix,
                                                                     numSorts, minNumReplicates, dbBenign, dbPath, DEBUG)
    # Calculate benign threshold parameters
    benignMean, benignLwrCI, lowerPercBenign, upperPercBenign, logFile = calcThresholds(benignSet, LABELS.LOWER_CONF_INT_COL,
                                                                                        LABELS.PATH_STRING, logFile)
    # Calculate pathogenic threshold parameters
    pathMean, pathUpCI, lowerPercPath, upperPercPath, logFile = calcThresholds(pathSet, LABELS.UPPER_CONF_INT_COL,
                                                                               LABELS.BENIGN_STRING, logFile)
    # Generate strings for labels/plots
    likelyString = ""
    filtDbString = ""
    if includeLikely:
        likelyString = "_likely"
    if (dbBenign != "") | (dbPath != ""):
        filtDbString = "_filtDb"
    logFn = outputPrefix+likelyString+filtDbString+LABELS.LOG_EXT
    outputFile = outputPrefix+likelyString+filtDbString+LABELS.TSV_EXT
    histSepFile = outputPrefix+likelyString+filtDbString+LABELS.PDF_EXT

    # Plot data
    labels = labels.split(',')
    colors = colors.split(',')
    if plotAllHist:
        # All data
        dataRef = []
        datalib = [combineReplicatesCleaned[LABELS.MEAN_ACROSS_REPS_COL]]
    else:
        # Truth and library set
        dataRef = [benignSet[LABELS.MEAN_ACROSS_REPS_COL], pathSet[LABELS.MEAN_ACROSS_REPS_COL]]
        datalib = [librarySet[LABELS.MEAN_ACROSS_REPS_COL]]
    plotter.plotSeperateHist(dataRef, datalib, labels, colors, N_bin, lowerPercBenign, histSepFile, fixedYAxis, yMaxOffset)
    # Classify variants
    combineReplicatesCleaned, logFile = classifyVariants(combineReplicatesCleaned, lowerPercBenign, upperPercPath,
                                                         includeLikely, numPath, numBenign, logFile)
    combineReplicatesCleaned.to_csv(outputFile, sep="\t", index=False)
    logFile.to_csv(logFn, index=False, header=None)
    return


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", "--pathTruthFile", type=str,
                        help="Mutant variant control file")
    parser.add_argument("-w", "--benignTruthFile", type=str,
                        help="Wild-type variant control file")
    parser.add_argument("-o", "--outputPrefix", type=str,
                        help="Output directory")
    parser.add_argument("-r", "--replicateTableFile", type=str,
                        help="Replicate perVariantScore file (substitute replicate number with %s)")
    parser.add_argument("-s", "--numSorts", type=int,
                        help="Number of replicate perVariantScore (sort) files")
    parser.add_argument("-n", "--minNumReplicates", type=int,
                        help="Minimum number of replicates per variant allowed above threshold filter")
    parser.add_argument("-db", "--dbBenign", type=str,
                        default="",
                        help="Wild-type variants databases used in analysis (Default: "" [all variants])")
    parser.add_argument("-dm", "--dbPath", type=str,
                        default="",
                        help="Mutant variants databases used in analysis (Default: "" [all variants])")
    parser.add_argument("-pl", "--labels", type=str,
                        help="Plot legend labels (Separate by commas)")
    parser.add_argument("-pc", "--colors", type=str,
                        help="Plot colors (Separate by commas)")
    parser.add_argument("-b", "--N_bin", type=int,
                        default=50,
                        help="Number of bins used for plotting histograms")
    parser.add_argument("-y", "--yMaxOffset", type=float,
                        default=1.2,
                        help="Scalar multiplier of the histogram y-axis")
    parser.add_argument('--fixedYAxis', dest='fixedYAxis', action='store_true',
                        help="Sets Y-axis as fixed the maximum histogram value")
    parser.set_defaults(fixedYAxis=False)
    parser.add_argument('--includeLikely', dest='includeLikely', action='store_true',
                        help="Include likely pathogenic variants in OddsPath score")
    parser.set_defaults(includeLikely=False)
    parser.add_argument('--plotAllHist', dest='plotAllHist', action='store_true',
                        help="Plot all variants in single histogram (requires only one color/labels)")
    parser.set_defaults(plotAllHist=False)
    parser.add_argument('--debug', dest='DEBUG', action='store_true')
    parser.set_defaults(DEBUG=False)

    start = datetime.datetime.now()
    args = parser.parse_args()
    print("Starting processing %s" % start)
    print(args)

    run(args.pathTruthFile, args.benignTruthFile, args.replicateTableFile, args.outputPrefix, args.numSorts,
        args.minNumReplicates, args.dbBenign, args.dbPath, args.labels, args.colors, args.N_bin, args.fixedYAxis,
        args.yMaxOffset, args.includeLikely, args.plotAllHist, args.DEBUG)

    done = datetime.datetime.now()
    elapsed = done - start
    duration = ':'.join(str(elapsed).split(':')[1:])
    print("The duration was: %s (%s seconds)" % (duration, elapsed.total_seconds()))
    print("Finished processing %s" % done)


if __name__ == '__main__':
    main()
