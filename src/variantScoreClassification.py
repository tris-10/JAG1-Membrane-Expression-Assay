# Name:         variantScoreClassification.py
# Purpose:      Read in bin and variant control files and generate sort files containing variant scores
#
# Author:       Tristan J. Hayeck, Christopher J. Sottolano, and Melissa A. Gilbert
#
# Created:      2024
# Copyright:    (c) Tristan J. Hayeck 2024
# Licence:      <your licence>
# -------------------------------------------------------------------------------

import argparse
import datetime
import numpy as np
import pandas as pd
import labelsAndConstants as LABELS


def createCombinedBinTable(binFiles, depthScaleFactor = 100000, DEBUG=False):
    """
    Reads in the bin files and generates the combined bin table that adds in a combined ID
    variant column (ie chr20_10651694_G_A), the number of bins column 'bin',
    the adjusted alt allele count to deal with depth, and the frequency score
         Chr       pos    Depth...      variant             bin Adj_Alt_Allele_Count  FreqScore
    0  chr20  10651694  1585250...      chr20_10651694_G_A   1           217.631288   0.000329
    1  chr20  10651694  1585250...      chr20_10651694_G_C   1           110.708090   0.000167
    2  chr20  10651694  1585250...      chr20_10651694_G_T   1           114.492982   0.000173
    :param binFiles: Library containing mutant and wild-type bin files
    :param binWeights: Library containing weight coefficients assigned to mutant and wild-type bin files
    :param depthScaleFactor:
    :param DEBUG: Toggle debug
    :return: combinedBinsTable
    """
    binTables = []
    binIndex = 1
    for curBinFile in binFiles:
        curBinTable = pd.read_csv(curBinFile, sep=",")
        if DEBUG:
            print(f"Read in {curBinFile}")
            print(curBinTable.head())
        if not LABELS.ALT_COL in curBinTable.columns:
            print(f"Alt column is missing, ie no {LABELS.ALT_COL}, possibly it's Allele and should be renamed")
            raise ValueError

        curBinTable[LABELS.VARIANT_COL] = curBinTable[LABELS.CHROM_COL].astype(str) + "_" + \
                                          curBinTable[LABELS.POS_COL].astype(str) + "_" + \
                                          curBinTable[LABELS.REF_COL] + "_" + curBinTable[LABELS.ALT_COL]
        curBinTable[LABELS.BIN_INDEX_COL] = binIndex
        curBinTable[LABELS.ADJ_ALT_COUNT_COL] = depthScaleFactor * curBinTable[LABELS.ALT_COUNT_COL] / \
                                                curBinTable[LABELS.DEPTH_COL]
        curBinTable[LABELS.FREQ_SCORE_COL] = curBinTable[LABELS.ALT_COUNT_COL]/curBinTable[LABELS.ALT_COUNT_COL].sum()
        binIndex += 1
        binTables.append(curBinTable)
    combinedBinsTable = pd.concat(binTables)
    return combinedBinsTable


def constructPerVariantTable(combinedBinsTable, binWeights, DEBUG):
    """
    Calculated the weightedAve for each variant that will be used to get the variant score in the next steps.
    From the combined bins table and the bin weights generates the per variant table which is the variant, number bins,
    and weighted average. Looking something like this:
                     variant  numBinsVarintIn  weightedAve
            0  chr20_10651694_G_A              2.0     0.522629
            1  chr20_10651694_G_C              2.0     0.584930
            2  chr20_10651694_G_T              2.0     0.608740
            3  chr20_10652131_C_A              2.0     0.655050
            4  chr20_10652131_C_G              2.0     0.406153
    Later this will be merged with other data.
    :param combinedBinsTable: Combined mutant and wild-type variant bin table
    :param binWeights: Library containing weight coefficients assigned to mutant and wild-type bin files
    :return: perVariantTable
    """
    uniqueVariants = combinedBinsTable[LABELS.VARIANT_COL].unique()
    perVariantTable = pd.DataFrame(uniqueVariants, columns=[LABELS.VARIANT_COL])
    perVariantTable[LABELS.NUM_BINS_VAR_IN_COL] = np.nan
    perVariantTable[LABELS.WEIGHTED_AVE_COL] = np.nan
    # For each variant calculate the weighted average across bins, i.e.:
    # Wv = sum over bins w_b*F_bv / sum over bins F_bv
    for curVariant in uniqueVariants:
        curVariantTable = combinedBinsTable[combinedBinsTable[LABELS.VARIANT_COL] == curVariant]
        varIndex = perVariantTable[perVariantTable[LABELS.VARIANT_COL] == curVariant].index
        if varIndex.shape[0] != 1:
            print(f"ERROR: Multiple variants have the same ID {curVariant}")
            print(perVariantTable[perVariantTable[LABELS.VARIANT_COL] == curVariant])
            if curVariant is np.nan:
                print(curVariant)
                print("Multiple null values found. Check input.")
                break
            raise ValueError
        curVarIndex = varIndex[0]
        perVariantTable.at[curVarIndex, LABELS.NUM_BINS_VAR_IN_COL] = curVariantTable[LABELS.BIN_INDEX_COL].unique().shape[0]
        numerator = 0
        denomenator = 0
        for index, row in curVariantTable.iterrows():
            if row[LABELS.BIN_INDEX_COL] in binWeights:
                numerator += binWeights[row[LABELS.BIN_INDEX_COL]] * row[LABELS.FREQ_SCORE_COL]
                denomenator += row[LABELS.FREQ_SCORE_COL]
            else:
                print(f"WARNING: Cannot find bin index {row[LABELS.BIN_INDEX_COL]}.")
            if DEBUG:
                print("%s %s | num += %s * %s = %s" %
                      (curVariant, row[LABELS.BIN_INDEX_COL], binWeights[row[LABELS.BIN_INDEX_COL]],
                       row[LABELS.FREQ_SCORE_COL], numerator))
                print("\t | den += %s = %s " % (row[LABELS.FREQ_SCORE_COL], denomenator))
        perVariantTable.at[curVarIndex, LABELS.WEIGHTED_AVE_COL] = numerator/denomenator
    return perVariantTable


def formatAndFilterVariants(combinedBinsTable, perVariantTable):
    """
    Read in and combine bin files
    :param combinedBinsTable: Combined mutant and wild-type variant bin table
    :param perVariantTable: Variant table for individual replicate (sort)
    :return: perVariantMergeTable
    """
    # Extract specific columns from combinedBinsTable
    colsToGet = [LABELS.VARIANT_COL, LABELS.CHROM_COL, LABELS.POS_COL, LABELS.REF_COL, LABELS.ALT_COL,
                 LABELS.CODONS_COL, LABELS.AA_COL, LABELS.CONSEQ_COL, LABELS.PASS_COL, LABELS.AA_POSITION_COL,
                 LABELS.REF_CODON_COL, LABELS.NEW_CODON_COL, LABELS.CADD_PHRED_COL, LABELS.CADD_RAW_COL]
    combinedSummary = combinedBinsTable[colsToGet].groupby(LABELS.VARIANT_COL).first().reset_index()
    # Merge combinedBins and perVariant tables
    perVariantMergeTable = perVariantTable.merge(combinedSummary[colsToGet], on=LABELS.VARIANT_COL, how='left')
    # Check to see if combinedBins and perVariant tables are the same size
    if perVariantMergeTable.shape[0] != perVariantTable.shape[0]:
        print("Error merge issue %s != %s" % (perVariantMergeTable.shape[0], perVariantTable.shape[0]))
        raise ValueError
    # Filter out variants that do not meet PASS criteria
    print(f"There are {perVariantMergeTable[perVariantMergeTable[LABELS.PASS_COL].isna()].shape[0] } missing "
          f"{LABELS.PASS_CRITERIA} criteria")
    variantsBeforePassFilter = perVariantMergeTable[LABELS.VARIANT_COL].unique().shape[0]
    variantsAfterPassFilter = perVariantMergeTable[perVariantMergeTable[LABELS.PASS_COL] == LABELS.PASS_CRITERIA].shape[0]
    print(f"Before using pass criteria, there were {variantsBeforePassFilter}; after there are {variantsAfterPassFilter}")
    perVariantMergeTable = perVariantMergeTable[perVariantMergeTable[LABELS.PASS_COL] == LABELS.PASS_CRITERIA]
    return perVariantMergeTable


def formatTrainingVariants(perVariantMergeTable, knownMutFile, benignFile):
    """
    Read in and format the training data
    :param knownMutFile: Mutant variant control file
    :param benignFile: Wild-type variant control file
    :param DEBUG: Toggle debug
    :return: perVariantMergeTable, missingTable
    """
    # Read in known mutants file and assign variant IDs
    knownMuts = pd.read_csv(knownMutFile, sep="\t")
    knownMuts[LABELS.VARIANT_COL] = "chr"+knownMuts[LABELS.CHROM_COL].astype(str) + "_" + \
                                    knownMuts[LABELS.POS_COL].astype(str) + "_" + knownMuts[LABELS.REF_COL] + "_" + \
                                    knownMuts[LABELS.ALT_COL]
    numMuts = knownMuts.shape[0]
    # Assign mutant label to variants in perVariantMergeTable that are also in knownMuts
    perVariantMergeTable[LABELS.CLASSIFICATION_CAT_COL] = 0
    perVariantMergeTable.loc[perVariantMergeTable[LABELS.VARIANT_COL].isin(knownMuts[LABELS.VARIANT_COL].unique().tolist()),
                             LABELS.CLASSIFICATION_CAT_COL] = LABELS.MUT_CAT
    # Check to see if number of variants filtered perVariantMergeTable equals the number of variants in knownMuts
    numRunMuts = len(perVariantMergeTable[perVariantMergeTable[LABELS.CLASSIFICATION_CAT_COL] == LABELS.MUT_CAT])
    if numRunMuts != numMuts:
        # Retrieve IDs of variants missing from perVariantMergeTable
        missingSet = list(set(knownMuts[LABELS.VARIANT_COL].unique()) - set(perVariantMergeTable[LABELS.VARIANT_COL].unique()))
        print(f"WARNING: Only {numRunMuts} out of {knownMuts.shape[0]} known mutants found.")
        print(f"Mutants lost:")
        for var in missingSet:
            print(var)
        # Generate table for missing variants
        missingTable = pd.concat([pd.DataFrame(missingSet, columns=[LABELS.VARIANT_COL]),
                                  pd.DataFrame([LABELS.MUT_CAT]*len(missingSet), columns=[LABELS.CLASSIFICATION_CAT_COL])],
                                 axis=1)
    # Read in benign variants file and assign variant IDs
    benignTable = pd.read_csv(benignFile, sep='\t')
    benignTable[LABELS.VARIANT_COL] = "chr"+benignTable[LABELS.CHROM_COL].astype(str) + "_" + \
                                      benignTable[LABELS.POS_COL].astype(str) + "_" + benignTable[LABELS.REF_COL] + "_" + \
                                      benignTable[LABELS.ALT_COL]
    numBenign = benignTable.shape[0]
    # Assign benign label to variants in perVariantMergeTable that are also in benignTable
    benignVariants = benignTable[LABELS.VARIANT_COL].unique()
    perVariantMergeTable.loc[perVariantMergeTable[LABELS.VARIANT_COL].isin(benignVariants),
                             LABELS.CLASSIFICATION_CAT_COL] = LABELS.WT_CAT
    # Check to see if number of variants filtered perVariantMergeTable equals the number of variants in benignTable
    numRunBenign = len(perVariantMergeTable[perVariantMergeTable[LABELS.CLASSIFICATION_CAT_COL] == LABELS.WT_CAT])
    if numRunBenign != numBenign:
        # Retrieve IDs of variants missing from perVariantMergeTable
        missingSet = list(set(benignTable[LABELS.VARIANT_COL].unique()) - set(perVariantMergeTable[LABELS.VARIANT_COL].unique()))
        print(f"WARNING: Only {numRunBenign} out of {benignTable.shape[0]} known benign variants found.")
        print(f"Benign variants lost:")
        for var in missingSet:
            print(var)
        # Generate table for missing variants
        missingTable = pd.concat([missingTable, pd.concat([pd.DataFrame(missingSet, columns=[LABELS.VARIANT_COL]),
                                  pd.DataFrame([LABELS.WT_CAT]*len(missingSet), columns=[LABELS.CLASSIFICATION_CAT_COL])],
                                 axis=1)], axis=0)
    return perVariantMergeTable, missingTable


def calcVariantScores(perVariantMergeTable, wtCalibrationSet, mutantCalibrationSet):
    """
    Calculate variant scores
    :param perVariantMergeTable: Variant table for individual replicate (sort)
    :param wtCalibrationSet: Wild-type calibration set (variant IDs seperated by commas)
    :param mutantCalibrationSet: Mutant calibration set (variant IDs seperated by commas)
    :return: perVariantMergeTable
    """
    wtCalibrationSet = wtCalibrationSet.split(',')
    # Extract rows for select benign calibration set
    wtCalibrationSetTable = perVariantMergeTable[perVariantMergeTable[LABELS.VARIANT_COL].isin(wtCalibrationSet)]
    # Check whether chosen benign calibration set is available in perVariantMergeTable
    if wtCalibrationSetTable.shape[0] != len(wtCalibrationSet):
        print("Incorrect number of calibration variants found")
        print(wtCalibrationSetTable)
        raise ValueError

    # Get weighted average mean of benign calibration variants
    wildTypeVarScoreMean = wtCalibrationSetTable[LABELS.WEIGHTED_AVE_COL].mean()
    # Extract rows for select mutant calibration set
    mutantCalibrationSet = mutantCalibrationSet.split(',')
    mutantCalibrationSetTable = perVariantMergeTable[perVariantMergeTable[LABELS.VARIANT_COL].isin(mutantCalibrationSet)]
    # Check whether chosen mutant calibration set is available in perVariantMergeTable
    if mutantCalibrationSetTable.shape[0] != len(mutantCalibrationSet):
        print("Incorrect number of mutant calibration variants found")
        print(mutantCalibrationSetTable)
        raise ValueError

    # Get weighted average mean of mutant calibration variants
    mutantVarScoreMean = mutantCalibrationSetTable[LABELS.WEIGHTED_AVE_COL].mean()
    # Calculate variant score
    perVariantMergeTable[LABELS.VARIANT_SCORE_COL] = (perVariantMergeTable[LABELS.WEIGHTED_AVE_COL] - mutantVarScoreMean) / \
                                                     (wildTypeVarScoreMean - mutantVarScoreMean)
    return perVariantMergeTable


def run(bin1File, bin2File, binWeight1, binWeight2, knownMutFile, benignFile, outputTable, wtCalibrationSet,
        mutantCalibrationSet, DEBUG):
    """
    Read in bin and variant control files and generate sort files containing variant scores
    :param bin1File: Mutant bin file
    :param bin2File: Wild-type bin file
    :param binWeight1: Mutant bin weight
    :param binWeight2: Wild-type bin weight
    :param knownMutFile: Mutant variant control file
    :param benignFile: Wild-type variant control file
    :param outputTable: Output directory
    :param wtCalibrationSet: Wild-type calibration set (variant IDs seperated by commas)
    :param mutantCalibrationSet: Mutant calibration set (variant IDs seperated by commas)
    :param DEBUG: Toggle debug
    :return: (None)
    """
    # Generate bin file list and weight dictionary
    binFiles=[bin1File, bin2File]
    binWeights = {1: binWeight1, 2: binWeight2}
    # Generate combined bin table
    combinedBinsTable = createCombinedBinTable(binFiles, DEBUG)
    # Calculate weighted averages for each variant
    perVariantTable = constructPerVariantTable(combinedBinsTable, binWeights, DEBUG)
    perVariantMergeTable = formatAndFilterVariants(combinedBinsTable, perVariantTable)
    # Set the categories
    perVariantMergeTable, missingTable = formatTrainingVariants(perVariantMergeTable, knownMutFile, benignFile)
    # Calculate variant scores
    perVariantMergeTable = calcVariantScores(perVariantMergeTable, wtCalibrationSet, mutantCalibrationSet)
    # Output variant table
    perVariantMergeTable.to_csv(outputTable+LABELS.TSV_EXT, sep="\t", index=False)
    # Output missing variant table
    missingTable.to_csv(outputTable+f"_missingVar{LABELS.LOG_EXT}", sep="\t", index=False)
    return

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("-bm", "--bin1File", type=str,
                        help="Input bin file containing mutant allele counts")
    parser.add_argument("-bw", "--bin2File", type=str,
                        help="Input bin file containing wild-type allele counts")
    parser.add_argument("-m", "--knownMutFile", type=str,
                        help="Mutant variant control file")
    parser.add_argument("-w", "--benignFile", type=str,
                        help="Wild-type variant control file")
    parser.add_argument("-o", "--outputTable", type=str,
                        help="Output directory")
    parser.add_argument("-w1", "--binWeight1", type=float,
                        default=0.1,
                        help="Scalar weight applied mutant variant scores")
    parser.add_argument("-w2", "--binWeight2", type=float,
                        default=0.9,
                        help="Scalar weight applied wild-type variant scores")
    parser.add_argument("-mc", "--mutantCalibrationSet", type=str,
                        help="Variants used for mutant calibration set (variant IDs separated by commas)")
    parser.add_argument("-wc", "--wtCalibrationSet", type=str,
                        help="Variants used for wild-type calibration set (variant IDs separated by commas)")
    parser.add_argument('--debug', dest='DEBUG', action='store_true')
    parser.set_defaults(DEBUG=False)

    start = datetime.datetime.now()
    args = parser.parse_args()
    print("Starting processing %s" % start)
    print(args)

    run(args.bin1File, args.bin2File, args.binWeight1, args.binWeight2, args.knownMutFile, args.benignFile,
        args.outputTable, args.wtCalibrationSet, args.mutantCalibrationSet, args.DEBUG)

    done = datetime.datetime.now()
    elapsed = done - start
    duration = ':'.join(str(elapsed).split(':')[1:])
    print("The duration was: %s (%s seconds)" % (duration, elapsed.total_seconds()))
    print("Finished processing %s" % done)


if __name__ == '__main__':
    main()
