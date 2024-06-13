# Bin file columns and criteria
CHROM_COL = "Chr"
POS_COL = "pos"
DEPTH_COL = "Depth"
REF_COL = "Ref_Allele"
REF_COUNT_COL = "Ref_Allele_Count"
ALT_COL = "Alt_Allele"
ALT_COUNT_COL = "Alt_Allele_Count"
AA_POSITION_COL = "AA_Position"
CODONS_COL = "Codons"
AA_COL = "Amino_acids"
CONSEQ_COL = "Consequence"
PASS_COL = "Twist_PassFail"
PASS_CRITERIA = "PASS"
CADD_PHRED_COL = "CADD_PHRED"
CADD_RAW_COL = "CADD_RAW"

# Consequence labels
NONSENSE_VAR = "stop_gained"
SYNON_CONSQ = 'synonymous_variant'
SYN_OR_SPLICE_CONSQ = "splice_region_variant,synonymous_variant"
MISSENSE_CONSQ = 'missense_variant'

# Variant control file columns
DATABASE_COL = "Database"
SOURCE_COL = "Source"

# Generated columns for calculation of variant score
ADJ_ALT_COUNT_COL = "Adj_Alt_Allele_Count"
VARIANT_COL = "variant"
BIN_INDEX_COL = "bin"
FREQ_SCORE_COL = "FreqScore"
NUM_BINS_VAR_IN_COL = "numBinsVarintIn"
WEIGHTED_AVE_COL = "weightedAve"
VARIANT_SCORE_COL = "VarScore"
FINAL_VARIANT_SCORE_COL = "FinalVarScore"

# Generated columns for combining across replicates
NULL_COUNT_COL = "nullCount"
MEAN_ACROSS_REPS_COL = "meanAcrossReps"
STD_ACROSS_REPS_COL = "stdAcrossReps"
UPPER_CONF_INT_COL = "upperConfidenceInterval"
LOWER_CONF_INT_COL = "lowerConfidenceInterval"
REF_CODON_COL = "ref_codon"
ALT_CODON_COL = "alt_codon"
NEW_CODON_COL = "new_codon"

# Categories codes and strings
CLASSIFICATION_CAT_COL = "Category"
MUT_CAT = 1
WT_CAT = -1
LIBRARY_CAT = 0
CLASSIFICATION_PERC_COL = "Classification"
PATH_CAT = -2
POSS_PATH_CAT = -1
AMBIG_CAT = 0
POSS_BENIGN_CAT = 1
BENIGN_CAT = 2
classicationCatToLabel = {MUT_CAT: "Mutant", WT_CAT: "WildType", LIBRARY_CAT: "Library"}
classicationCatToColor = {MUT_CAT: "red", WT_CAT: "blue", LIBRARY_CAT: "grey"}
PATH_STRING = "Pathogenic"
BENIGN_STRING = "Benign"

#File extensions
TSV_EXT = ".tsv"
LOG_EXT = ".log"
PDF_EXT = ".pdf"
