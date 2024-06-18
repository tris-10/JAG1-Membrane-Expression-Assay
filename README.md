# JAG1-Membrane-Expression-Assay
## Functional characterization of 2,832 JAG1 variants supports reclassification for Alagille syndrome and improves guidance for clinical variant interpretation

*Authors:* Melissa A. Gilbert, Ernest Keefer-Jacques, Tanaya Jadhav, Daniel Antfolk, Qianqian Ming, Nicolette Valente, Grace Tzun-Wen Shaw, Christopher J. Sottolano, Grace Matwijec, Vincent C. Luca, Kathleen M. Loomes, Ramakrishnan Rajagopalan, Tristan J. Hayeck, and Nancy B. Spinner

*Abstract:* Pathogenic variants in the *JAG1* gene are a primary cause of the multi-system disorder Alagille syndrome. Although variant detection rates are high for this disease, there is uncertainty associated with the classification of missense variants that leads to reduced diagnostic yield. Consequently, up to 85% of reported JAG1 missense variants have uncertain or conflicting classifications. We generated a library of 2,832 JAG1 nucleotide variants within exons 1-7, a region with a high number of reported missense variants, and designed a high-throughput assay to measure JAG1 membrane expression, a requirement for normal function. After calibration using a set of 175 known or predicted pathogenic and benign variants included within the variant library, 486 variants were characterized as functionally abnormal (n=277 abnormal and n=209 likely abnormal), of which 439 (90.3%) were missense. We identified divergent membrane expression occurring at specific residues, indicating that loss of the WT residue itself does not drive pathogenicity, a finding supported by structural modeling data and with broad implications for clinical variant classification both for Alagille syndrome and globally across other disease genes. Of 144 uncertain variants reported in patients undergoing clinical or research testing, 27 had functionally abnormal membrane expression and inclusion of our data resulted in the reclassification of 26 to likely pathogenic. Functional evidence augments the classification of genomic variants, reducing uncertainty and improving diagnostics. Inclusion of this repository of functional evidence during JAG1 variant reclassification will significantly affect resolution of variant pathogenicity, making a critical impact on the molecular diagnosis of Alagille syndrome.

![Figure 1](https://github.com/tris-10/JAG1-Membrane-Expression-Assay/blob/main/figures/Graphical%20Abstract1_rev.tif)

**Figure 1: The impact on membrane expression of 2,832 JAG1 variants**

-----------------------------------------------------

## Running Variant Score Classification
Python 2.11.5 and libraries required: `argparse, datetime, numpy, pandas, labelsAndConstants`

Sample input files and examples coming.

To run variantScoreClassification.py:

```
usage: variantScoreClassification.py [-h] [-bm BIN1FILE] [-bw BIN2FILE] [-m KNOWNMUTFILE] [-w BENIGNFILE] [-o OUTPUTTABLE] [-w1 BINWEIGHT1] [-w2 BINWEIGHT2]
                                     [-mc MUTANTCALIBRATIONSET] [-wc WTCALIBRATIONSET] [--debug]                                                            
                                                                                                                                                            
options:                                                                                                                                                    
  -h, --help            show this help message and exit                                                                                                     
  -bm BIN1FILE, --bin1File BIN1FILE                                                                                                                         
                        Input bin file containing mutant allele counts                                                                                      
  -bw BIN2FILE, --bin2File BIN2FILE                                                                                                                         
                        Input bin file containing wild-type allele counts                                                                                   
  -m KNOWNMUTFILE, --knownMutFile KNOWNMUTFILE
                        Mutant variant control file
  -w BENIGNFILE, --benignFile BENIGNFILE
                        Wild-type variant control file
  -o OUTPUTTABLE, --outputTable OUTPUTTABLE
                        Output directory
  -w1 BINWEIGHT1, --binWeight1 BINWEIGHT1
                        Scalar weight applied mutant variant scores
  -w2 BINWEIGHT2, --binWeight2 BINWEIGHT2
                        Scalar weight applied wild-type variant scores
  -mc MUTANTCALIBRATIONSET, --mutantCalibrationSet MUTANTCALIBRATIONSET
                        Variants used for mutant calibration set (variant IDs separated by commas)
  -wc WTCALIBRATIONSET, --wtCalibrationSet WTCALIBRATIONSET
                        Variants used for wild-type calibration set (variant IDs separated by commas)
  --debug
```


### Running Classify Cross Replicates
Python 2.11.5, and libraries required: `argparse, datetime, numpy, pandas, scipy.stats, labelsAndConstants, plotter`

Sample input files and examples coming.

To run classifyCrossReplicates.py:

```
usage: classifyCrossReplicates.py [-h] [-m PATHTRUTHFILE] [-w BENIGNTRUTHFILE] [-o OUTPUTPREFIX] [-r REPLICATETABLEFILE] [-s NUMSORTS] [-n MINNUMREPLICATES] [-db DBBENIGN]
                                  [-dm DBPATH] [-pl LABELS] [-pc COLORS] [-b N_BIN] [-y YMAXOFFSET] [--fixedYAxis] [--includeLikely] [--plotAllHist] [--debug]

options:
  -h, --help            show this help message and exit
  -m PATHTRUTHFILE, --pathTruthFile PATHTRUTHFILE
                        Mutant variant control file
  -w BENIGNTRUTHFILE, --benignTruthFile BENIGNTRUTHFILE
                        Wild-type variant control file
  -o OUTPUTPREFIX, --outputPrefix OUTPUTPREFIX
  -dm DBPATH, --dbPath DBPATH
                        Mutant variants databases used in analysis (Default: [all variants])
  -pl LABELS, --labels LABELS
                        Plot legend labels (Separate by commas)
  -pc COLORS, --colors COLORS
                        Plot colors (Separate by commas)
  -b N_BIN, --N_bin N_BIN
                        Number of bins used for plotting histograms
  -y YMAXOFFSET, --yMaxOffset YMAXOFFSET
                        Scalar multiplier of the histogram y-axis
  --fixedYAxis          Sets Y-axis as fixed the maximum histogram value
  --includeLikely       Include likely pathogenic variants in OddsPath score
  --plotAllHist         Plot all variants in single histogram (requires only one color/labels)
  --debug
```