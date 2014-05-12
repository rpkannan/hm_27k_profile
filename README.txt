README

DATA ACCESS

Access data from TCGA Data Matrix: https://tcga-data.nci.nih.gov/tcga/dataAccessMatrix.htm
Filter Settings:
Select a disease: Select one at a time: BRCA,COAD,GBM,KIRC,LUAD,LUSC,OV,UCEC
Common to all cancers:
Data Type: DNA Methylation
Batch Number: All
Data Level: Level 3
Availability: Available
Preservation: Frozen
Center/Platform: JHU_USC (HumanMethylation27)

Select all cells and click “Build Archive”. Use the default file structure they provide (select all files to include in your archive).

DEFAULT FILE STRUCTURE

These scripts obey a strict file structure. All scripts and annotation files must be kept in the same directory. Data input (see DATA ACCESS) must obey TCGA file structure. Furthermore, data must be stored in a directory “Input” with a sub-directory for each cancer: “brca”,”coad”,etc. Generated data will be stored in a directory “Output” with a sub-directory for each cancer.

ORDER OF EXECUTION

Scripts should be run in the order presented below. 

FILE DESCRIPTIONS
Scripts:

from_source_hm27.py
Summary: primary script for extracting values from TCGA data matrix file downloads.
Functions: 
For each cancer, you must build a manifest that maps the sample barcode to its corresponding file name. 
Next, you can build a probe_dict. This is a dictionary with keys: barcodes and values: list of beta values associated with probes. Obeys ordering of probes_27.txt. This probe_dict can be discretized or contain raw values.
You can also build a list of all beta values across all samples that are located in CGIs and a list of all beta values across all samples that are located in non-CGI regions.

plots.py
Summary: using betas_raw_island.p and betas_raw_nonisland.p, generate histograms of these beta value distributions.
You can also ascertain the total number of probes that are hyper- and hypo-methylated.

to_cluster_in_R.py
Summary: used to write to file data stored in Python variables so that it can be manipulated in R. 
For each cancer:
Input: probe_dict_dis.p 
Output: probe_dict_dis.csv,hamming_samples.csv
The probe_dict_dis dictionary is converted to a Pandas dataframe. Any columns (ie. probes IDs) that are entirely NAs (TCGA does this when the probe coincides with an SNP) are removed. This dataframe is then written to a tab delimited file with rows indexed by barcodes and columns indexed by probe IDs.
This file also generates a full hamming distance matrix (distance calculated between samples) and writes it to a tab delimited file.

cluster_methods.py
Summary: this is an accessory script that holds functions used in clustering.
Functions included convert dictionaries to data matrices, calculate full and lower hamming distance matrices and write them to tab delimited files or pickle them, and methods for clustering and visualizing dendrograms in python.

cluster_by_cancer.R
Summary: this script uses the output from to_cluster_in_R.py to cluster samples by cancer, produce dendrograms, and determine members of cancer sub-types.
For each cancer:
Input: probe_dict_dis.p, hamming_samples.p
Output: by_sample_dendro.png,(for some cancers) cut1.csv…cut4.csv
Data matrix and distance matrix are used along with the hclust() function to produce a dendrogram. Use Ward’s clustering option:
hclust(distance_matrix,method=“ward”)

resolve.R
Summary: for each cancer, probe beta values are resolved to their corresponding gene beta values
For each cancer:
Input: probe_dict_dis.p
Output: gene_dict.p

sub_cancers.py
Summary: partition probe_dict_dis.p into sub-cancer dictionaries. Generate consensus sequences for these sub-cancers and calculate the hamming distances between these sub-cancers. Write consensus sequences and hamming distance matrix to file for use in R.

cluster_sub_cancers.R
Summary: Uses input from sub_cancers.py to cluster sub-cancers. Uses same methods and parameters as cluster_by_cancer.R.

hypermeth_44.py
Summary: given hypermeth_targets.txt, this script finds the probes these genes map to and builds a consensus sequence of this information by cancer.

Annotation Files:

adf_27.p
Object type: pickled Pandas dataframe that includes relevant annotation information obtained from TCGA: https://tcga-data.nci.nih.gov/tcga/tcgaPlatformDesign.jsp
Rows: probe IDs
Columns: [Chr, MapInfo, Next_Base, TSS_Coordinate, Gene_ID, Symbol, Synonym, Accession, GID, Annotation, Product, Distance_to_TSS, CPG_ISLAND, CPG_ISLAND_LOCATIONS]

probes_27.txt
Object type: list of all probe names in sorted order. This is the primary indexing tool for accessing data in probe_dict objects

island_status.txt
Object type: list indexed by probes_27.txt that has value “True” if corresponding probe is found in a CGI and “False” otherwise

