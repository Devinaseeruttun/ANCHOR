ANCHOR processes paired-end reads, assembles and annotates them.

What is needed:

Machine: ANCHOR currently runs on linux-like machines

dependency:
 - Mothur (used in assembling contigs. See: https://www.mothur.org/wiki/Installation)
 - BLAST (see: https://www.ncbi.nlm.nih.gov/books/NBK279671)
 - usearch9 (used for chimera detection. See: https://drive5.com/usearch/download.html)
 - R should be previously installed on the machine (see for example: https://www.r-bloggers.com/how-to-install-r-on-linux-ubuntu-16-04-xenial-xerus)
 - python 2.7 should be previously installed on the machine (https://docs.python-guide.org/starting/install/linux)
 - pandas library
 
 Downloads:
1. download (or clone) ANCHOR_v1.0 in github (https://github.com/gonzalezem/ANCHOR/tree/master/ANCHOR_v1.0)
2. If not already within your system, create a link (or copy) of mothur main file into ANCHOR_v1.0/pipelineScripts/mothur. ANCHOR will look for a file called simply mothur within ANCHOR_v1.0/pipelineScripts/mothur/.
3. Create a link (or copy) of usearch9 main file (usearch9) into ANCHOR_v1.0/pipelineScripts/usearch9 . ANCHOR will look for a file called simply usearch9 within ANCHOR_v1.0/pipelineScripts/usearch9/. 
4. Build (or link) database(s) BLAST index into ANCHOR_v1.0/db folder. Note that NCBI 16S microbial database index is included in ANCHOR download. The name should be: databasename_index (ex: 16SMicrobial_index, nt_index, rdp_index, silva_index) (see how to build an index: https://www.ncbi.nlm.nih.gov/books/NBK279688)
 

ANCHORS needs a few files and folders:
 	1. A folder containing Illumina reads (ex. PEread1_R1.fastq.gz, PEread1_R2.fastq.gz, etc.)
 	2. A design file containing at least 2 columns: Samples and any_Condition_Name. Example:
 Samples	myCondition
PEread1	Condition1
PEread2	Condition1
PEread3	Condition1
PEread4	Condition2
PEread4	Condition2
PEread6	Condition2

Before running ANCHOR, prepare room for it:
To do so, you just have to run the preparation_script.sh from within ANCHOR folder. This script will check for dependencies and files. It needs 3 arguments:
		- argument 1: raw read location (full path)
		- argument 2: folder from where ANCHOR will be run (full path)
		- argument 3: design file (full path)

Example:
cd mycomputer/myfolder/ANCHOR_v1.0
bash preparation_script.sh myIlluminaFiles/my_raw_reads mycomputer/myExperiment mycomputer/myfolder/myconditions.txt

You can run the script multiple times until there is no more error message. You will have errors if the depencies are not met


Running ANCHOR:
If preparation_script.sh didn't retrun an error, you're good to go. The last line of preparation_script.sh will tell you what to do (basically customizing ANCHOR to your needsand running the main script)


Output:
When anchor is done a folder Results_a_b_c_d will be created (a-d values depend on user's input from metadata/pipe.ini)
A few folders are produced:
	- Summary (some summary files from ANCHOR run)
	- STAMP (inut for STAMP software)
	- Phyloseq (input for Phyloseq)
	- MicrobiomeAnalyst (input for microbiomeanalyst.ca)
	- metagenomeSeq (input for metagenomeSeq)
	- Excel (OTU table in excel format)
and files:
	- OTU and anchor sequences (fasta files)
	- OTU and anchor tables (txt files)


