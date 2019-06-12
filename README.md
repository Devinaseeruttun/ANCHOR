# ANCHOR

Anchor is a 16S rRNA gene amplicon pipeline for microbial analysis of multiple
environmental samples.

## Getting Started


### Installation with Docker

#### Docker

A pre-build docker image is available that includes all dependencies with the
exception of usearch, which has a licence that precludes binary distribution. We
recommend using our image and your own usearch binary to build a complete anchor
container.

First, register your email addresss at https://drive5.com/usearch/download.html
for USEARCH v9.2.64 and you will recieve an email with the download URL.

```sh
# Download the usearch binary
curl $URL_FROM_EMAIL > usearch9.2.64_i86linux32

# Make a tiny one-line Dockerfile
echo "FROM robsyme/anchor-onbuild" > Dockerfile.anchor-onbuild

# Build your Anchor Docker image.
docker build -t anchor -f Dockerfile.anchor-onbuild .
```

The `robsyme/anchor-onbuild` binary expects to find the file
`./usearch9.2.64_i86linux32` when it is used as a base image to build a new
container. It will grab the binary and fold it into your new docker image.

### Installation Locally

#### Dependencies

Anchor requires the following tools:
 - Mothur (used in assembling contigs. See:
   https://www.mothur.org/wiki/Installation)
 - BLAST (see: https://www.ncbi.nlm.nih.gov/books/NBK279671)
 - usearch9 (used for chimera detection. See:
   https://drive5.com/usearch/download.html)
 - python 2.7 should be already installed on the machine
   (https://docs.python-guide.org/starting/install/linux)

... and the following Python libraries:
- numpy
- pandas
- matplotlib
- seaborn
- openpyxl
- biopython

#### Anchor

1. Download (or clone) ANCHOR_v1.0 from github
   (https://github.com/gonzalezem/ANCHOR/tree/master/ANCHOR_v1.0)
2. If not already within your system, create a link (or copy) of mothur main
   file into `ANCHOR_v1.0/pipelineScripts/mothur`. ANCHOR will look for a file
   called simply mothur within `ANCHOR_v1.0/pipelineScripts/mothur/`.
3. Create a link (or copy) of usearch9 main file (usearch9) into
   `ANCHOR_v1.0/pipelineScripts/usearch9`. ANCHOR will look for a file called
   simply usearch9 within `ANCHOR_v1.0/pipelineScripts/usearch9/`
4. Build (or link) database(s) BLAST index into `ANCHOR_v1.0/db` folder. Note
   that NCBI 16S microbial database index is included in ANCHOR download. The
   name should be: databasename_index (ex: 16SMicrobial_index, nt_index,
   rdp_index, silva_index) (see how to build an index:
   https://www.ncbi.nlm.nih.gov/books/NBK279688)

### Running Anchor

#### Inputs

ANCHORS needs a few files and folders:
- A folder containing Illumina reads (ex. `PEread1_R1.fastq.gz`,
  `PEread1_R2.fastq.gz`, etc.)
- A design file containing at least 2 columns: Samples and Condition name.
  Example:

  | Samples       | myCondition   |
  | ------------- | ------------- |
  | PEread1       | Condition1    |
  | PEread2       | Condition1    |
  | PEread3       | Condition1    |
  | PEread4       | Condition2    |
  | PEread5       | Condition2    |
  | PEread6       | Condition2    |

Before running ANCHOR, prepare some room for it. The script
`preparation_script.sh` from within the ANCHOR folder will do this. This script
will check for dependencies and required files. It needs 3 arguments to be able
to run:
-  argument 1: raw read location (full path)
-  argument 2: folder from where ANCHOR will be run (full path)
-  argument 3: design file (full path)

Example:
```sh
cd mycomputer/myfolder/ANCHOR_v1.0
bash preparation_script.sh myIlluminaFiles/my_raw_reads mycomputer/myExperiment mycomputer/myfolder/myconditions.txt
```

The script can be run multiple times until there is no more error message.

#### Run

If running the script preparation_script.sh didn't retrun an error, you're good to go. The last output lines from preparation_script.sh run will tell you what to do (basically customizing ANCHOR to your needs and running the main script).

```sh
bash bashMe.sh
```

#### Output

When anchor is done a folder `Results_a_b_c_d` will be created (a-d values depend on user's input from `metadata/pipe.ini`).

A few folders are produced:
-  Summary (some summary files from ANCHOR run)
-  STAMP (inut for STAMP software)
-  Phyloseq (input for Phyloseq)
-  MicrobiomeAnalyst (input for microbiomeanalyst.ca)
-  metagenomeSeq (input for metagenomeSeq)
-  Excel (OTU table in excel format) and files
-  OTU and anchor sequences (fasta files)
-  OTU and anchor tables (txt files)


#### Test Run
Go inside `test_run` folder and run the following command (it takes around 2
minutes to run):

```sh
bash run_test.sh
```

### Figures From the Article
Figures [here](article).
