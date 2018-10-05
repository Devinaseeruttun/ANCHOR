#!/bin/bash
set -e


: <<'END'
What you need:
1. anchor_table.txt from chimeraFlag.sh

what it does:
1. creates the OTU count + taxonomy file


END




dir0="YourMainPath"
iniFile="${dir0}/metadata/pipe.ini"
dirpipe="${dir0}/pipelineScripts"
runDir="${dir0}/run"
successDir="${runDir}/successfulRuns"
expName=$(grep "^expName" ${iniFile} | cut -d"=" -f2)
cutoff=$(grep "^cutoff" ${iniFile} | cut -d"=" -f2)
AnchorMinBlastIdentity=$(grep "^AnchorMinBlastIdentity" ${iniFile} | cut -d"=" -f2)
lowCountSeqThreshold=$(grep "^lowCountSeqThreshold" ${iniFile} | cut -d"=" -f2)
primerSelectionBypass=$(grep "^primerSelectionBypass" ${iniFile} | cut -d"=" -f2)


echo -e "\n-----------------------------------\nPARAMETERS\n"
echo -e "Experiment name: ${expName}"
echo -e "Main directory: $dir0"
echo -e "Parameter file: ${iniFile}"
echo -e "Pipeline folder: $dirpipe"
echo -e "Database list: $databaseList"
echo -e "Design file: $Conditions"
echo -e "Anchor count threshold: $cutoff"
echo -e "Anchor identity threshold for BLASTn: $AnchorMinBlastIdentity"
echo -e "Primer selection: $primerSelectionBypass"
echo -e "Identity and coverage thresholds for BLASTn: $lowCountSeqThreshold"
echo -e "-----------------------------------\n\n"



mkdir -p ${dir0}/Results_${expName}_${cutoff}_${AnchorMinBlastIdentity}_${lowCountSeqThreshold}/Excel
cd ${dir0}/Results_${expName}_${cutoff}_${AnchorMinBlastIdentity}_${lowCountSeqThreshold}
cp ${runDir}/OTUandCounts/OTU_table.txt ./
#transform files into Excel-ready files
python ${dirpipe}/python/tsvToxlsx.py -i OTU_table.txt
mv  OTU_table.xlsx ./Excel
cp -r ${runDir}/Summary/ ./
cp ${runDir}/chimeraFlag/anchor_table.txt ./
python ${dirpipe}/python/tsvToxlsx.py -i anchor_table.txt
mv  anchor_table.xlsx ./Excel
#create a fasta file from anchor_table.txt
sequenceColNumber=$(head -n1 anchor_table.txt | sed "s/\t/\n/g" | grep -n "sequence" | cut -d":" -f1)
cut -f1,${sequenceColNumber} anchor_table.txt | awk 'FNR>1' | sed "s/^/>/" | sed "s/\t/\n/" > anchors.fasta
#create a fasta file from OTU_table.txt
sequenceColNumber=$(head -n1 OTU_table.txt | sed "s/\t/\n/g" | grep -n "sequence" | cut -d":" -f1)
cut -f1,${sequenceColNumber} OTU_table.txt | awk 'FNR>1' | sed "s/^/>/" | sed "s/\t/\n/" > OTU.fasta


echo -e "Create microbiomeanalyst.ca-ready files"
mkdir -p ${dir0}/Results_${expName}_${cutoff}_${AnchorMinBlastIdentity}_${lowCountSeqThreshold}/MicrobiomeAnalyst
cd ${dir0}/Results_${expName}_${cutoff}_${AnchorMinBlastIdentity}_${lowCountSeqThreshold}/MicrobiomeAnalyst
ln -nsf ../OTU_table.txt anchorOTU_table
#OTU TABLE
#counts start right after sequence column and ends before totalcounts column
sequenceCol=$(head -n1 anchorOTU_table | sed "s/\t/\n/g" | grep -n "sequence" | cut -d":" -f1)
totalcountsCol=$(head -n1 anchorOTU_table | sed "s/\t/\n/g" | grep -n "totalcounts" | cut -d":" -f1)
startCountCols=$((${sequenceCol} + 1))
endCountCols=$((${totalcountsCol} - 1))
cut -f1,${startCountCols}-${endCountCols} anchorOTU_table | sed "1s/^OTU/#NAME/" > OTU_table.txt
#TAXONOMY TABLE
#Taxonomy starts at Domain column and ends at Species column
domainCol=$(head -n1 anchorOTU_table | sed "s/\t/\n/g" | grep -n "Domain" | cut -d":" -f1)
speciesCol=$(head -n1 anchorOTU_table | sed "s/\t/\n/g" | grep -n "Species" | cut -d":" -f1)
cut -f1,${domainCol}-${speciesCol} anchorOTU_table | sed "1s/^OTU/#TAXONOMY/" > Taxonomy_table.txt
#METADATA FILE
sed "1s/^[^\t]*\t/#NAME\t/" ${dir0}/metadata/design.txt > Metadata_file.txt
rm -f anchorOTU_table



echo -e "Create phyloseq-ready files"
ln -nsf ../OTU_table.txt anchorOTU_table
mkdir -p ${dir0}/Results_${expName}_${cutoff}_${AnchorMinBlastIdentity}_${lowCountSeqThreshold}/Phyloseq
cd ${dir0}/Results_${expName}_${cutoff}_${AnchorMinBlastIdentity}_${lowCountSeqThreshold}/Phyloseq
ln -nsf ../OTU_table.txt anchorOTU_table
#Sample Data
sed "1s/.[^\t]*\t/\t/" ${dir0}/metadata/design.txt > sample_data.txt
#Taxonomy table
cut -f1,${domainCol}-${speciesCol} anchorOTU_table | sed "1s/^OTU\t/\t/" > taxonomy_table.txt
#Dada2 expect Kingdom instead of Domain:
sed -i "1s/Domain/Kingdom/" taxonomy_table.txt
sed -i "1s/^OTU//" taxonomy_table.txt


echo -e "\n\nCreating OTU TABLE"
echo -e "Species"
#OTU table
cut -f1,${startCountCols}-${endCountCols} anchorOTU_table | sed "1s/^OTU\t/\t/" > otu_table.txt
sed -i "1s/^OTU//" otu_table.txt
rm -f anchorOTU_table



#CLEAN WORKDIR
mkdir -p ${successDir}
touch ${successDir}/compileResults.ok


echo -e "\n---\ncompileResults.sh is exiting normally.\nWell done!"

