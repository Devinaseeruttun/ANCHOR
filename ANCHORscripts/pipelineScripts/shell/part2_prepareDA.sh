#!/bin/bash
set -ex


dir0="YourMainPath"
iniFile="${dir0}/metadata/pipe.ini"
iniFile2="${dir0}/metadata/pipe_part2.ini"
dirpipe="${dir0}/pipelineScripts"
runDir="${dir0}/run"
successDir="${runDir}/successfulRuns"
expName=$(grep "^expName" ${iniFile} | cut -d"=" -f2)
databaseList=$(grep "^databaseList" ${iniFile} | cut -d"=" -f2)
AnchorMinBlastIdentity=$(grep "^AnchorMinBlastIdentity" ${iniFile} | cut -d"=" -f2)
lowCountSeqThreshold=$(grep "^lowCountSeqThreshold" ${iniFile} | cut -d"=" -f2)
primerSelectionBypass=$(grep "^primerSelectionBypass" ${iniFile} | cut -d"=" -f2)
cutoff=$(grep "^cutoff" ${iniFile} | cut -d"=" -f2)
design=${dir0}/metadata/design.txt
minReplicates=$(grep "^minReplicates" ${iniFile2} | cut -d"=" -f2)
colMinReplicates=$(grep "^colMinReplicates" ${iniFile2} | cut -d"=" -f2)


echo -e "\n-----------------------------------\nPARAMETERS\n"
echo -e "Experiment name: ${expName}"
echo -e "Main directory: $dir0"
echo -e "Parameter file: ${iniFile}"
echo -e "Pipeline folder: $dirpipe"
echo -e "Database list: $databaseList"
echo -e "Primer selection: $primerSelectionBypass"
echo -e "Identity threshold for BLASTn: $AnchorMinBlastIdentity"
echo -e "Identity and coverage thresholds for BLASTn: $lowCountSeqThreshold"
echo -e "Anchor threshold: $cutoff"
echo -e "-----------------------------------\n\n"



if ! [[ -e "${design}" && -s "${design}" ]];
then 
	echo -e "\n-----\nERROR\n${design} does not exist or is empty.\nYou need a design file (name of samples (first column) and conditions(other columns)) in the following folder:${dir0}/metadata/design.txt" 
fi


if ! [[ -e "${iniFile2}" && -s "${iniFile2}" ]];
then 
	echo -e "\n-----\nERROR\n${iniFile2} does not exist or is empty. Copy the template (${dir0}/pipelineScripts/iniFileTemplates/pipe_part2.ini) and place it in the metadata folder ${dir0}/metadata "; 
fi

conditions=$(head -n1 ${design} | cut -f2- | sed "s/\t/\n/g" | grep -v -e "Samples" -e "Replicates" | sed "s/\n/ /")


mkdir -p ${dir0}/part2_prepareDA/phyloseq_output
cd ${dir0}/part2_prepareDA
ln -nsf ${dir0}/Results_${expName}_${cutoff}_${AnchorMinBlastIdentity}_${lowCountSeqThreshold}/Phyloseq phyloseq_input


echo -e "---\nR suite part1: SPARSITY AND RAREFACTION\n"
cp ${dirpipe}/R/part1_sparsity_and_rarefaction_template.R ${dir0}/part2_prepareDA/part1_sparsity_and_rarefaction.R
#A is is the folder with the phyloseq inputfile
#B is the output folder from phyloseq
#C is the OTU TABLE
#D is SAMPLE TABLE
#E is TAX TABLE
sed -i "s#AAAAA#${dir0}/part2_prepareDA/phyloseq_input#g" part1_sparsity_and_rarefaction.R
sed -i "s#BBBBB#${dir0}/part2_prepareDA/phyloseq_output#g" part1_sparsity_and_rarefaction.R
sed -i "s/CCCCC/otu_table.txt/g" part1_sparsity_and_rarefaction.R
sed -i "s/DDDDD/sample_data.txt/g" part1_sparsity_and_rarefaction.R
sed -i "s/EEEEE/taxonomy_table.txt/g" part1_sparsity_and_rarefaction.R
sed -i "s/FFFFF/${minReplicates}/" part1_sparsity_and_rarefaction.R
sed -i "s/GGGGG/${colMinReplicates}/" part1_sparsity_and_rarefaction.R
sed -i "1s/^/.libPaths\( c\( \.libPaths\(\), \"~\/R_custom_libraries\"\)\)\n/" part1_sparsity_and_rarefaction.R
echo -e "Job part1_sparsity_and_rarefaction"
Rscript part1_sparsity_and_rarefaction.R



touch ${successDir}/part2_prepareDA.ok
mv ${dir0}/part2_prepareDA ${runDir}

echo -e "\n---\npart10_create_phyloseq_files successfully finished.\n\nNow check the files:\n\t\tlibrarySizes_Samples.pdf\n\t\trarefaction_curves.pdf\n\t\tSparsity_test.pdf\n\nWhen you analyzed the results, use your own cutoff values in part11_run_phyloseq\n" && exit

echo -e "\n---\nWell done dude!" && exit 0
