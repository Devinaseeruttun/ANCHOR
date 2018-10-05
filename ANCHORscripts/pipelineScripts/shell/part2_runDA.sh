#!/bin/bash
set -e





dir0="/media/emmanuel/storage4/ISS_testrun"
iniFile="${dir0}/metadata/pipe.ini"
iniFile2="${dir0}/metadata/pipe_part2.ini"
dirpipe="${dir0}/pipelineScripts"
runDir="${dir0}/run"
successDir="${runDir}/successfulRuns"
design="${dir0}/metadata/design.txt"
expName=$(grep "^expName" ${iniFile} | cut -d"=" -f2)
databaseList=$(grep "^databaseList" ${iniFile} | cut -d"=" -f2)
AnchorMinBlastIdentity=$(grep "^AnchorMinBlastIdentity" ${iniFile} | cut -d"=" -f2)
lowCountSeqThreshold=$(grep "^lowCountSeqThreshold" ${iniFile} | cut -d"=" -f2)
primerSelectionBypass=$(grep "^primerSelectionBypass" ${iniFile} | cut -d"=" -f2)
cutoff=$(grep "^cutoff" ${iniFile} | cut -d"=" -f2)

rankingFile=$(grep "^rankingFile" ${iniFile2} | cut -d"=" -f2)
phylumRankingFile=$(grep "^phylumRankingFile" ${iniFile2} | cut -d"=" -f2)
minReplicates=$(grep "^minReplicates" ${iniFile2} | cut -d"=" -f2)
colMinReplicates=$(grep "^colMinReplicates" ${iniFile2} | cut -d"=" -f2)
countPerSampleCutoff=$(grep "^countPerSampleCutoff" ${iniFile2} | cut -d"=" -f2)
countPerTaxonCutoff=$(grep "^countPerTaxonCutoff" ${iniFile2} | cut -d"=" -f2)
sampleListToRemove=$(grep "^sampleListToRemove" ${iniFile2} | cut -d"=" -f2)
dataTransformation=$(grep "^dataTransformation" ${iniFile2} | cut -d"=" -f2)
effectSize=$(grep "^effectSize" ${iniFile2} | cut -d"=" -f2)
sparsityControl=$(grep "^sparsityControl" ${iniFile2} | cut -d"=" -f2)
FDRThreshold=$(grep "^FDRThreshold" ${iniFile2} | cut -d"=" -f2)
minimumOTUOcurrence=$(grep "^minimumOTUOcurrence" ${iniFile2} | cut -d"=" -f2)
minimumNumberOfSampleOcurrence=$(grep "^minimumNumberOfSampleOcurrence" ${iniFile2} | cut -d"=" -f2)


if [ "${rankingFile}"  == "" ]; then
    rankingFile="None"
fi
if [ "${phylumRankingFile}"  == "" ]; then
    phylumRankingFile="None"
fi


if [ "${countPerSampleCutoff}"  == "" ]; then
    countPerSampleCutoff=1000
fi

if [ "${countPerTaxonCutoff}"  == "" ]; then
    countPerTaxonCutoff=$((${minReplicates}*3))
fi

if [ "${taxaListToRemove}"  == "" ]; then
    taxaListToRemove="None"
fi

if [ "${sampleListToRemove}"  == "" ]; then
    sampleListToRemove="None"
fi

if [ "${dataTransformation}"  == "" ]; then
    dataTransformation="RLOGBLIND"
fi

if [ "${effectSize}"  == "" ]; then
    effectSize=1
fi

if [ "${sparsityControl}"  == "" ]; then
    sparsityControl=0.9
fi

if [ "${FDRThreshold}"  == "" ]; then
    FDRThreshold=0.05
fi

if [ "${minimumOTUOcurrence}"  == "" ]; then
    minimumOTUOcurrence="None"
    minimumNumberOfSampleOcurrence="None"
fi


echo -e "\n-----------------------------------\nPARAMETERS\n"
echo -e "dir0: $dir0"
echo -e "Pipeline scripts folder : $dirpipe"
echo -e "Anchor BLASTn threshold: $highCountersThreshold"
echo -e "Low count sequences BLASTn threshold: $lowCountersThreshold"
echo -e "Minimum number of replicates: $minReplicates"
echo -e "Column name with minimum number of replicates: $colMinReplicates"
echo -e "Conditions: $design"
echo -e "Minimum OTU counts per sample: $countPerSampleCutoff"
echo -e "Minimum count of a given taxon across samples: $countPerTaxonCutoff"
echo -e "List of taxa to remove: $taxaListToRemove"
echo -e "List of samples to remove: $sampleListToRemove"
echo -e "Data transformation : $dataTransformation"
echo -e "Effect size: fold change of ${effectSize}"
echo -e "Sparsity control parameter  : $sparsityControl"
echo -e "Anchor cutoff: $anchorCutoff"
echo -e "DA plot ranking file: ${rankingFile}"
echo -e "DA plot phylum ranking file: ${phylumRankingFile}"
if [ "${minimumOTUOcurrence}"  != "None" ]; then
	perc=$(echo "scale=3; ${minimumNumberOfSampleOcurrence}*100" | bc)
	echo -e "Variance filter ON: remove taxa not seen more than ${minimumOTUOcurrence} times in at least ${perc}% of ${minReplicates} samples (i.e. the lowest number of replicated samples"
fi
echo -e "-----------------------------------\n\n"


conditions=$(head -n1 ${design} | cut -f2- | sed "s/\t/\n/g" | grep -v -e "Samples" -e "Replicates" | sed "s/\n/ /")


mkdir -p ${dir0}/part2_runDA/phyloseq_output
cd ${dir0}/part2_runDA
ln -nsf ${dir0}/Results_${expName}_${cutoff}_${AnchorMinBlastIdentity}_${lowCountSeqThreshold}/Phyloseq phyloseq_input

if [ "${rankingFile}"  != "None" ]; then
		echo -e "Ranking file for DA plot exists."
		cd ${dir0}/part2_runDA/phyloseq_input
		ln -nsf ${rankingFile} ranking.txt
		cd ${dir0}/part2_runDA
fi

if [ "${phylumRankingFile}"  != "None" ]; then
		echo -e "Phylum ranking file for DA plot exists."
		cd ${dir0}/part2_runDA/phyloseq_input
		ln -nsf ${phylumRankingFile} phylum_ranking.txt
		cd ${dir0}/part2_runDA
fi

FILE="${dir0}/part2_runDA/phyloseq_output/part2_denoised_data_species.ok"
if [ -e "${FILE}" ];then
    echo -e "\n-----\nSkipping DE NOISING DATA step!"
else
	echo -e "---\nR suite part2: DE NOISING DATA\n"
	cp ${dirpipe}/R/part2_denoised_data.R ${dir0}/part2_runDA/part2_denoised_data_species.R
	#A is is the folder with the phyloseq inputfile
	#B is the output folder from phyloseq
	#C is the OTU TABLE
	#D is SAMPLE TABLE
	#E is TAX TABLE
	sed -i "s#AAAAA#${dir0}/part2_runDA/phyloseq_input#g"  part2_denoised_data_species.R
	mkdir -p ${dir0}/part2_runDA/phyloseq_output/part2_denoised_data
	sed -i "s#BBBBB#${dir0}/part2_runDA/phyloseq_output/part2_denoised_data#g" part2_denoised_data_species.R
	sed -i "s/CCCCC/otu_table.txt/g" part2_denoised_data_species.R
	sed -i "s/DDDDD/sample_data.txt/g" part2_denoised_data_species.R
	sed -i "s/EEEEE/taxonomy_table.txt/g" part2_denoised_data_species.R
	sed -i "s/FFFFF/${minReplicates}/" part2_denoised_data_species.R ### This is the condition of reference for sparsity graph (use the condition you want, it needs replicates or it is useless)
	sed -i "s/IIIII/${colMinReplicates}/" part2_denoised_data_species.R
	
	
	#REMOVE SPARSITY AND LOW SAMPLE COUNT
	
	
	#remove taxa. Provide a list of taxa to remove (from Potential_noise_OTU_from_sparsity_test.txt)
	
	if [ "${taxaListToRemove}"  != "None" ]; then
		echo -e "Removing samples in ${taxaListToRemove}."
		sed -i "s #tax2Remove_path..* tax2Remove_path=\"${taxaListToRemove}\" g" part2_denoised_data_species.R
		sed -i "s/#TAXA2REMOVE//" part2_denoised_data_species.R
	
	else
		echo -e "Lookig at sparsity and remove sparse OTU if they exist" >&2
		checkSparseOTU=$(cat ${runDir}/part2_prepareDA/phyloseq_output/Potential_noise_OTU_from_sparsity_test.txt | wc -l)
		if [ ${checkSparseOTU} -gt 1 ]; then
			ln -nsf ${runDir}/part2_prepareDA/phyloseq_output/Potential_noise_OTU_from_sparsity_test.txt
			#Filter the taxas that have counts in less than 3 samples (absolute minimum number of replicates)
			awk 'FNR>1'  Potential_noise_OTU_from_sparsity_test.txt | awk -F'\t' '$6<3' | cut -f1 > _liste1
			#Filter the taxas which counts have a sparsity ratio higher than sparsityControl (default is 90%) 
			awk 'FNR>1'  Potential_noise_OTU_from_sparsity_test.txt | awk -F'\t' '$6>=3' |  awk -F'\t' -v var="${sparsityControl}" '($5>=var)' | cut -f1 > _liste2
			cat _liste1 _liste2 | sort | uniq > ${dir0}/part2_runDA/phyloseq_output/part11_tax2remove.txt
			rm -f _liste* Potential_noise_OTU_from_sparsity_test.txt
			sed -i "s #tax2Remove_path..* tax2Remove_path=\"${dir0}/part2_runDA/phyloseq_output/part11_tax2remove.txt\" g" part2_denoised_data_species.R
			sed -i "s/#TAXA2REMOVE//" part2_denoised_data_species.R
		fi
	fi
	
	if [ "${sampleListToRemove}"  != "None" ]; then
		sed -i "s #samples2Remove_path..* samples2Remove_path=\"${sampleListToRemove}\" g" part2_denoised_data_species.R
		sed -i "s/#SAMPLES2REMOVE//" part2_denoised_data_species.R
	fi
	
	#VARIANCE FILTER
	if [ "${minimumOTUOcurrence}"  != "None" ]; then
		sed -i "s/^#varianceFilter//" part2_denoised_data_species.R
		sed -i "s/JJJJJ/${minimumOTUOcurrence}/" part2_denoised_data_species.R
		sed -i "s/KKKKK/${minimumNumberOfSampleOcurrence}/" part2_denoised_data_species.R

	fi


	#countPerSampleCutoff : If needed, replace 1 by a value determined by you in the following line (or leave 1)
	sed -i "s/GGGGG/${countPerSampleCutoff}/" part2_denoised_data_species.R
	
	#countPerTaxonCutoff : If needed, replace 1 by a value determined by you in the following line (or leave 1)
	sed -i "s/HHHHH/${countPerTaxonCutoff}/" part2_denoised_data_species.R
	
	echo -e "Job part2_denoised_data_species"
	Rscript part2_denoised_data_species.R
	touch ${dir0}/part2_runDA/phyloseq_output/part2_denoised_data_species.ok
fi




FILE="${dir0}/part2_runDA/phyloseq_output/part3_DEseq.ok"
if [ -e "${FILE}" ];then
    echo -e "\n-----\nSkipping DA ANALYSIS step!"
else
	echo -e "---\nR suite part3: DESEQ2\n"
	cp ${dirpipe}/R/part3_DEseq_template.R ${dir0}/part2_runDA/part3_DEseq.R
	cd ${dir0}/part2_runDA
	sed -i "s#AAAAA#${dir0}/part2_runDA/phyloseq_input#g" part3_DEseq.R
	sed -i "s#BBBBB#${dir0}/part2_runDA/phyloseq_output#g" part3_DEseq.R
	sed -i "s/EEEEE/taxonomy_table.txt/g" part3_DEseq.R
	#Data tranformation
	sed -i "s/#${dataTransformation}//g" part3_DEseq.R
	#FDR cutoff
	sed -i "s/CCCCC/${FDRThreshold}/" part3_DEseq.R
	#remove the logfoldchnge cutoff
	sed -i "s/DDDDD/${effectSize}/" part3_DEseq.R
	#Remove cook's cutoof because I control the outliers
	#sed -i "s/diagdds, contrast/diagdds,cooksCutoff=FALSE, contrast/g" part3_DEseq.R
	
	for cond in ${conditions}
	do
		echo -e "Condition: ${cond}"
		sed "s/XXXXX/${cond}/g" part3_DEseq.R > part3_DEseq_${cond}.R 
		Rscript part3_DEseq_${cond}.R
	
	done
	touch ${dir0}/part2_runDA/phyloseq_output/part3_DEseq.ok
fi


FILE="${dir0}/part2_runDA/phyloseq_output/part4_ordination_plots.ok"
if [ -e "${FILE}" ];then
    echo -e "\n-----\nSkipping ORDINATION step!"
else
	echo -e "---\nR suite part4: ORDINATION\n"
	cp ${dirpipe}/R/part4_ordination_plots_template.R ${dir0}/part2_runDA/part4_ordination_plots.R
	cd ${dir0}/part2_runDA
	sed -i "s#AAAAA#${dir0}/part2_runDA/phyloseq_input#g" part4_ordination_plots.R
	sed -i "s#BBBBB#${dir0}/part2_runDA/phyloseq_output#g" part4_ordination_plots.R
	sed -i "s/DDDDD/sample_data.txt/g" part4_ordination_plots.R
	sed -i "s/EEEEE/taxonomy_table.txt/g" part4_ordination_plots.R
	#Data tranformation
	sed -i "s/#${dataTransformation}//g" part4_ordination_plots.R
	
	for cond in ${conditions}
	do
		echo -e "Condition: ${cond}"
		sed "s/XXXXX/${cond}/g" part4_ordination_plots.R > part4_ordination_plots_${cond}.R
		Rscript part4_ordination_plots_${cond}.R
	done
	touch ${dir0}/part2_runDA/phyloseq_output/part4_ordination_plots.ok
fi



FILE="${dir0}/part2_runDA/phyloseq_output/part5_richness_plots.ok"
if [ -e "${FILE}" ];then
    echo -e "\n-----\nSkipping RICHNESS step!"
else
	echo -e "---\nR suite part5: RICHNESS\n"
	cp ${dirpipe}/R/part5_richness_plots_template.R ${dir0}/part2_runDA/part5_richness_plots.R
	cd ${dir0}/part2_runDA
	sed -i "s#AAAAA#${dir0}/part2_runDA/phyloseq_input#g" part5_richness_plots.R
	sed -i "s#BBBBB#${dir0}/part2_runDA/phyloseq_output#g" part5_richness_plots.R
	sed -i "s/CCCCC/otu_table.txt/g" part5_richness_plots.R
	sed -i "s/DDDDD/sample_data.txt/g" part5_richness_plots.R
	sed -i "s/EEEEE/taxonomy_table.txt/g" part5_richness_plots.R
	
	for cond in ${conditions}
	do
		echo -e "Condition: ${cond}"
		sed "s/XXXXX/${cond}/g" part5_richness_plots.R > part5_richness_plots_${cond}.R 
		Rscript part5_richness_plots_${cond}.R
	done
	touch ${dir0}/part2_runDA/phyloseq_output/part5_richness_plots.ok
fi






FILE="${dir0}/part2_runDA/phyloseq_output/part6_heatmaps.ok"
if [ -e "${FILE}" ];then
    echo -e "\n-----\nSkipping HEATMAPS step!"
else
	echo -e "---\nR suite part6: HEATMAPS\n"
	cp ${dirpipe}/R/part6_heatmaps_template.R ${dir0}/part2_runDA/part6_heatmaps.R
	cd ${dir0}/part2_runDA
	sed -i "s#AAAAA#${dir0}/part2_runDA/phyloseq_input#g" part6_heatmaps.R
	sed -i "s#BBBBB#${dir0}/part2_runDA/phyloseq_output#g" part6_heatmaps.R
	
	sed -i "s/DDDDD/sample_data.txt/g" part6_heatmaps.R
	sed -i "s/EEEEE/taxonomy_table.txt/g" part6_heatmaps.R
	#Data tranformation
	sed -i "s/#${dataTransformation}//g" part6_heatmaps.R
	for cond in ${conditions}
	do
		echo -e "\n\n\n\n\nCondition: ${cond}"
		sed "s/XXXXX/${cond}/g" part6_heatmaps.R > part6_heatmaps_${cond}.R
		Rscript part6_heatmaps_${cond}.R
	done
	touch ${dir0}/part2_runDA/phyloseq_output/part6_heatmaps.ok
fi




FILE="${dir0}/part2_runDA/phyloseq_output/part7_stackedBars.ok"
if [ -e "${FILE}" ];then
    echo -e "\n-----\nSkipping STACKED BARS step!"
else
	echo -e "---\nR suite part7: STACKED BARS\n"
	cp ${dirpipe}/R/part7_stackedBars_template.R ${dir0}/part2_runDA/part7_stackedBars.R
	cd ${dir0}/part2_runDA
	sed -i "s#AAAAA#${dir0}/part2_runDA/phyloseq_input#g" part7_stackedBars.R
	sed -i "s#BBBBB#${dir0}/part2_runDA/phyloseq_output#g" part7_stackedBars.R
	
	sed -i "s/DDDDD/sample_data.txt/g" part7_stackedBars.R
	sed -i "s/EEEEE/taxonomy_table.txt/g" part7_stackedBars.R
	#adjust character  size
	sed -i "s/base_size=12/base_size=90/" part7_stackedBars.R
	sed -i "s/height=20/height=90/" part7_stackedBars.R
	#Data tranformation
	sed -i "s/#${dataTransformation}//g" part7_stackedBars.R
	
	for cond in ${conditions}
	do
		echo -e "Condition: ${cond}"
		sed "s/XXXXX/${cond}/g" part7_stackedBars.R > part7_stackedBars_${cond}.R 
		Rscript part7_stackedBars_${cond}.R
	
	done
	touch ${dir0}/part2_runDA/phyloseq_output/part7_stackedBars.ok
fi



FILE="${dir0}/part2_runDA/phyloseq_output/DA_multiplots.ok"
if [ -e "${FILE}" ];then
    echo -e "\n-----\nSkipping DA MULTIPLOT step!"
else
	echo -e "---\nR suite part8: DA MULTIPLOT\n"
	cp ${dirpipe}/R/DA_multiplots.R ${dir0}/part2_runDA/DA_multiplots.R
	cd ${dir0}/part2_runDA
	sed -i "s#AAAAA#${dir0}/part2_runDA/phyloseq_input#g" DA_multiplots.R
	sed -i "s#BBBBB#${dir0}/part2_runDA/phyloseq_output#g" DA_multiplots.R
	sed -i "s/EEEEE/taxonomy_table.txt/g" DA_multiplots.R
	#Data tranformation
	sed -i "s/#${dataTransformation}//g" DA_multiplots.R
	#FDR cutoff
	sed -i "s/CCCCC/${FDRThreshold}/" DA_multiplots.R
	#Effect size cutoff
	sed -i "s/DDDDD/${effectSize}/" DA_multiplots.R
	#Remove cook's cutoof because I control the outliers
	#sed -i "s/diagdds, contrast/diagdds,cooksCutoff=FALSE, contrast/g" DA_multiplots.R
	
	for cond in ${conditions}
	do
		echo -e "Condition : ${cond}"
		sed "s/XXXXX/${cond}/g" DA_multiplots.R > DA_multiplots_${cond}.R
		Rscript DA_multiplots_${cond}.R
	done
	touch ${dir0}/part2_runDA/phyloseq_output/DA_multiplots.ok
fi



touch ${successDir}/part2_runDA.ok
rm -f ${dir0}/part2_runDA/phyloseq_output/*.ok
mv ${dir0}/part2_runDA ${runDir}


echo -e "\n---\npart11_run_phyloseq ran smoothly.\nWell done dude!" && exit 0
