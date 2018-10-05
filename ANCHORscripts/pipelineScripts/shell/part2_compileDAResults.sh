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
design=${dir0}/metadata/design.txt
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

conditions=$(head -n1 ${design} | cut -f2- | sed "s/\t/\n/g" | grep -v -e "Samples" -e "Replicates" | sed "s/\n/ /")


############################################ SPARSITY ####################################################################3

mkdir -p ${dir0}/Results_${expName}_${cutoff}_${AnchorMinBlastIdentity}_${lowCountSeqThreshold}/DA_analysis/sparsityControl
cd ${dir0}/Results_${expName}_${cutoff}_${AnchorMinBlastIdentity}_${lowCountSeqThreshold}/DA_analysis/sparsityControl
cp ${runDir}/part2_runDA/phyloseq_output/part2_denoised_data/*.pdf ./
cp ${runDir}/part2_runDA/phyloseq_output/part2_denoised_data/*.png ./




############################################ DA ANALYSIS ####################################################################3

echo -e "DA Analysis"
for cond in ${conditions}
do
	echo -e "Condition: ${cond}"
	mkdir -p ${dir0}/Results_${expName}_${cutoff}_${AnchorMinBlastIdentity}_${lowCountSeqThreshold}/DA_analysis/${cond}/DA
	cd ${dir0}/Results_${expName}_${cutoff}_${AnchorMinBlastIdentity}_${lowCountSeqThreshold}/DA_analysis/${cond}/DA
	cp ${runDir}/part2_runDA/phyloseq_output/${cond}/part3_DESeq2/DeSeq2_*_Signif.txt ./
	cp ${runDir}/part2_runDA/phyloseq_output/${cond}/part3_DESeq2/DeSeq2_*_all.txt ./
	cp ${runDir}/part2_runDA/phyloseq_output/${cond}/DA_multiplots/DA_multiplot*.pdf ./
	cp ${runDir}/part2_runDA/phyloseq_output/${cond}/part3_DESeq2/rlog_transformed_otu_table_${cond}.txt ./
	for i in *all.txt; 
	do 
		#Adding OTU_table.txt columns to all.txt file
		ln -nsf ${dir0}/Results_${expName}_${cutoff}_${AnchorMinBlastIdentity}_${lowCountSeqThreshold}/OTU_table.txt
		# 1. remove Taxonomy from DESeq2 output
		sed -i "1s/Kingdom/Domain/" ${i}
		domainCol=$(head -n1 ${i} | sed "s/\t/\n/g" | grep -n "Domain" | cut -d":" -f1)
		colToRemove=$((${domainCol} - 1))
		cut -f-${colToRemove} ${i}  | awk 'FNR>1' > __all
		# 2. join with counts.txt
		join  -1 1 -2 1 -t $'\t' -a 1 <(sort -t $'\t' -k1,1 __all) <(sort -t $'\t' -k1,1 OTU_table.txt) > _temp1
		#creating headers
		head -n1 ${i} | cut -f-${colToRemove} | sed "1s/^taxon/OTU/" > _headers_all
		join  -1 1 -2 1 -t $'\t' _headers_all <(head -n1 OTU_table.txt) > __headers
		#construct file
		cat __headers _temp1 > ${i}
		rm -f __all _temp1 __headers _headers_all
		echo -e "file: ${i}"
		new=$(basename ${i} .txt)
		FILE="${new}.xlsx"
		if  [ ! -e "$FILE" ]
		then 
			python ${dirpipe}/python/tsvToxlsx.py -i ${i}
			rm -f ${i}
		fi
	done
	set +e
	checkDE=$(ls -1 *Signif.txt)
	set -e
	if [ -n "${checkDE}" ]; then
		echo -e "Adding OTU ambiguity to DA table"
		for i in *Signif.txt
		do 
			# 1. remove Taxonomy from DESeq2 output (from kingdom column to the end)
			sed -i "1s/Kingdom/Domain/" ${i}
			domainCol=$(head -n1 ${i} | sed "s/\t/\n/g" | grep -n "Domain" | cut -d":" -f1)
			colToRemove=$((${domainCol} - 1))
			cut -f-${colToRemove} ${i}  | awk 'FNR>1' > __all
			# 2. join with counts.txt
			join  -1 1 -2 1 -t $'\t' -a 1 <(sort -t $'\t' -k1,1 __all) <(sort -t $'\t' -k1,1 OTU_table.txt) > _temp1
			#creating headers
			head -n1 ${i} | cut -f-${colToRemove} | sed "1s/^/OTU/" > _headers_all
			join  -1 1 -2 1 -t $'\t' _headers_all <(head -n1 OTU_table.txt) > __headers
			#construct file
			cat __headers _temp1 > ${i}
			rm -f __all _temp1 __headers _headers_all
			echo -e "file: ${i}"
			new=$(basename ${i} .txt)
			FILE="${new}.xlsx"
			set +e
			if  [ ! -e "$FILE" ]
			then 
				python ${dirpipe}/python/tsvToxlsx.py -i ${i}
				rm -f ${i}
			fi
			set -e
		done
	else
		echo -e "\nNo DE in ${cond}\n"
	fi
	FILE="rlog_transformed_otu_table_${cond}.xlsx"
	if  [ ! -e "$FILE" ]
	then 
		python ${dirpipe}/python/tsvToxlsx.py -i rlog_transformed_otu_table_${cond}.txt
		rm -f rlog_transformed_otu_table_${cond}.txt
	fi
	rm -f OTU_table.txt
done



############################################ ORDINATION ####################################################################3

echo -e "Ordination"
for cond in ${conditions}
do
	echo -e "Condition: ${cond}"
	mkdir -p ${dir0}/Results_${expName}_${cutoff}_${AnchorMinBlastIdentity}_${lowCountSeqThreshold}/DA_analysis/${cond}/Ordination
	cd ${dir0}/Results_${expName}_${cutoff}_${AnchorMinBlastIdentity}_${lowCountSeqThreshold}/DA_analysis/${cond}/Ordination
	cp ${runDir}/part2_runDA/phyloseq_output/${cond}/part4_ordination/* ./
done


#Parse ordination statistical test tables
for cond in ${conditions}
do
	cd ${dir0}/Results_${expName}_${cutoff}_${AnchorMinBlastIdentity}_${lowCountSeqThreshold}/DA_analysis/${cond}/Ordination
	echo -e "Statistical test\tsignificance probability" > statistical_analysis_${cond}.txt
	rm -f Analysis_of_Similarity_rlog_*.txt
	for statTest in PERMANOVA_DESEQ2NormCounts*.txt
	do
		#I'll remove all lines including and before "adonis2"
		adonis2LineNumber=$(grep -n "adonis2" ${statTest} | cut -d":" -f1)
		awk -v var="${adonis2LineNumber}" 'FNR>var' ${statTest} > _temp
		sed -i "s/  */\t/g" _temp
		pvalCol=$(head -n1 _temp | sed "s/\t/\n/g" | grep -n "Pr(>F)" | cut -d":" -f1)
		probability=$(cut -f${pvalCol} _temp | perl -ne '(!/^\s+$/)&&print' | tail -n1)
		echo -e "Permanova\t${probability}" >> statistical_analysis_${cond}.txt
		rm -f _temp ${statTest}
	done
	for statTest in Anova_RDA_DESEQ2NormCounts*.txt
	do
		#I'll remove all lines including and before "adonis2"
		anovaLineNumber=$(grep -n "rda(formula =" ${statTest} | cut -d":" -f1)
		awk -v var="${anovaLineNumber}" 'FNR>var' ${statTest} > _temp
		sed -i "s/  */\t/g" _temp
		probability=$(cut -f${pvalCol} _temp | perl -ne '(!/^\s+$/)&&print' | tail -n1)
		echo -e "Anova_RDA\t${probability}" >> statistical_analysis_${cond}.txt
		rm -f _temp ${statTest}
	done
	for statTest in Anova_CAP_DESEQ2NormCounts*.txt
	do
		#I'll remove all lines including and before "adonis2"
		anovaLineNumber=$(grep -n "capscale(formula =" ${statTest} | cut -d":" -f1)
		awk -v var="${anovaLineNumber}" 'FNR>var' ${statTest} > _temp
		sed -i "s/  */\t/g" _temp
		probability=$(cut -f${pvalCol} _temp | perl -ne '(!/^\s+$/)&&print' | tail -n1)
		echo -e "Anova_CAP\t${probability}" >> statistical_analysis_${cond}.txt
		rm -f _temp ${statTest}
	done
	for statTest in Analysis_of_Similarity_DESEQ2NormCounts*.txt
	do
		probability=$(grep "Significance:" ${statTest} | sed "s/ $//" | rev | cut -d" " -f1 | rev)
		echo -e "ANOSIM\t${probability}" >> statistical_analysis_${cond}.txt
		rm -f  ${statTest}
	done
	#reorder statistical_analysis
	cat <(head -n1 statistical_analysis_${cond}.txt) <(awk 'FNR>1' statistical_analysis_${cond}.txt | sort -t $'\t' -k1,1) > _temp
	mv _temp statistical_analysis_${cond}.txt
	python ${dirpipe}/python/tsvToxlsx.py -i statistical_analysis_${cond}.txt
	rm -f statistical_analysis_${cond}.txt
done



############################################ RICHNESS ####################################################################3

echo -e "Richness"
for cond in ${conditions}
do
	echo -e "Condition: ${cond}"
	mkdir -p ${dir0}/Results_${expName}_${cutoff}_${AnchorMinBlastIdentity}_${lowCountSeqThreshold}/DA_analysis/${cond}/Richness
	cd ${dir0}/Results_${expName}_${cutoff}_${AnchorMinBlastIdentity}_${lowCountSeqThreshold}/DA_analysis/${cond}/Richness
	cp ${runDir}/part2_runDA/phyloseq_output/${cond}/part5_richness_plots/* ./
done

#Parse richness t-test tables
for cond in ${conditions}
do
	cd ${dir0}/Results_${expName}_${cutoff}_${AnchorMinBlastIdentity}_${lowCountSeqThreshold}/DA_analysis/${cond}/Richness
	for ttest in t_test_richness_*.txt
	do
		pvalCol=$(head -n1 ${ttest} | sed "s/\t/\n/g" | grep -n "p.value" | cut -d":" -f1)
		cut -f1-${pvalCol} ${ttest} > _temp
		sed -i "1s/^\t/alpha-diversity measures\t/" _temp
		sed -i "s/estimate\.mean in group/richness estimate in/g" _temp
		mv _temp ${ttest}
		python ${dirpipe}/python/tsvToxlsx.py -i ${ttest}
		rm -f ${ttest}
	done

done


############################################ KRONA ####################################################################3

echo -e "Relative_abundance"
for cond in ${conditions}
do
	echo -e "Condition: ${cond}"
	mkdir -p ${dir0}/Results_${expName}_${cutoff}_${AnchorMinBlastIdentity}_${lowCountSeqThreshold}/DA_analysis/${cond}/Krona
	cd ${dir0}/Results_${expName}_${cutoff}_${AnchorMinBlastIdentity}_${lowCountSeqThreshold}/DA_analysis/${cond}/Krona
	cp ${runDir}/part2_runDA/phyloseq_output/${cond}/part7_stackedBars/* ./
done




############################################ SUMMARY OF DA ANALYSIS ####################################################################3
echo -e "Summary ofdifferntial abundance analysis"
mkdir -p ${dir0}/Results_${expName}_${cutoff}_${AnchorMinBlastIdentity}_${lowCountSeqThreshold}/Summary
cd ${dir0}/Results_${expName}_${cutoff}_${AnchorMinBlastIdentity}_${lowCountSeqThreshold}/Summary

rm -f __*

for cond in ${conditions}
do
	cd ${runDir}/part2_runDA/phyloseq_output/${cond}/part3_DESeq2
	set +e
	checkDE=$(ls -1 *Signif.txt)
	set -e
	if [ -n "${checkDE}" ]; then
		echo -e "Adding OTU ambiguity to DA table"
		for i in *Signif.txt
		do
			regulationCol=$(head -n1 ${i} | sed "s/\t/\n/g" | grep -n "Final_Regulation" | cut -d":" -f1)
			cut -f1 ${i} | awk 'FNR>1' >> ${dir0}/Results_${expName}_${cutoff}_${AnchorMinBlastIdentity}_${lowCountSeqThreshold}/Summary/__allDA
			cut -f1,${regulationCol} ${i} | awk 'FNR>1' > ${dir0}/Results_${expName}_${cutoff}_${AnchorMinBlastIdentity}_${lowCountSeqThreshold}/Summary/__${cond}_DA_list
		done
	fi
done

cd ${dir0}/Results_${expName}_${cutoff}_${AnchorMinBlastIdentity}_${lowCountSeqThreshold}/Summary

#unique OTU
sort __allDA | uniq > _temp
mv _temp __allDA


echo -ne "OTU" > summaryDA.txt
for cond in ${conditions}
do
	echo -ne "\t${cond}" >> summaryDA.txt
done
echo -e "" >> summaryDA.txt




while read otu
do
	echo -ne "${otu}">> summaryDA.txt
	for cond in ${conditions}
	do
		regulation=$(cat __${cond}_DA_list | grep "^${otu}" | cut -f2)
		echo -ne "\t${regulation}" >> summaryDA.txt
	done
	echo -e "" >> summaryDA.txt
done<__allDA

#clean table
sed -i "s/\t\t/\t-\t/g" summaryDA.txt
sed -i "s/\t$/\t-/" summaryDA.txt
python ${dirpipe}/python/tsvToxlsx.py -i summaryDA.txt

rm -f __* summaryDA.txt



touch ${successDir}/part2_compileDAResults.ok

echo -e "\n---\npart2_compileDAResults.sh ran smoothly.\nWell done dude!" && exit 0












