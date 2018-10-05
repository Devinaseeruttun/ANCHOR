#!/bin/bash
set -e


: <<'END'
What I do:
	The goal is to keep sequences that contain the exact primer in the beginning of the read sequence (amibiguous nucleotide are dealt with) and to remove added bases in case of staggered primers

1. calculate the length of all primers
2. sort them from smallest to largest length
3. Extract the first one as a full primer (assuming the staggrered primers will have one or more ncleotide added to it)
4. Put it in a file (staggeredGroup_${i}.txt)
5. Extract the staggered primers from that sequence and put them into the same file
6. Create a new primer file (all primers - the sequences in staggeredGroup_${i}.txt)
7. do 1 to 6 until the new primer file is empty
8. For each of the staggeredGroup_${i}.txt file:
	a. sort smallest primer length to largest
	b. extract the smallest value
	c. calculate how many bases are added to the smallest primer (that will command how many bases to trim on the next script)
9. Unzip fastq if needed
10. For each read and each group of primer, a script is called (extract_reads_with_primers.sh) that will keep the reads that have the primer(s) at the very beginning of their sequences
11. Clean PE reads (fastq_PE_fixer.sh) by removing any single reads created by the previous step


END



dir0="YourMainPath"
dirpipe="${dir0}/pipelineScripts/shell"
dirMeta="${dir0}/metadata"
dirreads="${dir0}/raw_reads"
successDir="${dir0}/run/successfulRuns"
primers="${dirMeta}/primers.txt"
primerSelectionBypass=$(grep "^primerSelectionBypass" $1 | cut -d"=" -f2)


echo -e "\n-----------------------------------\nPARAMETERS\n"
echo -e "Main directory: $dir0"
echo -e "Pipeline scripts directory: $dirpipe"
echo -e "Metadata directory: $dirMeta"
echo -e "Raw reads directory: ${dirreads}"
echo -e "Bypass primer selection: $primerSelectionBypass"
echo -e "-----------------------------------\n\n"

mkdir -p ${dir0}/fastqsWithPrimers
cd ${dir0}/fastqsWithPrimers

echo -e "Checking input before running the script"
FILE="${dirpipe}/extract_reads_with_primers.sh"
if [ ! -e "$FILE" ];
then
   echo -e "Shell script $FILE does not exist. You need it to run. Find it and place it here: ${dirpipe}/" >&2
   exit 1
fi


FILE="${dirpipe}/fastq_PE_fixer.sh"
if [ ! -e "$FILE" ];
then
   echo -e "Shell script $FILE does not exist. You need it to run. Find it and place it here: ${dirpipe}/" >&2
   exit 1
fi


if [ ! -d "$dirreads" ]; then
  echo -e "Couldn't find the reads directory (${dirreads}). Did you run part1_remove_phix.pbs.sh before running this script?"
fi

if [ "${primerSelectionBypass}"  != "YES" ]; then
	FILE="${primers}"
	if [ -e "$FILE" ];
	then
		# 1. Get forward and reverse column numbers 2. for each lin in the primer file (there could be several primers) get sequence information and create a file ${dirMeta}/primers_parsed.txt 
		forwardColNumber=$(head -n1 ${primers} | sed "s/\t/\n/g" | grep -n "forward" | cut -d":" -f1)
		reverseColNumber=$(head -n1 ${primers} | sed "s/\t/\n/g" | grep -n "reverse" | cut -d":" -f1)
		if ([ -n "${forwardColNumber}" ] && [ -n "${reverseColNumber}" ]); then
			rm -f __parsedPrimersF __parsedPrimersR forward_primers.txt reverse_primers.txt ${dirMeta}/primers_parsed.txt
			awk 'FNR>1' ${primers} > __tempPrimers
			while read row
			do
				#echo -e "ROW:     ${row}"
				#extract the forward primers: 1. cut the forward column and replace ; by new lines (in case of staggered primers)  2. add a sequences lengths  3. extract the smallest primer and put it in a file primers_parsed.txt
				echo -e "${row}" | cut -f${forwardColNumber} | sed "s/;/\n/g"  > __forward
				paste __forward <(awk '{ print length($0); }' __forward) | sort -n -k2,2  > _forward_primers.txt
				cat _forward_primers.txt >> forward_primers.txt
				head -n1 _forward_primers.txt | cut -f1 >> __parsedPrimersF
				#Do the same with reverse primers
				echo -e "${row}" | cut -f${reverseColNumber} | sed "s/;/\n/g"  > __reverse
				paste __reverse <(awk '{ print length($0); }' __reverse) | sort -n -k2,2  > _reverse_primers.txt
				cat _reverse_primers.txt >> reverse_primers.txt
				head -n1 _reverse_primers.txt | cut -f1 >> __parsedPrimersR
			done<__tempPrimers
			sort forward_primers.txt | uniq |  perl -ne '(!/^\s+$/)&&print' | sort -n -k2,2 | sed "1s/^/forwardPrimerSequence\tlength\n/" > _temp
			mv _temp forward_primers.txt
			sort reverse_primers.txt | uniq |  perl -ne '(!/^\s+$/)&&print' | sort -n -k2,2 | sed "1s/^/reversePrimerSequence\tlength\n/" > _temp
			mv _temp reverse_primers.txt
			paste __parsedPrimersF __parsedPrimersR | sort | uniq |  perl -ne '(!/^\s+$/)&&print' | sed "1s/^/forward\treverse\n/" > ${dirMeta}/primers_parsed.txt
			rm -f __parsedPrimersF __parsedPrimersR __forward __reverse __tempPrimers _forward_primers.txt _reverse_primers.txt
		else
			echo -e "Couldn't find the forward and/or reverse columns in the file: ${primers}. Be sure to have the file with the right format.\nExample:\nforward\treverse\nGTGCCAGCMGCCGCGGTAA;TGTGCCAGCMGCCGCGGTAA;ACGTGCCAGCMGCCGCGGTAA;CTAGTGCCAGCMGCCGCGGTAA\tGGACTACHVGGGTWTCTAAT;TGGACTACHVGGGTWTCTAAT;ACGGACTACHVGGGTWTCTAAT;CTAGGACTACHVGGGTWTCTAAT\n"
		fi

	else
	   echo -e "\n-----ERROR\nCouldn't find the primer file: ${dirMeta}/primers.txt.\nIf you don't have the primers, change the pipe.ini file, otherwise put all primers sequences into a file (see path and file name above). Example:\nforward\treverse\nGTGCCAGCMGCCGCGGTAA;TGTGCCAGCMGCCGCGGTAA;ACGTGCCAGCMGCCGCGGTAA;CTAGTGCCAGCMGCCGCGGTAA\tGGACTACHVGGGTWTCTAAT;TGGACTACHVGGGTWTCTAAT;ACGGACTACHVGGGTWTCTAAT;CTAGGACTACHVGGGTWTCTAAT\n" >&2
	   exit 1
	fi
fi


cd ${dirreads}
#unzip raw reads if needed
checkFastqs=$(ls -1 *.* | rev | cut -d"." -f1 | rev | grep "fastq" | wc -l)
if [ ${checkFastqs} -eq 0 ]; then
	for i in *.gz
	do
		removeFastq="YES"
		newfile=$(basename $i .fastq.gz)
		gunzip -c $i > $newfile.fastq
	done
fi

mkdir -p ${dir0}/run
cd ${dir0}/run

if [ "${primerSelectionBypass}"  == "YES" ]; then
	echo -e "By-passing primer removal step."
	mkdir -p ${dir0}/run/fastqsWithPrimers
	cd ${dir0}/raw_reads
	#remove fastq in raw reads directory if they were zipped
	if [ "${removeFastq}"  == "YES" ]; then
		mv ${dir0}/raw_reads/*.fastq ${dir0}/run/fastqsWithPrimers
	else
		cd ${dir0}/raw_reads
		cp *.fastq ${dir0}/run/fastqsWithPrimers
		for i in *.fastq; do echo -e "zipping --> ${i}";gzip $i; done
	fi
	mkdir -p ${successDir}
	mkdir -p ${dir0}/run/Summary
	touch ${successDir}/primerSelection.ok
	echo -e "\n---\nNo primer were removed at the user's choice. Script is exiting normally!"
	exit 0
fi


#Separate primers in groups (in case you have staggered primers):
cd ${dir0}/fastqsWithPrimers
#gather all primers into one 
cat <(awk 'FNR>1' forward_primers.txt) <(awk 'FNR>1' reverse_primers.txt) | perl -ne '(!/^\s+$/)&&print' | sort -n -k2,2  > _temp1
i=1
#check whether the sequence observed is one of the main primer (meaning not stagerred)
awk 'FNR>1' ${dirMeta}/primers_parsed.txt | sed "s/\t/\n/" | perl -ne '(!/^\s+$/)&&print' >__mainPrimers
awk 'FNR>1' ${primers} | sed "s/\t/\n/" | perl -ne '(!/^\s+$/)&&print' > __allPrimers
while read seq length
do
	checkIfSeqInFile=$(grep "^${seq}$" __mainPrimers | wc -l)
	if [ ${checkIfSeqInFile} -gt 0 ];then
		rm -f _restOfSequences
		
		#I'll extract the sequences from primers.txt file corresponding to that main primer. There are 3 posiibilities (1 primer per line, primer in the beginning of a line or primer in between ;)
		grep "^${seq}$" __allPrimers >>_restOfSequences |true
		grep "^${seq};" __allPrimers >>_restOfSequences |true
		grep ";${seq};" __allPrimers >>_restOfSequences |true
		#remove ; character in case of staggered primers
		sed -i "s/;/\n/g" _restOfSequences
		echo -e "${seq}" > staggeredGroup_${i}.txt
		grep -v "^${seq}$" _restOfSequences | grep "${seq}$" >> staggeredGroup_${i}.txt | true
		((i = i + 1))
	fi
done<_temp1

rm -f _temp1 _restOfSequences __mainPrimers __allPrimers

#Now we have several group of staggered primers (or not)
groupList=$(ls -1 staggeredGroup_*.txt)
for group in ${groupList}
do
	paste ${group} <(awk '{ print length($0); }' ${group}) | sort -n -k2,2 > _temp1
	smallest=$(sort -n -k2,2 _temp1 | cut -f2 | head -n1)
	rm -f _temp2
	while read seq length
	do
		numberOfAddedBases=$(echo "scale=3; ${length} - ${smallest}" | bc)
		echo -e "${seq}\t${length}\t${numberOfAddedBases}" >> _temp2
	done<_temp1
	mv _temp2 parsed_${group}
done


ls -1 ${dirreads}/*.fastq | sort > __readsList


if [ ! -s "__readsList" ];
then
	echo -e "It seems there are no fastq in ${dirreads}" >&2
	exit 1
fi


mkdir -p ${dir0}/fastqsWithPrimers
mkdir -p ${dir0}/fastqsWithoutExactPrimers
cd ${dir0}/fastqsWithPrimers
rm -f *_noPrimerFound.fastq *_withPrimers.fastq primerSelection_summary.log
echo -ne "Read_Name\tTotal_reads" > headers_primerSelection_summary.log
rm -f primerSelection_summary.log

while read readfile;
do
	fullname=$(basename "$readfile")
	filename="${fullname##*/}"
	sampleReadName=$(basename "$filename" | cut -d. -f1)
	echo -ne "${sampleReadName}\t" >>primerSelection_summary.log
	ln -nsf ${readfile} ${filename}
	echo -e "Inspecting file: ${filename}"
	i=1
	initialSeqNumber=$(awk '{s++}END{print s/4}' ${sampleReadName}.fastq)
	echo -ne "${initialSeqNumber}\t" >>primerSelection_summary.log
	mv ${filename} _inputSeqs
	for group in ${dir0}/fastqsWithPrimers/parsed_staggeredGroup_*.txt
	do
		while read seq length trim
		do
			echo -e "primer: ${seq}    |     base(s) to trim: ${trim}"
			bash ${dirpipe}/extract_reads_with_primers.sh _inputSeqs ${seq} ${trim}
			#echo -e "NORMAL EXIT" && exit
			if [[ -s _inputSeqs_noPrimerFound.fastq ]];then
				mv _inputSeqs_noPrimerFound.fastq _inputSeqs
			fi
			if [[ -s _inputSeqs_withPrimers.fastq ]];then 
				NumbSeqs=$(awk '{s++}END{print s/4}' _inputSeqs_withPrimers.fastq)
				echo -ne "${NumbSeqs}\t">>primerSelection_summary.log
				echo -ne "\t${seq}" >> headers_primerSelection_summary.log
				mv _inputSeqs_withPrimers.fastq _${i}_withPrimers.fastq
			else
				echo -ne "0\t" >>primerSelection_summary.log
				echo -ne "\t${seq}" >> headers_primerSelection_summary.log
			fi
			((i = i + 1))
			unset NumbSeqs
		done<${group}
	done
	cat _*_withPrimers.fastq > ${sampleReadName}_withPrimers.fastq |true
	rm -f _*_withPrimers.fastq _temp1 _restOfSequences
	if [[ -s ${sampleReadName}_withPrimers.fastq ]];then
		NumbSeqs=$(awk '{s++}END{print s/4}' ${sampleReadName}_withPrimers.fastq)
		percentageKept=$(echo "scale=3; ${NumbSeqs} / ${initialSeqNumber} * 100" | bc)
		echo -ne "${NumbSeqs}\t${percentageKept}">>primerSelection_summary.log
		echo -ne "\tReadsWithExactPrimer\tPercentageWithExactPrimers" >> headers_primerSelection_summary.log
		echo -e "" >> headers_primerSelection_summary.log
	else
		echo -ne "0\t0" >>primerSelection_summary.log
		echo -ne "\ttReadsWithExactPrimer\tPercentageWithExactPrimers" >> headers_primerSelection_summary.log
		echo -e "" >> headers_primerSelection_summary.log
	fi
	echo -e "" >> primerSelection_summary.log
	mv _inputSeqs ${dir0}/fastqsWithoutExactPrimers/${sampleReadName}_noExactPrimerFound.fastq
done<__readsList

rm -f parsed_staggeredGroup_*.txt staggeredGroup_*.txt

head -n1 headers_primerSelection_summary.log > _temp
cat _temp primerSelection_summary.log > _temp2
mv _temp2 primerSelection_summary.log

rm -f *tatistics.log _temp _temp2 headers_primerSelection_summary.log



#CLEAN PE READS (for example if a fwd read had a primer but not the reverse, we need to remove that read)

cd ${dir0}/fastqsWithPrimers
echo -e "\n-----\nCleaning the Paired-end fastqs"
for i in *R1_withPrimers.fastq; 
do 
	core=$(basename $i R1_withPrimers.fastq)
	bash ${dirpipe}/fastq_PE_fixer.sh ${core}R1_withPrimers.fastq ${core}R2_withPrimers.fastq
done



cd ${dir0}/fastqsWithoutExactPrimers
echo -e "\n-----\nCleaning the Paired-end fastqs"
for i in *R1_noExactPrimerFound.fastq; 
do 
	core=$(basename $i R1_noExactPrimerFound.fastq)
	#echo -e "${i}     |    ${core}"
	bash ${dirpipe}/fastq_PE_fixer.sh ${core}R1_noExactPrimerFound.fastq ${core}R2_noExactPrimerFound.fastq
	mv ${core}R1_noExactPrimerFound_ordered.fastq ${core}R1.fastq
	mv ${core}R2_noExactPrimerFound_ordered.fastq ${core}R2.fastq
	rm -f ${core}R1_noExactPrimerFound.fastq ${core}R2_noExactPrimerFound.fastq
	gzip ${core}R1.fastq
	gzip ${core}R2.fastq
done




cd ${dir0}/fastqsWithPrimers
rm -f *_withPrimers.fastq

echo -e "Add statistics of processed files"
sed -i "1s/$/\tPEFixing\tPercentageReadsKept/" primerSelection_summary.log
i=2
while read readfile;
do
	fullname=$(basename "$readfile")
	filename="${fullname##*/}"
	outputname=$(basename "$filename" | cut -d. -f1)
	NumbSeqs=$(awk '{s++}END{print s/4}' ${outputname}_withPrimers_ordered.fastq)
	#now the percentage of reads kept
	initialSeqNumber=$(awk '{s++}END{print s/4}' ${dir0}/raw_reads/${outputname}.fastq)
	percentageKept=$(echo "scale=3; ${NumbSeqs} / ${initialSeqNumber} * 100" | bc)
	mv ${outputname}_withPrimers_ordered.fastq ${outputname}.fastq
	sed -i "${i}s/$/\t${NumbSeqs}\t${percentageKept}/" primerSelection_summary.log
	unset NumbSeqs
	((i = i + 1))
done<__readsList



rm -f __readsList


#CLEAN WORKDIR

mv ${dir0}/fastqsWithPrimers ${dir0}/run/
mv ${dir0}/fastqsWithoutExactPrimers ${dir0}/run/
mkdir -p ${successDir}
mkdir -p ${dir0}/run/Summary
cp ${dir0}/run/fastqsWithPrimers/primerSelection_summary.log ${dir0}/run/Summary/primer_selection.txt
touch ${successDir}/primerSelection.ok
#remove fastq in raw reads directory if they were zipped
if [ "${removeFastq}"  == "YES" ]; then
	rm -f ${dir0}/raw_reads/*.fastq
else
	cd ${dir0}/raw_reads
	for i in *.fastq; do echo -e "zipping --> ${i}";gzip $i; done
fi


echo -e "\n---\nprimerSelection.sh exiting normally.\nWell done!"

