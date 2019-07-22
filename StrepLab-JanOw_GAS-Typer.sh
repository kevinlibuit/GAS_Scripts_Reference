#!/bin/bash -l

temp_path=$(pwd)
export PATH=$PATH:$temp_path

###Load Modules###
source /usr/share/modules/init/bash 
module load perl/5.22.1
module load ncbi-blast+/2.2.29
module load BEDTools/2.17.0
module load freebayes/0.9.21
module load prodigal/2.60
#module load cutadapt/1.8
module load cutadapt/1.8.3
module load srst2/0.1.7

###This script is called for each job in the qsub array. The purpose of this code is to read in and parse a line of the job-control.txt file
###created by 'StrepLab-JanOw_GAS-wrapr.sh' and pass that information, as arguments, to other programs responsible for various parts of strain
###characterization (MLST, emm type and antibiotic drug resistance prediction).

sample_name=$1
r1=$2
r2=$3
out_dir=$4
allDB_dir=/GAS_Scripts_Reference/GAS_Reference_DB

###Pre-Process Paired-end Reads###
mkdir -p $out_dir/cut_adapt_output/$sample_name
r1_trimd=cutadapt_${sample_name}_R1.fastq.gz
r2_trimd=cutadapt_${sample_name}_R2.fastq.gz
cutadapt -b AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -q 20 --minimum-length 50 --paired-output ${out_dir}/cut_adapt_output/${sample_name}/temp2.fastq.gz -o ${out_dir}/cut_adapt_output/${sample_name}/temp1.fastq.gz $r1 $r2
cutadapt -b AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -q 20 --minimum-length 50 --paired-output ${out_dir}/cut_adapt_output/${sample_name}/$r1_trimd -o ${out_dir}/cut_adapt_output/${sample_name}/$r2_trimd ${out_dir}/cut_adapt_output/${sample_name}/temp2.fastq.gz ${out_dir}/cut_adapt_output/${sample_name}/temp1.fastq.gz
rm ${out_dir}/cut_adapt_output/${sample_name}/temp1.fastq.gz
rm ${out_dir}/cut_adapt_output/${sample_name}/temp2.fastq.gz

###Call MLST###
export SRST2_SAMTOOLS="/samtools/samtools-0.1.18/samtools" && srst2 --samtools_args "\-A" --mlst_delimiter '_' --input_pe "${out_dir}/cut_adapt_output/${sample_name}/$r1_trimd" "${out_dir}/cut_adapt_output/${sample_name}/$r2_trimd" --forward "_R1" --reverse "_R2" --output "${out_dir}/GAS_output/${sample_name}/${sanple_name}_MLST" --save_scores --mlst_db "$allDB_dir/Streptococcus_pyogenes.fasta" --mlst_definitions "$allDB_dir/spyogenes.txt" --min_coverage 99.999
###Check and extract new MLST alleles###
MLST_allele_checkr.pl "$out_nameMLST"__mlst__Streptococcus_pyogenes__results.txt "$out_nameMLST"__*.Streptococcus_pyogenes.sorted.bam "$allDB_dir/Streptococcus_pyogenes.fasta"

###Call emm Type###
module unload perl/5.22.1
module load perl/5.16.1-MT
emm_typer.pl -1 "${out_dir}/cut_adapt_output/${sample_name}/$r1_trimd" -2 "${out_dir}/cut_adapt_output/${sample_name}/$r2_trimd" -r "$allDB_dir" -n "$sample_name" -o "${out_dir}"
PBP-Gene_Typer.pl -1 "$readPair_1" -2 "$readPair_2" -r "$allDB_dir/GAS_bLactam_Ref.fasta" -n "$just_name" -s GAS -p 2X
module unload perl/5.16.1-MT
module load perl/5.22.1

###Call GAS Misc Resistance###
GAS_Res_Typer.pl -1 "${out_dir}/cut_adapt_output/${sample_name}/$r1_trimd" -2 "${out_dir}/cut_adapt_output/${sample_name}/$r2_trimd" -d "$allDB_dir" -r GAS_Res_Gene-DB_Final.fasta -n "$sample_name" -o "${out_dir}"
#GAS_Target2MIC.pl TEMP_Res_Results.txt "$just_name" TEMP_pbpID_Results.txt

###Type Surface and Secretory Proteins###
#GAS_Features_Typer.pl -1 "${out_dir}/cut_adapt_output/${sample_name}/$r1_trimd" -2 "${out_dir}/cut_adapt_output/${sample_name}/$r2_trimd" -d "$allDB_dir" -f GAS_features_Gene-DB_Final.fasta -n "$sample_name" -o "${out_dir}/GAS_output/${sample_name}"


###Output the emm type/MLST/drug resistance data for this sample to it's results output file###
#tabl_out="TABLE_Isolate_Typing_results.txt"
#bin_out="BIN_Isolate_Typing_results.txt"
#contamination_level=10
#printf "$just_name\t" >> "$tabl_out"
#printf "$just_name," >> "$bin_out"
###EMM TYPE OUTPUT###
#emm_out="NF"
#while read -r line
#do
#    if [[ -n "$line" ]]
#    then
#        justTarget=$(echo "$line" | awk -F"\t" '{print $1}')
#        if [[ "$emm_out" == "NF" ]]
#        then
#            emm_out="$justTarget#"
        #else
        #    emm_out="$emm_out;$justTarget"
        #fi
    #fi
#done <<< "$(sed 1d *__emm-Type__Results.txt)"
#printf "$emm_out\t" >> "$tabl_out"
#printf "$emm_out," >> "$bin_out"
####MLST OUTPUT###
#sed 1d "$out_nameMLST"__mlst__Streptococcus_pyogenes__results.txt | while read -r line
#do
#    MLST_tabl=$(echo "$line" | cut -f2-9)
#    echo "MLST line: $MLST_tabl\n";
#    printf "$MLST_tabl\t" >> "$tabl_out"
#    MLST_val=$(echo "$line" | awk -F" " '{print $2}')
#    printf "$MLST_val," >> "$bin_out"
#done #< "$out_nameMLST"__mlst__Streptococcus_pyogenes__results.txt
##tail -n+2 "$out_nameMLST"__mlst__Streptococcus_pyogenes__results.txt | cut -f2-9 >> "$tabl_out"

###Features Targets###
#while read -r line
#do
#    FEAT_targ=$(echo "$line" | cut -f2)
#    printf "$FEAT_targ\t" >> "$tabl_out"
#done < TEMP_protein_Results.txt
#
####PBP_ID Output###
##justPBPs="NF"
#sed 1d TEMP_pbpID_Results.txt | while read -r line
#do
#    if [[ -n "$line" ]]
#    then
#	justPBPs=$(echo "$line" | awk -F"\t" '{print $2}')
#    fi
#    printf "$justPBPs\t" >> "$tabl_out"
#done

###Resistance Targets###
#while read -r line
#do
#    #RES_targ=$(echo "$line" | cut -f2)
#    #printf "$RES_targ\t" >> "$tabl_out"
#    printf "$line\t" | tr ',' '\t' >> "$tabl_out"
#done < RES-MIC_"$just_name"
#printf "\n" >> "$tabl_out"

#cat BIN_Features_Results.txt | sed 's/$/,/g' >> "$bin_out"
#cat BIN_Res_Results.txt >> "$bin_out"
#printf "\n" >> "$bin_out"

###Remove Temporary Files###
#rm cutadapt*.fastq
#rm *.pileup
#rm *.bam
#rm *.sam
#rm TEMP*

###Unload Modules###
#module unload perl/5.22.1
#module unload ncbi-blast+/2.2.29
#module unload BEDTools/2.17.0
#module unload freebayes/0.9.21
#module unload prodigal/2.60
##module unload cutadapt/1.8
#module unload cutadapt/1.8.3
#module unload srst2/0.1.7
