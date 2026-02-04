#!/bin/bash

#############################################################################################################################
# $1:     [ Designated_samples | All_samples ]
# $2:     read_file_type [ sra | fastq ]
# $3:     raw_data_path
# $4:     cellranger index
# $5:     cellranger localcores
# $6:     cellranger localmema
# $7:     Selected sample list [ sample1-sample2-sample3 ]
# $8:     10X Genomics library construction strategies [ 3_end_v1 | 3_end_v2 | 3_end_v3 | 3_end_v3.1 | 5_end_v1 | 5_end_v2 ]
#############################################################################################################################

# Build the path of intermediate output file and expression matrixes.
single_project_path=$(dirname "$3")
intermediate_output=${single_project_path}/2_output
matrix_path=${single_project_path}/3_expression_result
mkdir -p $intermediate_output && mkdir -p $matrix_path

# The path of core scripts.
CORE_SCRIPTS_PATH=$(cd $(dirname $0) && pwd)

# # Build the barcode array used to distinguish R1 (e.g. 150bp) and R2 (e.g. 150bp) according to different library construction strategies.
# if [[ "$8" == "3_end_v3" ]] || [[ "$8" == "3_end_v3.1" ]];then
#     # barcodes for Single Cell 3' v3, or Single Cell 3' v3.1, or Single Cell 3' HT v3.1
#     barcode_array=($(cat ${CORE_SCRIPTS_PATH}/config/10X_Genomics_Barcodes/3M-february-2018.txt | awk '{print $1}'))
# elif [[ "$8" == "3_end_v2" ]] || [[ "$8" == "5_end_v1" ]] || [[ "$8" == "5_end_v2" ]];then
#     # barcodes for Single Cell 3' v2, or Single Cell 5' v1 or v2, or Single Cell 5' HT v2
#     barcode_array=($(cat ${CORE_SCRIPTS_PATH}/config/10X_Genomics_Barcodes/737K-august-2016.txt | awk '{print $1}'))
# elif [[ "$8" == "3_end_v1" ]];then
#     # barcodes for Single Cell 3' v1
#     barcode_array=($(cat ${CORE_SCRIPTS_PATH}/config/10X_Genomics_Barcodes/737K-april-2014_rc.txt | awk '{print $1}'))
#     # The library construction of 10X genomics is not clear. Maybe it is 3' or 5'.
# elif [[ "$8" == "3_end" ]] || [[ "$8" == "5_end" ]];then
#     barcode_array=($(cat ${CORE_SCRIPTS_PATH}/config/10X_Genomics_Barcodes/*.txt | awk '{print $1}'))
# else
#     echo "Check the 10X Genomics library strategy you used or the version of cellranger."
# fi

# The library construction of 10X genomics is not clear. Maybe it is 3' or 5'.
barcode_array=($(cat ${CORE_SCRIPTS_PATH}/config/10X_Genomics_Barcodes/*.txt | awk '{print $1}'))

# When you run all samples in your project, this control flow is chosen.
if [[ "$1" == "All_samples" ]];then
    for Sample in $3/*
    do
        # Seize the sample name.
        sample_name=$(echo "$Sample" | awk -F "/" '{print $NF}')
        echo "The sample $sample_name exists."

        if [ -d $Sample ];then
            for Read in $Sample/*
            do
                if [ -f "$Read" ] && [ "$2" = "sra" ];then
                    # Seize the read file name.
                    run_file=$(echo "$Read" | awk -F "/" '{print $NF}')
                    echo "The read file $run_file exists."

                    # Check whether SRR files exist, and no operations are performed on fastq files.
                    if [ "${Read##*.}"x != "fastq"x ] || [ "${Read##*.}"x != "fq"x ];then
                        read_file=$(echo "$Read" | awk -F "/" '{print $NF}')
               
                        # Convert sra file to fastq file.
                        fasterq-dump --split-files --include-technical -e $5 -O $Sample -t $Sample $Read
                    fi

                elif [ -f "$Read" ] && [ "$2" = "fastq" ];then
                    # Convert .fq into .fastq
                    if [ "${Read##*.}"x = "fq"x ];then
                        mv $Read ${Read%.*}.fastq
                    fi
                fi  
            done

            # The fastq files whose name have been changed will be put in the dir.
            mkdir -p $intermediate_output/$sample_name

            if [[ "$2" == "sra" ]];then
                for fastq_file in `ls -1 $Sample/*.fastq`
                do
                    if [ -f "$fastq_file" ];then
                        head -100 $fastq_file > ${fastq_file}.head

                        # Seize the prefix of the fastq.
                        sra_file=${fastq_file##*/}
                        fastq_prefix=${sra_file%.*}  

                        # Seize a read in each fastq file.
                        reads=`sed -n 2p $Sample/${sra_file}.head`

                        # Seize the length of the read.
                        read_length=`sed -n 2p $Sample/${sra_file}.head | wc -m`

                        # Seize the first 16 bases, which may be barcode
                        pseudo_barcode=${reads:0:16}
			echo $pseudo_barcode

                        # Sezie the 15 bases after the first base, for the first base may be 'N' sometimes.
                        # sub_pseudo_barcode=${reads:1:15}
			            # echo $sub_pseudo_barcode

                        # # Seize the first base, which may be 'N' sometimes.
                        # first_base=${reads:0:1}

                        # Determine whether a or some 'N' base(s) in the barcode of R1.
                        N_BASE="N"
                        
                        IF_N_BASE=$(echo $pseudo_barcode | grep "${N_BASE}")

                        if [ -n "$IF_N_BASE" ];then
                            # There is a 'N' base in the barcode of read1, determine whether the read is I1, I2, R1, R2

                            # Get the index of the "N" base.
                            N_BASE_INDEX=`echo $pseudo_barcode | awk -v param=${N_BASE} '{printf("%d\n",match($0,param))}'`

                            # Change the base into "N" base in barcode configuration files according to the "N" base in the barcode in fastq file (With the same index).
                            barcode_array_contain_N=($(cat ${CORE_SCRIPTS_PATH}/config/10X_Genomics_Barcodes/*.txt | awk -F "" -v col=$N_BASE_INDEX 'BEGIN{OFS=""}{$col="N"; print $0}'))

                            if [[ $read_length -ge 20 ]] && [[ ! "${barcode_array_contain_N[@]}" =~ "$pseudo_barcode" ]];then
                                mv $Sample/$sra_file $Sample/${fastq_prefix}_S1_L001_R2_001.fastq
                            elif [[ $read_length -ge 20 ]] && [[ "${barcode_array_contain_N[@]}" =~ "$pseudo_barcode" ]];then
                                mv $Sample/$sra_file $Sample/${fastq_prefix}_S1_L001_R1_001.fastq
                            else
                                mv $Sample/$sra_file $Sample/${fastq_prefix}_S1_L001_I1_001.fastq
                            fi

                        else
                            # There is no 'N' bases in the barcode of read1, determine whether the read is I1, I2, R1, R2
                            if [[ $read_length -ge 20 ]] && [[ ! "${barcode_array[@]}" =~ "$pseudo_barcode" ]];then
                                mv $Sample/$sra_file $Sample/${fastq_prefix}_S1_L001_R2_001.fastq
                            elif [[ $read_length -ge 20 ]] && [[ "${barcode_array[@]}" =~ "$pseudo_barcode" ]];then
                                mv $Sample/$sra_file $Sample/${fastq_prefix}_S1_L001_R1_001.fastq
                            else
                                mv $Sample/$sra_file $Sample/${fastq_prefix}_S1_L001_I1_001.fastq
                            fi
                        fi

                    fi
                    
                done

            elif [[ "$2" == "fastq" ]];then
                echo "Make sure that the name of fastq files as shown below: "
                echo "MySample_S1_L001_I1_001.fastq (optional)"
                echo "MySample_S1_L001_I2_001.fastq (optional)"
                echo "MySample_S1_L001_R1_001.fastq"
                echo "MySample_S1_L001_R2_001.fastq"

            else
                echo "Please check the parameter --ReadType/-rt, there is something wrong!"
            fi

            # If there is more one sra file in a Sample path, combine them into one.
            cat ${Sample}/*_S1_L001_R2_001.fastq > $intermediate_output/$sample_name/${sample_name}_S1_L001_R2_001.fastq
            cat ${Sample}/*_S1_L001_R1_001.fastq > $intermediate_output/$sample_name/${sample_name}_S1_L001_R1_001.fastq

	    # If I1 and I2 exist at the same time, it is difficult to distinguish them. Fortunately, they are optional and we decide to delete them.
            fastq_I1_file=$(ls ${Sample}/*_S1_L001_I1_001.fastq 2> /dev/null | wc -l)
            if [ "$fastq_I1_file" != "0" ];then
		rm ${Sample}/*_S1_L001_I1_001.fastq
		#cat ${Sample}/*_S1_L001_I1_001.fastq > $intermediate_output/$sample_name/${sample_name}_S1_L001_I1_001.fastq
            fi

            # If I1 and I2 exist at the same time, it is difficult to distinguish them. Fortunately, they are optional and we decide to delete them.
           # if [[ ! -s $intermediate_output/$sample_name/${sample_name}_S1_L001_I1_001.fastq ]];then
           #     rm $intermediate_output/$sample_name/${sample_name}_S1_L001_I1_001.fastq
           # fi

            # Enter the matrix path.
            cd $matrix_path

            # Generating the matrix...
            cellranger count \
            --id=${sample_name} \
            --sample=${sample_name} \
            --transcriptome=$4 \
            --fastqs=$intermediate_output/${sample_name} \
            --localcores=$5 \
            --localmem=$6

        fi
    done
    
# If some samples run failed because of the limitation of cores, or mem, or something else, you can re-run the failed samples only.
# And the failed samples in "03.expression" should be deleted before re-runnning.
elif [[ "$1" == "Designated_samples" ]];then

    # Enter the matrix path.
    cd $matrix_path

    # Construct the array of failed samples.
    samples=(`echo "$7" | awk '{len=split($0,sample_list,"-");for(i=1;i<=len;i++) print sample_list[i]}'`)

    for each_sample in "${samples[@]}"
    do
        cellranger count \
        --id=${each_sample} \
        --sample=${each_sample} \
        --transcriptome=$4 \
        --fastqs=$intermediate_output/${each_sample} \
        --localcores=$5 \
        --localmem=$6
    done
fi
