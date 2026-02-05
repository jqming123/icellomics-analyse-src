#!/bin/bash

####################################################################
# $1:     build_library_type [ Bulk | Smart-seq2 ]                                       
# $2:     [ Designated_samples | All_samples ]
# $3:     Selected sample list [ sample1-sample2-sample3 ]
# $4:     sequencing type [ single | pair ]
# $5:     read file type [ sra | fastq ]
# $6:     raw data path
# $7:     hisat2 index
# $8:     RSEM index
# $9:     bed file
# =================================================================
#                       fasterq-dump parameters
# ${10}:  fasterq-dump_thread_num, default: 12
# =================================================================
#                          fastp parameters
# ${11}:  -q default: 20
# ${12}:  -u default: 40
# ${13}:  -l default: 50
# ${14}:  -W default: 4
# ${15}:  -M default: 20
# ${16}:  -w default: 12
# Note: Other parameters are default.                               
# =================================================================
#                          hisat2 parameters
# ${17}:  -p default: 12
# =================================================================
#                          samtools parameters
# ${18}:  -@ default: 12
# =================================================================
#                  rsem-calculate-expression parameters
# ${19}:  -p default: 12
# ${20}:  --star-path
# ${21}:  Whether generating matrixes when "Designated_samples"
#         [ g_matrix | n_matrix]
####################################################################

# Build the path of intermediate output file and expression matrixes.
single_project_path=$(dirname "$6")
intermediate_output=${single_project_path}/2_output
matrix_path=${single_project_path}/3_expression_result
mkdir -p $intermediate_output
mkdir -p $matrix_path

# When the samples are too many, or several samples needs to be re-run, or something else, you can choose specific to run, not all samples.
if [[ "$2" == "Designated_samples" ]];then
    for Sample in $6/*
    do
        # Seize each sample name in the whole sample list.
        sample_name=$(echo "$Sample" | awk -F "/" '{print $NF}')

        # Determine if the sample list you selected in the whole sample list.
        samples=(`echo "$3" | awk '{len=split($0,sample_list,"-");for(i=1;i<=len;i++) print sample_list[i]}'`)
        for each_sample in "${samples[@]}"
        do
            if [ "$sample_name" = "$each_sample" ];then
                echo "The $sample_name you selescted is in the sample list."
                if [ -d "$Sample" ];then
                    echo "The sample $sample_name exists."

                    for Read in $Sample/*
                    do
                        # If the read file is "sra" format, convert it to "fastq" format.
                        if [[ -f "$Read" ]] && [[ "$5" == "sra" ]];then
                            if [ "${Read##*.}"x = "sra"x ];then
                                read_file=$(echo "$Read" | awk -F "/" '{print $NF}')
                                echo "The read file $read_file exist."
                                fasterq-dump --split-3 -e ${10} -O $Sample -t $Sample $Read
                            fi

                        elif [ -f "$Read" ] && [ "$5" = "fastq" ];then

                            if [ "${Read##*.}"x = "fq"x ];then
                                mv $Read ${Read%.*}.fastq
                            fi
                        fi
                    done

                    mkdir -p $intermediate_output/$sample_name/
                    # Handle paired-end reads.
                    if [[ "$4" == "pair" ]];then
                        # If there are more than two read1 files (fastq format) or two read2 files (fastq format) in a sample, combine them into a single read1 file (fastq format) or read2 file (fastq format).
                        cat $Sample/*1.fastq > $intermediate_output/$sample_name/${sample_name}_1.fastq
                        cat $Sample/*2.fastq > $intermediate_output/$sample_name/${sample_name}_2.fastq
                        
                        # Quality Control (QC) by "fastp".
                        fastp -q ${11} -u ${12} -l ${13} -g -x -r -W ${14} -M ${15} -w ${16} -i $intermediate_output/$sample_name/${sample_name}_1.fastq -o $intermediate_output/$sample_name/${sample_name}_clean_1.fastq -I $intermediate_output/$sample_name/${sample_name}_2.fastq -O $intermediate_output/$sample_name/${sample_name}_clean_2.fastq -h $intermediate_output/$sample_name/${sample_name}_report.html -j $intermediate_output/$sample_name/${sample_name}_fastp.json 2> $intermediate_output/$sample_name/${sample_name}_fastp _report.txt

                        # Quality Control (QC) by "hisat2".
                        hisat2 -p ${17} -x $7 -1 $intermediate_output/$sample_name/${sample_name}_clean_1.fastq -2 $intermediate_output/$sample_name/${sample_name}_clean_2.fastq -S $intermediate_output/$sample_name/${sample_name}_Alignment-unsorted.sam 2> $intermediate_output/$sample_name/${sample_name}_hisat2_Mapping_Rate.txt
                    
                    # Handle single-end reads.
                    elif [[ "$4" == "single" ]];then
                        cat $Sample/*.fastq > $intermediate_output/$sample_name/${sample_name}.fastq

                        fastp -q ${11} -u ${12} -l ${13} -g -x -r -W ${14} -M ${15} -w ${16} -i $intermediate_output/$sample_name/${sample_name}.fastq -o $intermediate_output/$sample_name/${sample_name}_clean.fastq -h $intermediate_output/$sample_name/${sample_name}_report.html -j $intermediate_output/$sample_name/${sample_name}_fastp.json 2> $intermediate_output/$sample_name/${sample_name}_fastp _report.txt

                        hisat2 -p ${17} -x $7 -U $intermediate_output/$sample_name/${sample_name}_clean.fastq -S $intermediate_output/$sample_name/${sample_name}_Alignment-unsorted.sam 2> $intermediate_output/$sample_name/${sample_name}_hisat2_Mapping_Rate.txt
                    fi

                    # Convert sam file to bam file.
                    samtools view -b -@ ${18} -S $intermediate_output/$sample_name/${sample_name}_Alignment-unsorted.sam -o $intermediate_output/$sample_name/${sample_name}_Alignment-unsorted.bam

                    # Infer "Strand-Specific" by RSeQC.
                    infer_experiment.py -r $9 -i $intermediate_output/$sample_name/${sample_name}_Alignment-unsorted.bam > $intermediate_output/$sample_name/${sample_name}_RseQC.txt

                    parameter_1=`sed -n 5p $intermediate_output/$sample_name/${sample_name}_RseQC.txt | awk '{print \$7}'`
                    parameter_2=`sed -n 6p $intermediate_output/$sample_name/${sample_name}_RseQC.txt | awk '{print \$7}'`

                    # In the RSeQC output file, the threshold for the difference between the two strand is set to 0.2
                    decimal=`awk -v x="$parameter_1" -v y="$parameter_2" 'BEGIN{printf "%.0f\n",(x-y)*100}'`
                    echo "$decimal"
                    if [ $decimal -gt 20 ]; then
                        value=1
                    elif [ $decimal -gt -20 ]; then
                        value=0.5
                    else
                        value=0
                    fi
                    echo "The value is $value!"

                    # Use rsem-calculate-expression to claculate gene expression (paired-end or single end).
                    if [[ "$1" == "Bulk" ]] && [[ "$4" == "pair" ]];then
                        rsem-calculate-expression \
                        --keep-intermediate-files --temporary-folder $intermediate_output/$sample_name/${sample_name}_STAR_rsem \
                        --paired-end --forward-prob=$value -p ${19} --time \
                        --star --star-path ${20} \
                        --append-names --output-genome-bam \
                        --sort-bam-by-coordinate \
                        $intermediate_output/$sample_name/${sample_name}_clean_1.fastq $intermediate_output/$sample_name/${sample_name}_clean_2.fastq \
                        $8 \
                        $intermediate_output/$sample_name/${sample_name}_star_RSEM
                    elif [[ "$1" == "Smart-seq2" ]] && [[ "$4" == "pair" ]];then
                        rsem-calculate-expression  --single-cell-prior \
                        --keep-intermediate-files --temporary-folder $intermediate_output/$sample_name/${sample_name}_STAR_rsem \
                        --paired-end --forward-prob=$value -p ${19} --time \
                        --star --star-path ${20} \
                        --append-names --output-genome-bam \
                        --sort-bam-by-coordinate \
                        $intermediate_output/$sample_name/${sample_name}_clean_1.fastq $intermediate_output/$sample_name/${sample_name}_clean_2.fastq \
                        $8 \
                        $intermediate_output/$sample_name/${sample_name}_star_RSEM
                    elif [[ "$1" == "Bulk" ]] && [[ "$4" == "single" ]];then
                        rsem-calculate-expression \
                        --keep-intermediate-files --temporary-folder $intermediate_output/$sample_name/${sample_name}_STAR_rsem \
                        --forward-prob=$value -p ${19} --time \
                        --star --star-path ${20} \
                        --append-names --output-genome-bam \
                        --sort-bam-by-coordinate \
                        $intermediate_output/$sample_name/${sample_name}_clean.fastq \
                        $8 \
                        $intermediate_output/$sample_name/${sample_name}_star_RSEM
                    elif [[ "$1" == "Smart-seq2" ]] && [[ "$4" == "single" ]];then
                        rsem-calculate-expression --single-cell-prior \
                        --keep-intermediate-files --temporary-folder $intermediate_output/$sample_name/${sample_name}_STAR_rsem \
                        --forward-prob=$value -p ${19} --time \
                        --star --star-path ${20} \
                        --append-names --output-genome-bam \
                        --sort-bam-by-coordinate \
                        $intermediate_output/$sample_name/${sample_name}_clean.fastq \
                        $8 \
                        $intermediate_output/$sample_name/${sample_name}_star_RSEM
                    fi
                fi
            fi
        done
    done

    if [[ "${21}" == "g_matrix" ]];then
        # Generate gene expression matrixes of all samples by RSEM.
        PathSplit=(${6//\// })
        ProjectName=${PathSplit[-2]}
        echo "$ProjectName is partial prefix of matrix output files."

        rsem-generate-data-matrix $intermediate_output/*/*.genes.results > $matrix_path/${ProjectName}_GeneMat_rawCounts.txt
        rsem-generate-data-matrix $intermediate_output/*/*.isoforms.results > $matrix_path/${ProjectName}_TransMat_rawCounts.txt
        rsem-generate-data-matrix-TPM $intermediate_output/*/*.genes.results > $matrix_path/${ProjectName}_GeneMat_TPM.txt
        rsem-generate-data-matrix-TPM $intermediate_output/*/*.isoforms.results > $matrix_path/${ProjectName}_TransMat_TPM.txt
        rsem-generate-data-matrix-FPKM $intermediate_output/*/*.genes.results > $matrix_path/${ProjectName}_GeneMat_FPKM.txt
        rsem-generate-data-matrix-FPKM $intermediate_output/*/*.isoforms.results > $matrix_path/${ProjectName}_TransMat_FPKM.txt
    fi

# You can run all the samples at once.
elif [[ "$2" == "All_samples" ]];then
    for Sample in $6/*
    do
        sample_name=$(echo "$Sample" | awk -F "/" '{print $NF}')
        
        if [ -d "$Sample" ];then
            echo "The sample $sample_name exists."
            echo "$Sample"

            for Read in $Sample/*
            do
                if [ -f "$Read" ] && [ "$5" = "sra" ];then

                    if [ "${Read##*.}"x = "sra"x ];then
                        read_file=$(echo "$Read" | awk -F "/" '{print $NF}')
                        echo "The read file $read_file exist."

                        fasterq-dump --split-3 -e ${10} -O $Sample -t $Sample $Read
                    fi

                elif [ -f "$Read" ] && [ "$5" = "fastq" ];then

                    if [ "${Read##*.}"x = "fq"x ];then
                        mv $Read ${Read%.*}.fastq
                    fi
		        fi
            done

            mkdir -p $intermediate_output/$sample_name

            if [ "$4" = "pair" ];then
                cat $Sample/*1.fastq > $intermediate_output/$sample_name/${sample_name}_1.fastq
                cat $Sample/*2.fastq > $intermediate_output/$sample_name/${sample_name}_2.fastq

                fastp -q ${11} -u ${12} -l ${13} -g -x -r -W ${14} -M ${15} -w ${16} -i $intermediate_output/$sample_name/${sample_name}_1.fastq -o $intermediate_output/$sample_name/${sample_name}_clean_1.fastq -I $intermediate_output/$sample_name/${sample_name}_2.fastq -O $intermediate_output/$sample_name/${sample_name}_clean_2.fastq -h $intermediate_output/$sample_name/${sample_name}_report.html -j $intermediate_output/$sample_name/${sample_name}_fastp.json 2> $intermediate_output/$sample_name/${sample_name}_fastp _report.txt

                hisat2 -p ${17} -x $7 -1 $intermediate_output/$sample_name/${sample_name}_clean_1.fastq -2 $intermediate_output/$sample_name/${sample_name}_clean_2.fastq -S $intermediate_output/$sample_name/${sample_name}_Alignment-unsorted.sam 2> $intermediate_output/$sample_name/${sample_name}_hisat2_Mapping_Rate.txt
            
            elif [ "$4" = "single" ];then
                cat $Sample/*.fastq > $intermediate_output/$sample_name/${sample_name}.fastq

                fastp -q ${11} -u ${12} -l ${13} -g -x -r -W ${14} -M ${15} -w ${16} -i $intermediate_output/$sample_name/${sample_name}.fastq -o $intermediate_output/$sample_name/${sample_name}_clean.fastq -h $intermediate_output/$sample_name/${sample_name}_report.html -j $intermediate_output/$sample_name/${sample_name}_fastp.json 2> $intermediate_output/$sample_name/${sample_name}_fastp _report.txt

                hisat2 -p ${17} -x $7 -U $intermediate_output/$sample_name/${sample_name}_clean.fastq -S $intermediate_output/$sample_name/${sample_name}_Alignment-unsorted.sam 2> $intermediate_output/$sample_name/${sample_name}_hisat2_Mapping_Rate.txt
            fi

            samtools view -b -@ ${18} -S $intermediate_output/$sample_name/${sample_name}_Alignment-unsorted.sam -o $intermediate_output/$sample_name/${sample_name}_Alignment-unsorted.bam
            infer_experiment.py -r $9 -i $intermediate_output/$sample_name/${sample_name}_Alignment-unsorted.bam > $intermediate_output/$sample_name/${sample_name}_RseQC.txt

            parameter_1=`sed -n 5p $intermediate_output/$sample_name/${sample_name}_RseQC.txt | awk '{print \$7}'`
            parameter_2=`sed -n 6p $intermediate_output/$sample_name/${sample_name}_RseQC.txt | awk '{print \$7}'`

            decimal=`awk -v x="$parameter_1" -v y="$parameter_2" 'BEGIN{printf "%.0f\n",(x-y)*100}'`

            echo "$decimal"

            if [ $decimal -gt 20 ]; then
                value=1
            elif [ $decimal -gt -20 ]; then
                value=0.5
            else
                value=0
            fi

            echo "The value is $value!"

            if [[ "$1" == "Bulk" ]] && [[ "$4" == "pair" ]];then
                rsem-calculate-expression \
                --keep-intermediate-files --temporary-folder $intermediate_output/$sample_name/${sample_name}_STAR_rsem \
                --paired-end --forward-prob=$value -p ${19} --time \
                --star --star-path ${20} \
                --append-names --output-genome-bam \
                --sort-bam-by-coordinate \
                $intermediate_output/$sample_name/${sample_name}_clean_1.fastq $intermediate_output/$sample_name/${sample_name}_clean_2.fastq \
                $8 \
                $intermediate_output/$sample_name/${sample_name}_star_RSEM
            elif [[ "$1" == "Smart-seq2" ]] && [[ "$4" == "pair" ]];then
                rsem-calculate-expression  --single-cell-prior \
                --keep-intermediate-files --temporary-folder $intermediate_output/$sample_name/${sample_name}_STAR_rsem \
                --paired-end --forward-prob=$value -p ${19} --time \
                --star --star-path ${20} \
                --append-names --output-genome-bam \
                --sort-bam-by-coordinate \
                $intermediate_output/$sample_name/${sample_name}_clean_1.fastq $intermediate_output/$sample_name/${sample_name}_clean_2.fastq \
                $8 \
                $intermediate_output/$sample_name/${sample_name}_star_RSEM
            elif [[ "$1" == "Bulk" ]] && [[ "$4" == "single" ]];then
                rsem-calculate-expression \
                --keep-intermediate-files --temporary-folder $intermediate_output/$sample_name/${sample_name}_STAR_rsem \
                --forward-prob=$value -p ${19} --time \
                --star --star-path ${20} \
                --append-names --output-genome-bam \
                --sort-bam-by-coordinate \
                $intermediate_output/$sample_name/${sample_name}_clean.fastq \
                $8 \
                $intermediate_output/$sample_name/${sample_name}_star_RSEM
            elif [[ "$1" == "Smart-seq2" ]] && [[ "$4" == "single" ]];then
                rsem-calculate-expression --single-cell-prior \
                --keep-intermediate-files --temporary-folder $intermediate_output/$sample_name/${sample_name}_STAR_rsem \
                --forward-prob=$value -p ${19} --time \
                --star --star-path ${20} \
                --append-names --output-genome-bam \
                --sort-bam-by-coordinate \
                $intermediate_output/$sample_name/${sample_name}_clean.fastq \
                $8 \
                $intermediate_output/$sample_name/${sample_name}_star_RSEM
            fi

        fi

    done

    # Generate gene expression matrixes of all samples by RSEM.
    PathSplit=(${6//\// })
    ProjectName=${PathSplit[-2]}
    echo "$ProjectName is partial prefix of matrix output files."

    rsem-generate-data-matrix $intermediate_output/*/*.genes.results > $matrix_path/${ProjectName}_GeneMat_rawCounts.txt
    rsem-generate-data-matrix $intermediate_output/*/*.isoforms.results > $matrix_path/${ProjectName}_TransMat_rawCounts.txt
    rsem-generate-data-matrix-TPM $intermediate_output/*/*.genes.results > $matrix_path/${ProjectName}_GeneMat_TPM.txt
    rsem-generate-data-matrix-TPM $intermediate_output/*/*.isoforms.results > $matrix_path/${ProjectName}_TransMat_TPM.txt
    rsem-generate-data-matrix-FPKM $intermediate_output/*/*.genes.results > $matrix_path/${ProjectName}_GeneMat_FPKM.txt
    rsem-generate-data-matrix-FPKM $intermediate_output/*/*.isoforms.results > $matrix_path/${ProjectName}_TransMat_FPKM.txt
fi              
