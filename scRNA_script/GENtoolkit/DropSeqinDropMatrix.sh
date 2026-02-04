#!/bin/bash
#PBS -q 512G
#PBS -l mem=90gb,walltime=1000:00:00
#PBS -l nodes=1:ppn=$6
#HSCHED -s DropSeqinDrop+expression+matrix

###################################################################
# $1:     [ Drop-seq | inDrop_v1 | inDrop_v2 | inDrop_v3]
# $2:     [ Designated_samples | All_samples ]
# $3:     Selected sample list [ sample1-sample2-sample3 ]
# $4:     read_file_type [ sra | fastq ]
# $5:     raw_data_path
# $6:     dropTag_p
# $7:     star_index
# $8:     star_runThreadN
# $9:     dropEst_g
# ${10}:  dropReport_m
###################################################################

# Build the path of intermediate output file and expression matrixes.
single_project_path=$(dirname "$5")
intermediate_output=${single_project_path}/2_output
matrix_path=${single_project_path}/3_expression_result
mkdir -p $intermediate_output && mkdir -p $matrix_path

# Config files path.
script_path="$(cd "$(dirname "$0")" && pwd)"
#config_path="$(cd "$(dirname "$script_path")" && pwd)"
config_path="$(cd "$(dirname "$0")" && pwd)"

# When the samples are too many, or several samples needs to be re-run, or something else, you can choose specific to run, not all samples.
if [[ "$2" == "Designated_samples" ]];then

        for Sample in $5/* 
        do
                # Seize each sample name in the whole sample list.
                sample_name=$(echo "$Sample"| awk -F "/" '{print $NF}')
                echo "The sample $sample_name exists."

                # Determine if the sample list you selected in the whole sample list.
                samples=(`echo "$3" | awk '{len=split($0,sample_list,"-");for(i=1;i<=len;i++) print sample_list[i]}'`)
                for each_sample in "${samples[@]}"
                do
                        if [[ "$sample_name" =~ "$each_sample" ]];then
                                echo "The $sample_name you selescted is in the sample list."

                                if [ -d "$Sample" ];then 
                                        for Read in $Sample/*
                                        do
                                                if [ -f "$Read" ] && [ "$4" = "sra" ];then
                                                        if [ "${Read##*.}"x != "fastq"x ] || [ "${Read##*.}"x != "fq"x ];then
                                                                read_file=$(echo "$Read"| awk -F "/" '{print $NF}')
                                                                echo "The read file $read_file exist."

                                                                # Convert sra file to fastq file.
                                                                fastq-dump --split-files -O $Sample $Read
                                                        fi
                                                        
                                                elif [ -f "$Read" ] && [ "$4" = "fastq" ];then
                                                        if [ "${Read##*.}"x = "fq"x ];then
                                                                mv $Read ${Read%.*}.fastq
                                                        fi   
                                                fi
                                        done
					
                                        # Build directory of intermediate output.
                                        mkdir -p $intermediate_output/$sample_name

                                        # Build three output paths: dropTag, aligment, dropEst.
                                        mkdir -p -p $matrix_path/$sample_name/${sample_name}_01_dropTag $matrix_path/$sample_name/${sample_name}_02_alignment $matrix_path/$sample_name/${sample_name}_03_dropEst
                                        
                                        if [[ "$1" == "Drop-seq" ]] || [[ "$1" == "inDrop_v1" ]] || [[ "$1" == "inDrop_v2" ]];then
                                                # For Drop-seq, inDrop v1 and v2, combine more than one fastq file into one.
                                                cat $Sample/*_1.fastq > $intermediate_output/$sample_name/${sample_name}_1.fastq
                                                cat $Sample/*_2.fastq > $intermediate_output/$sample_name/${sample_name}_2.fastq

                                                for file in `ls -1 $intermediate_output/$sample_name/${sample_name}*.fastq`;do head -100 $file > ${file}.head;done

                                                #2.changeNames
                                                read1=`sed -n 2p $intermediate_output/$sample_name/${sample_name}_1.fastq.head | wc -m`
                                                read2=`sed -n 2p $intermediate_output/$sample_name/${sample_name}_2.fastq.head | wc -m`

                                                echo "parameter_1 is $read1"
                                                if [[ $read1 -gt 30 ]]; then
                                                        mv $intermediate_output/$sample_name/${sample_name}_1.fastq $intermediate_output/$sample_name/${sample_name}_gene.fastq
                                                else
                                                        mv $intermediate_output/$sample_name/${sample_name}_1.fastq $intermediate_output/$sample_name/${sample_name}_barcode.fastq
                                                fi

                                                echo "parameter_2 is $read2"
                                                if [[ $read2 -gt 30 ]]; then
                                                        mv $intermediate_output/$sample_name/${sample_name}_2.fastq $intermediate_output/$sample_name/${sample_name}_gene.fastq
                                                else
                                                        mv $intermediate_output/$sample_name/${sample_name}_2.fastq $intermediate_output/$sample_name/${sample_name}_barcode.fastq
                                                fi

                                                # dropTag
                                                cd $matrix_path/$sample_name/${sample_name}_01_dropTag

                                                if [[ "$1" == "Drop-seq" ]];then
                                                        droptag \
                                                        -c $config_path/config/drop_seq.xml -l $sample_name -p $6 -s -S -r 0 \
                                                        $intermediate_output/$sample_name/${sample_name}_barcode.fastq $intermediate_output/$sample_name/${sample_name}_gene.fastq
                                                elif [[ "$1" == "inDrop_v1" ]] || [[ "$1" == "inDrop_v2" ]];then
                                                        droptag \
                                                        -c $config_path/config/indrop_v1_2.xml -l $sample_name -p $6 -s -S -r 0 \
                                                        $intermediate_output/$sample_name/${sample_name}_barcode.fastq $intermediate_output/$sample_name/${sample_name}_gene.fastq
                                                else
                                                        echo "Check the library strategy of sequencing!"
                                                fi
                                        elif [[ "$1" == "inDrop_v3" ]];then
                                                # For inDrop v3, combine more than one fastq file into one.
                                                cat $Sample/*_1.fastq > $intermediate_output/$sample_name/${sample_name}_gene.fastq
                                                cat $Sample/*_2.fastq > $intermediate_output/$sample_name/${sample_name}_library_tags.fastq
                                                cat $Sample/*_3.fastq > $intermediate_output/$sample_name/${sample_name}_barcode1.fastq
                                                cat $Sample/*_4.fastq > $intermediate_output/$sample_name/${sample_name}_barcode2.fastq

                                                # dropTag
                                                cd $matrix_path/$sample_name/${sample_name}_01_dropTag

                                                droptag \
                                                -c $config_path/config/indrop_v3.xml -l $sample_name -p $6 -s -S -r 0 \
                                                $intermediate_output/$sample_name/${sample_name}_barcode1.fastq $intermediate_output/$sample_name/${sample_name}_barcode2.fastq $intermediate_output/$sample_name/${sample_name}_gene.fastq

                                        else
                                                echo "Check the library strategy of sequencing!"
                                        fi

                                        # alignment
                                        cd $matrix_path/$sample_name/${sample_name}_02_alignment

                                        nohup STAR \
                                        --genomeDir $7 --runThreadN $8 \
                                        --outFileNamePrefix $matrix_path/$sample_name/${sample_name}_02_alignment/${sample_name}. \
                                        --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
                                        --readFilesIn $matrix_path/$sample_name/${sample_name}_01_dropTag/${sample_name}_gene.fastq.tagged.fastq.gz

                                        # dropEst
                                        cd $matrix_path/$sample_name/${sample_name}_03_dropEst

                                        if [[ "$1" == "Drop-seq" ]];then
                                                dropest -m -L eEBA \
                                                -g $9 -c $config_path/config/drop_seq.xml \
                                                -r $matrix_path/$sample_name/${sample_name}_01_dropTag/${sample_name}_gene.fastq.tagged.params.gz \
                                                -o $matrix_path/$sample_name/${sample_name}_03_dropEst/${sample_name}_cell.counts.rds \
                                                $matrix_path/$sample_name/${sample_name}_02_alignment/${sample_name}.Aligned.sortedByCoord.out.bam
                                        elif [[ "$1" == "inDrop_v1" ]] || [[ "$1" == "inDrop_v2" ]];then
                                                dropest -m -L eEBA \
                                                -g $9 -c $config_path/config/indrop_v1_2.xml \
                                                -r $matrix_path/$sample_name/${sample_name}_01_dropTag/${sample_name}_gene.fastq.tagged.params.gz \
                                                -o $matrix_path/$sample_name/${sample_name}_03_dropEst/${sample_name}_cell.counts.rds \
                                                $matrix_path/$sample_name/${sample_name}_02_alignment/${sample_name}.Aligned.sortedByCoord.out.bam
                                        elif [[ "$1" == "inDrop_v3" ]];then
                                                dropest -m -L eEBA \
                                                -g $9 -c $config_path/config/indrop_v3.xml \
                                                -r $matrix_path/$sample_name/${sample_name}_01_dropTag/${sample_name}_gene.fastq.tagged.params.gz \
                                                -o $matrix_path/$sample_name/${sample_name}_03_dropEst/${sample_name}_cell.counts.rds \
                                                $matrix_path/$sample_name/${sample_name}_02_alignment/${sample_name}.Aligned.sortedByCoord.out.bam
                                        fi

                                        # dropReport
                                        dropReport.Rsc \
                                        -t $matrix_path/$sample_name/${sample_name}_01_dropTag/${sample_name}_gene.fastq.tagged.rds \
                                        -m ${10} \
                                        -o $matrix_path/$sample_name/${sample_name}_03_dropEst/${sample_name}_report.html \
                                        $matrix_path/$sample_name/${sample_name}_03_dropEst/${sample_name}_cell.counts.rds
                                fi
                        fi
                done
        done
elif [[ "$2" == "All_samples" ]];then

        for Sample in $5/* 
        do
		sample_name=$(echo "$Sample"| awk -F "/" '{print $NF}')
                if [ -d "$Sample" ];then 
                        for Read in $Sample/*
                                do
                                        if [ -f "$Read" ] && [ "$4" = "sra" ];then
                                                if [ "${Read##*.}"x != "fastq"x ] || [ "${Read##*.}"x != "fq"x ];then
                                                        read_file=$(echo "$Read"| awk -F "/" '{print $NF}')
                                                        echo "The read file $read_file exist."

                                                        # Convert sra file to fastq file.
                                                        fastq-dump --split-files -O $Sample $Read
                                                fi
                                        elif [ -f "$Read" ] && [ "$4" = "fastq" ];then
                                                if [ "${Read##*.}"x = "fq"x ];then
                                                        mv $Read ${Read%.*}.fastq
                                                fi   
                                        fi
                                done

                        # Build directory of intermediate output.
                        mkdir -p $intermediate_output/$sample_name

                        # Build three output paths: dropTag, aligment, dropEst.
                        mkdir -p -p $matrix_path/$sample_name/${sample_name}_01_dropTag $matrix_path/$sample_name/${sample_name}_02_alignment $matrix_path/$sample_name/${sample_name}_03_dropEst

                        if [[ "$1" == "Drop-seq" ]] || [[ "$1" == "inDrop_v1" ]] || [[ "$1" == "inDrop_v2" ]];then
                                # For Drop-seq, inDrop v1 and v2, combine more than one fastq file into one.
                                cat $Sample/*_1.fastq > $intermediate_output/$sample_name/${sample_name}_1.fastq
                                cat $Sample/*_2.fastq > $intermediate_output/$sample_name/${sample_name}_2.fastq

                                for file in `ls -1 $intermediate_output/$sample_name/${sample_name}*.fastq`;do head -100 $file > ${file}.head;done

                                #2.changeNames
                                read1=`sed -n 2p $intermediate_output/$sample_name/${sample_name}_1.fastq.head | wc -m`
                                read2=`sed -n 2p $intermediate_output/$sample_name/${sample_name}_2.fastq.head | wc -m`

                                echo "parameter_1 is $read1"
                                if [[ $read1 -gt 30 ]]; then
                                        mv $intermediate_output/$sample_name/${sample_name}_1.fastq $intermediate_output/$sample_name/${sample_name}_gene.fastq
                                else
                                        mv $intermediate_output/$sample_name/${sample_name}_1.fastq $intermediate_output/$sample_name/${sample_name}_barcode.fastq
                                fi

                                echo "parameter_2 is $read2"
                                if [[ $read2 -gt 30 ]]; then
                                        mv $intermediate_output/$sample_name/${sample_name}_2.fastq $intermediate_output/$sample_name/${sample_name}_gene.fastq
                                else
                                        mv $intermediate_output/$sample_name/${sample_name}_2.fastq $intermediate_output/$sample_name/${sample_name}_barcode.fastq
                                fi

                                # dropTag
                                cd $matrix_path/$sample_name/${sample_name}_01_dropTag

                                if [[ "$1" == "Drop-seq" ]];then
                                        droptag \
                                        -c $config_path/config/drop_seq.xml -l $sample_name -p $6 -s -S -r 0 \
                                        $intermediate_output/$sample_name/${sample_name}_barcode.fastq $intermediate_output/$sample_name/${sample_name}_gene.fastq
                                elif [[ "$1" == "inDrop_v1" ]] || [[ "$1" == "inDrop_v2" ]];then
                                        droptag \
                                        -c $config_path/config/indrop_v1_2.xml -l $sample_name -p $6 -s -S -r 0 \
                                        $intermediate_output/$sample_name/${sample_name}_barcode.fastq $intermediate_output/$sample_name/${sample_name}_gene.fastq
                                else
                                        echo "Check the library strategy of sequencing!"
                                fi
                                        
                        elif [[ "$1" == "inDrop_v3" ]];then
                                # For inDrop v3, combine more than one fastq file into one.
                                cat $Sample/*_1.fastq > $intermediate_output/$sample_name/${sample_name}_gene.fastq
                                cat $Sample/*_2.fastq > $intermediate_output/$sample_name/${sample_name}_library_tags.fastq
                                cat $Sample/*_3.fastq > $intermediate_output/$sample_name/${sample_name}_barcode1.fastq
                                cat $Sample/*_4.fastq > $intermediate_output/$sample_name/${sample_name}_barcode2.fastq

                                # dropTag
                                cd $matrix_path/$sample_name/${sample_name}_01_dropTag

                                droptag \
                                -c $config_path/config/indrop_v3.xml -l $sample_name -p $6 -s -S -r 0 \
                                $intermediate_output/$sample_name/${sample_name}_barcode1.fastq $intermediate_output/$sample_name/${sample_name}_barcode2.fastq $intermediate_output/$sample_name/${sample_name}_gene.fastq

                        else
                                echo "Check the library strategy of sequencing!"
                        fi

                        # alignment
                        cd $matrix_path/$sample_name/${sample_name}_02_alignment

                        nohup STAR \
                        --genomeDir $7 --runThreadN $8 \
                        --outFileNamePrefix $matrix_path/$sample_name/${sample_name}_02_alignment/${sample_name}. \
                        --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
                        --readFilesIn $matrix_path/$sample_name/${sample_name}_01_dropTag/${sample_name}_gene.fastq.tagged.fastq.gz

                        # dropEst
                        cd $matrix_path/$sample_name/${sample_name}_03_dropEst

                        if [[ "$1" == "Drop-seq" ]];then
                                dropest -m -L eEBA \
                                -g $9 -c $config_path/config/drop_seq.xml \
                                -r $matrix_path/$sample_name/${sample_name}_01_dropTag/${sample_name}_gene.fastq.tagged.params.gz \
                                -o $matrix_path/$sample_name/${sample_name}_03_dropEst/${sample_name}_cell.counts.rds \
                                $matrix_path/$sample_name/${sample_name}_02_alignment/${sample_name}.Aligned.sortedByCoord.out.bam
                        elif [[ "$1" == "inDrop_v1" ]] || [[ "$1" == "inDrop_v2" ]];then
                                dropest -m -L eEBA \
                                -g $9 -c $config_path/config/indrop_v1_2.xml \
                                -r $matrix_path/$sample_name/${sample_name}_01_dropTag/${sample_name}_gene.fastq.tagged.params.gz \
                                -o $matrix_path/$sample_name/${sample_name}_03_dropEst/${sample_name}_cell.counts.rds \
                                $matrix_path/$sample_name/${sample_name}_02_alignment/${sample_name}.Aligned.sortedByCoord.out.bam
                        elif [[ "$1" == "inDrop_v3" ]];then
                                dropest -m -L eEBA \
                                -g $9 -c $config_path/config/indrop_v3.xml \
                                -r $matrix_path/$sample_name/${sample_name}_01_dropTag/${sample_name}_gene.fastq.tagged.params.gz \
                                -o $matrix_path/$sample_name/${sample_name}_03_dropEst/${sample_name}_cell.counts.rds \
                                $matrix_path/$sample_name/${sample_name}_02_alignment/${sample_name}.Aligned.sortedByCoord.out.bam
                        fi

                        # dropReport
                        dropReport.Rsc \
                        -t $matrix_path/$sample_name/${sample_name}_01_dropTag/${sample_name}_gene.fastq.tagged.rds \
                        -m ${10} \
                        -o $matrix_path/$sample_name/${sample_name}_03_dropEst/${sample_name}_report.html \
                        $matrix_path/$sample_name/${sample_name}_03_dropEst/${sample_name}_cell.counts.rds
                fi
        done
fi
