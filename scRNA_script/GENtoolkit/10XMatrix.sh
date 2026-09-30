#!/bin/bash
set -eo pipefail

#############################################################################################################################
# $1: [ Designated_samples | All_samples ]
# $2: read_file_type [ sra | fastq ]
# $3: raw_data_path
# $4: cellranger index
# $5: cellranger localcores
# $6: cellranger localmem
# $7: Selected sample list [ sample1-sample2-sample3 ]
# $8: cellranger create-bam [ true | false ]
#############################################################################################################################

single_project_path=$(dirname "$3")
intermediate_output=${single_project_path}/2_output
matrix_path=${single_project_path}/3_expression_result
mkdir -p "$intermediate_output" "$matrix_path"

create_bam=${8:-true}

if [[ "$1" == "All_samples" ]]; then
    for Sample in "$3"/*; do
        [[ -d "$Sample" ]] || continue

        sample_name=$(basename "$Sample")
        echo "The sample $sample_name exists."

        mkdir -p "$intermediate_output/$sample_name"

        if [[ "$2" == "sra" ]]; then
            for Read in "$Sample"/*.sra; do
                [[ -f "$Read" ]] || continue
                run_file=$(basename "$Read")
                echo "The read file $run_file exists."

                fasterq-dump \
                    --split-files \
                    --include-technical \
                    -e "$5" \
                    -O "$Sample" \
                    -t "$Sample" \
                    "$Read"
            done

            for fastq_file in "$Sample"/*.fastq; do
                [[ -f "$fastq_file" ]] || continue

                sra_file=$(basename "$fastq_file")
                # 把根据barcode 推断的逻辑改成直接根据文件名判断 _1.fastq -> R1、_2.fastq -> R2
                # fasterq-dump --split-files normally creates:
                # SRRxxxx_1.fastq = R1
                # SRRxxxx_2.fastq = R2
                # SRRxxxx_3.fastq or more = technical/index reads, optional for cellranger
                if [[ "$sra_file" == *_1.fastq ]]; then
                    mv "$fastq_file" "$Sample/${sample_name}_S1_L001_R1_001.fastq"
                elif [[ "$sra_file" == *_2.fastq ]]; then
                    mv "$fastq_file" "$Sample/${sample_name}_S1_L001_R2_001.fastq"
                elif [[ "$sra_file" =~ _[3-9][0-9]*\.fastq$ ]]; then
                    echo "Skip extra FASTQ file: $sra_file"
                    mkdir -p "$Sample/extra_fastq"
                    mv "$fastq_file" "$Sample/extra_fastq/$sra_file"
                else
                    echo "Skip unrecognized FASTQ file: $sra_file"
                fi
            done

        elif [[ "$2" == "fastq" ]]; then
            for Read in "$Sample"/*.fq; do
                [[ -f "$Read" ]] || continue
                mv "$Read" "${Read%.fq}.fastq"
            done

            echo "Make sure FASTQ names look like:"
            echo "${sample_name}_S1_L001_R1_001.fastq"
            echo "${sample_name}_S1_L001_R2_001.fastq"

        else
            echo "Please check --ReadType/-rt: should be sra or fastq"
            exit 1
        fi

        r1_files=("$Sample"/*_S1_L001_R1_001.fastq)
        r2_files=("$Sample"/*_S1_L001_R2_001.fastq)

        if [[ -e "${r1_files[0]}" ]]; then
            cat "${r1_files[@]}" > "$intermediate_output/$sample_name/${sample_name}_S1_L001_R1_001.fastq"
        else
            echo "No R1 fastq files found for sample: $sample_name"
            continue
        fi

        if [[ -e "${r2_files[0]}" ]]; then
            cat "${r2_files[@]}" > "$intermediate_output/$sample_name/${sample_name}_S1_L001_R2_001.fastq"
        else
            echo "No R2 fastq files found for sample: $sample_name"
            continue
        fi

        i1_files=("$Sample"/*_S1_L001_I1_001.fastq)
        if [[ -e "${i1_files[0]}" ]]; then
            rm -f "${i1_files[@]}"
        fi

        cd "$matrix_path"

        cellranger count \
            --id="$sample_name" \
            --create-bam="$create_bam" \
            --sample="$sample_name" \
            --transcriptome="$4" \
            --fastqs="$intermediate_output/$sample_name" \
            --localcores="$5" \
            --localmem="$6"

    done

elif [[ "$1" == "Designated_samples" ]]; then
    cd "$matrix_path"

    IFS='-' read -r -a samples <<< "$7"

    for each_sample in "${samples[@]}"; do
        cellranger count \
            --id="$each_sample" \
            --create-bam="$create_bam" \
            --sample="$each_sample" \
            --transcriptome="$4" \
            --fastqs="$intermediate_output/$each_sample" \
            --localcores="$5" \
            --localmem="$6"
    done
else
    echo "Please check parameter 1: should be All_samples or Designated_samples"
    exit 1
fi