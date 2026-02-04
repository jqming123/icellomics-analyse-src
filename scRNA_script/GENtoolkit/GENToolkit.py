#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File  : GENToolkit.py
@Author: MING CHEN & Zhu TT
Contact: chenm@big.ac.cn
'''

import os
import argparse
import datetime

"""
    Note: The all directories and files of the whole generating matrix process as shown below:

    Project
        |--- Project Reference (For example, 01.reference)
                                |--- Species name (For example, Homo_sapiens)
                                                |--- reference genome fasta file (.fasta | .fa | .fna file)
                                                |--- reference genome annotation file 1 (.gtf file)
                                                |--- reference genome annotation file 2 (.bed file)
                                                |--- Hisat2 index
                                                                |--- hisat2 index output files, such as genome.1.ht2
                                                |--- RSEM index
                                                              |--- RSEM index output files
                                                |--- Cellranger index 
                                                                    |--- Species_name.genome
                                                                                           |--- fasta
                                                                                           |--- genes
                                                                                           |--- pickle
                                                                                           |--- star
        |--- Bulk/Samrt-seq2 Project 
                                    |--- Single Project Name (For example, GEND000254)
                                                            |--- Raw data (For example, 01.raw)
                                                                        |--- Sample name (For example, GENS00035294, GENS00035295, etc)
                                                                                        |--- sra files or fastq files
                                                            |--- Intermediate output (For example, 02.output)
                                                                                    |--- Sample name (For example, GENS00035294, GENS00035295, etc)
                                                                                                |--- Quality Control (QC) reports, bam files, etc.
                                                            |--- Matrix output (For example, 03.expression)
                                                                            |--- Project_GeneMat_FPKM.txt
                                                                            |--- Project_GeneMat_rawCounts.txt
                                                                            |--- Project_GeneMat_TPM.txt
                                                                            |--- Project_TransMat_FPKM.txt
                                                                            |--- Project_TransMat_rawCounts.txt
                                                                            |--- Project_TransMat_TPM.txt
                                                            |--- Visualization and further analysis
        |--- 10X Project
                        |--- Single Project Name (For example, GEND000326)
                                                |--- Raw data (For example, 01.raw)
                                                            |--- Sample name (For example, GENS00037886)
                                                                            |--- sra files or fastq files
                                                            |--- Intermediate output (For example, 02.output)
                                                                                    |--- Sample name (For example, GENS00037886)
                                                                                                |--- Sample_S1_L001_R2_001.fastq
                                                                                                |--- Sample_S1_L001_R1_001.fastq
                                                                                                |--- Sample_S1_L001_I1_001.fastq
                                                            |--- Matrix output (For example, 03.expression)
                                                                            |--- outs
                                                                            |--- other output files
                                                            |--- Visualization and further analysis
        |--- Drop-seq / inDrop
                    |--- Single Project Name (For example, GEND000158)
                                           |--- Raw data (For example, 01.raw)
                                                       |--- Sample name (For example, GENS00015351, GENS00015352, etc)
                                                                      |--- sra files or fastq files
                                                       |--- Intermediate output (For example, 02.output)
                                                                            |--- Sample name (For example, GENS00015351, GENS00015352, etc)
                                                                                           |--- Processed sra or fastq file
                                                       |--- Matrix output (For example, 03.expression)
                                                                        |--- Sample_01_dropTag
                                                                        |--- Sample_02_alignment
                                                                        |--- Samplpe_03_dropEst
"""

class BulkSingleCellReferenceGenomeIndex():
    """ Before you start generating expressions, you need to build reference genome index first."""

    def __init__(self, project_path, hisat2_thread_num, RSEM_thread_num, reference_genome_fasta, reference_genome_gtf, star_path):
        # The project_path means "Species name"
        self.project_path = project_path
        self.hisat2_thread_num = hisat2_thread_num
        self.RSEM_thread_num = RSEM_thread_num
        self.reference_genome_fasta = reference_genome_fasta
        self.reference_genome_gtf = reference_genome_gtf
        self.star_path = star_path

    def Hisat2Index(self):
        """
        Build bulk hisat2 refernece genome index (bulk).
        """
        os.chdir(self.project_path)
        species_name = os.getcwd().split('/')[-1]
        hisat2_output_path = os.getcwd() + '/' + species_name + '_hisat2'

        os.system('mkdir %s' %(hisat2_output_path))
        os.system('hisat2-build -p %s %s %s/%s' %(self.hisat2_thread_num, self.reference_genome_fasta, hisat2_output_path, "genome"))

    def RSEMIndex(self):
        """
        Build bulk RSEM reference genome index (bulk, Smart-seq2, Drop-seq and inDrop).
        """
        os.chdir(self.project_path)
        species_name = os.getcwd().split('/')[-1]
        RSEM_output_path = os.getcwd() + '/' + species_name + '_RSEM'
        
        os.system('mkdir %s' %(RSEM_output_path))
        os.system('rsem-prepare-reference --gtf %s -p %s --star --star-path %s %s %s/%s' \
            %(self.reference_genome_gtf, self.RSEM_thread_num, self.star_path, self.reference_genome_fasta, RSEM_output_path, species_name))

    def cellrangerIndex(self):
        """
        Build cellranger index reference genome (10X).
        """
        os.chdir(self.project_path)
        species_name = os.getcwd().split('/')[-1]
        cellranger_output_path = species_name + '.genome'
        os.system('cellranger mkref --genome=%s --gene=%s --fasta=%s' %(cellranger_output_path, self.reference_genome_gtf, self.reference_genome_fasta))

class BulkSmartSeqMatrix():
    """ Generate Bulk or Smart-seq2 matrixes."""

    def __init__(self, bulk_smart_seq, designated_all, sample_list, sequencing_type, read_type, raw_data_path, hisat2_index, RSEM_index, bed_file, fasterq_dump_thread_num, fastp_q, fastp_u, fastp_l, fastp_W, fastp_M, fastp_w, hisat2_p, samtools_thread_num, rsem_thread_num, star_path, generate_matrix):
        self.bulk_smart_seq = bulk_smart_seq
        self.designated_all = designated_all
        self.sample_list = sample_list
        self.sequencing_type = sequencing_type
        self.read_type = read_type
        self.raw_data_path = raw_data_path
        self.hisat2_index = hisat2_index
        self.RSEM_index = RSEM_index
        self.bed_file = bed_file
        self.fasterq_dump_thread_num = fasterq_dump_thread_num
        self.fastp_q = fastp_q
        self.fastp_u = fastp_u
        self.fastp_l = fastp_l
        self.fastp_W = fastp_W
        self.fastp_M = fastp_M
        self.fastp_w = fastp_w
        self.hisat2_p = hisat2_p
        self.samtools_thread_num = samtools_thread_num
        self.rsem_thread_num = rsem_thread_num
        self.star_path = star_path
        self.generate_matrix = generate_matrix

    def BulkSmartSeq(self):
        """
        Run the script on the local computer (independent node).
        """ 
        os.system('bash ./BulkSmartSeqMatrix.sh %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s' \
            %(self.bulk_smart_seq, self.designated_all, self.sample_list, self.sequencing_type, self.read_type, self.raw_data_path, \
                self.hisat2_index, self.RSEM_index, self.bed_file, self.fasterq_dump_thread_num, self.fastp_q, \
                    self.fastp_u, self.fastp_l, self.fastp_W, self.fastp_M, self.fastp_w, self.hisat2_p, self.samtools_thread_num, self.rsem_thread_num, self.star_path, self.generate_matrix))

class TenXMatrix():
    """
    Generate 10X matrixes.
    """
    def __init__(self, designated_all, read_type, raw_data_path, cellranger_index, cellranger_localcores, cellranger_localmem, sample_list):
        self.designated_all = designated_all
        self.read_type = read_type
        self.raw_data_path = raw_data_path
        self.cellranger_index = cellranger_index
        self.cellranger_localcores = cellranger_localcores
        self.cellranger_localmem = cellranger_localmem
        self.sample_list = sample_list

    def TenX(self):
        """
        Run the script on the local computer (independent node).
        """
        os.system('bash ./10XMatrix.sh %s %s %s %s %s %s %s' \
            %(self.designated_all, self.read_type, self.raw_data_path, self.cellranger_index, self.cellranger_localcores, self.cellranger_localcores, self.sample_list))

class DropSeqinDropMatrix():
    """ Generate Drop-seq or inDrop matrix."""
    def __init__(self, DropinDrop, designated_all, sample_list, read_type, raw_data_path, dropTag_p, star_index, star_runThreadN, dropEst_g, dropReport_m):
        self.DropinDrop = DropinDrop
        self.designated_all = designated_all
        self.sample_list = sample_list
        self.read_type = read_type
        self.raw_data_path = raw_data_path
        self.dropTag_p = dropTag_p
        self.star_index = star_index
        self.star_runThreadN = star_runThreadN
        self.dropEst_g = dropEst_g
        self.dropReport_m = dropReport_m
        

    def DropSeqinDrop(self):
        """
        Run the script on the local computer (independent node).
        """
        os.system('bash ./DropSeqinDropMatrix.sh %s %s %s %s %s %s %s %s %s %s' \
            %(self.DropinDrop, self.designated_all, self.sample_list, self.read_type, self.raw_data_path, self.dropTag_p, self.star_index, \
                self.star_runThreadN, self.dropEst_g, self.dropReport_m))

def main():
    """ Pass parameters and enter the matrix generating process."""

    # Set parameters
    parser = argparse.ArgumentParser(description = "Please input the right path and choose suitable parameters.")

    # The parameters of reference genome index building (bulk and single cell).
    parser.add_argument('--IndexBuild', '-ib', type = str, default = "index_build", help = "Building reference genome index depends on whether the index exists. [index_build | index_exist ]")
    parser.add_argument('--BuildLibraryType', '-blt', type = str, help = "Library building type, [ Bulk | 10X | Smart-seq2 | inDrop_v1 | inDrop_v2 | inDrop_v3 | Drop-seq ]")
    parser.add_argument('--IndexProjectPath', '-ipp', type = str, help = "The absolute path of the index-building project (Species name).")
    parser.add_argument('--Hisat2ThreadNum', '-htn', type = int, default = 12, help = "The thread number in hisat2 index-building.")
    parser.add_argument('--RSEMThreadNum', '-rtn', type = int, default = 12, help = "The thread number in RSEM index-building.")
    parser.add_argument('--ReferenceGenomeFasta', '-rgf', type = str, help = "The path (absolute path is OK) of the fasta file of reference genome.")
    parser.add_argument('--ReferenceGenomeGtf', '-rgg', type = str, help = "The path (absolute path is OK) of the gtf file of reference genome.")
    parser.add_argument('--StarPath', '-sp', help = "The absolute path of the software STAR.")

    # Partial parameters of bulk and smart-seq2, part of 10X's, Drop-seq's and inDrops.
    parser.add_argument('--DesignatedAll', '-da', type = str, default = "All_samples", help = "The sample list you want to run at once. [ Designated_samples | All_samples ]")
    parser.add_argument('--SampleList', '-sl', type = str, help = "The list samples you designated. [ Sample1-Sample2-Sample3 ] Note: The short dashes are needed.")
    parser.add_argument('--SequencingType', '-st', type = str, help  = "[ pair | single ]")
    parser.add_argument('--ReadType', '-rt', default = "fastq", type = str, help = "[ sra | fastq ]")
    parser.add_argument('--RawData', '-rd', type = str, help = "The absolute path of raw data.")
    parser.add_argument('--Hisat2Index', '-hi', type = str, help = "The absolute path of Hisat2 index, for example, ../Homo_Sapiens_hisat2/genome")
    parser.add_argument('--RSEMIndex', '-ri', type = str, help = "The absolute path of RSEM index, for example, ../Homo_Sapiens_RSEM/Homo_Spaiens")
    parser.add_argument('--BedFile', '-bf', type = str, help = "The absolute path of annotation file (.bed file).")
    parser.add_argument('--FasterqDumpThread', '-fdt', type = int, default = 12, help = "The thread number of fasterq-dump.")
    parser.add_argument('--Fastp_q', '-fq', type = int, default = 20, help = "The parameter \"-q\" in fastp.")
    parser.add_argument('--Fastp_u', '-fu', type = int, default = 20, help = "The parameter \"-u\" in fastp.")
    parser.add_argument('--Fastp_l', '-fl', type = int, default = 40, help = "The parameter \"-l\" in fastp.")
    parser.add_argument('--Fastp_W', '-fW', type = int, default = 4, help = "The parameter \"-W\" in fastp.")
    parser.add_argument('--Fastp_M', '-fM', type = int, default = 20, help = "The parameter \"-M\" in fastp.")
    parser.add_argument('--Fastp_w', '-fw', type = int, default = 12, help = "The parameter \"-w\" in fastp.")
    parser.add_argument('--Hisat2_p', '-hp', type = int, default = 12, help = "The thread number of hisat2.")
    parser.add_argument('--SamtoolsThread', '-std', type = int, default = 12, help = "The thread number of samtools.")
    parser.add_argument('--RSEMThread', '-rtd', type = int, default = 12, help = "The thread number of RSEM.")
    parser.add_argument('--GenerateMatrix', '-gm', type = str, default = "n_matrix", help = "Whether generating matrixes when \"Designated_samples\". [ g_matrix | n_matrix]")

    # Partial parameters of 10X.
    parser.add_argument('--CellrangerIndex', '-ci', type = str, help = "The absolute path of cellranger index, for example, ../Homo_Sapiens/Homo_Sapiens.genome")
    parser.add_argument('--CellrangerLocalCores', '-clc', type = int, default = 12, help = "Caution! Caution! Caution! The default localcores may not be appropriate all the time, you can adjust the localcores according to https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/count.")
    parser.add_argument('--CellrangerLocalMem', '-clm', type = int, default = 64, help = "Caution! Caution! Caution! The default localmem may not be enough all the time, you can adjust the localmem according to https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/count.")

    # Partial parameters of Drop-seq or inDrop.
    parser.add_argument('--DropTag_p', '-dp', type = int, default = 12, help = "The thread number of dropTag process.")
    parser.add_argument('--StarIndex', '-si', type = str, help = "The absolute path of STAR index.")
    parser.add_argument('--StarRunThreadN', '-srtn', type = int, default = 12, help = "The thread number of STAR alignment.")
    parser.add_argument('--DropReport_m', '-dm', type = str, help = "The reference organelle gene's rds file.")
   
    args = parser.parse_args()

    # The parameters of reference genome index building (bulk and single cell).
    indexBuild = args.IndexBuild
    buildLibraryType = args.BuildLibraryType
    indexProjectPath = args.IndexProjectPath
    hisat2ThreadNum = args.Hisat2ThreadNum
    rsemThreadNum = args.RSEMThreadNum
    referenceGenomeFasta = args.ReferenceGenomeFasta
    referenceGenomeGtf = args.ReferenceGenomeGtf
    starPath = args.StarPath

    # The parameters of bulk and smart-seq2, and part of 10X's.
    designatedAll = args.DesignatedAll
    sampleList = args.SampleList
    sequencingType = args.SequencingType
    readType = args.ReadType
    rawData = args.RawData
    hisat2Index = args.Hisat2Index
    rsemIndex = args.RSEMIndex
    bedFile = args.BedFile
    fasterqDumpThread = args.FasterqDumpThread
    FASTP_q = args.Fastp_q
    FASTP_u = args.Fastp_u
    FASTP_l = args.Fastp_l
    FASTP_W = args.Fastp_W
    FASTP_M = args.Fastp_M
    FASTP_w = args.Fastp_w
    HISAT2_p = args.Hisat2_p
    samtoolsThread = args.SamtoolsThread
    rsemThread = args.RSEMThread
    generateMatrix = args.GenerateMatrix

    # Partial parameters of 10X.
    cellrangerIndex = args.CellrangerIndex
    cellrangerLocalCores = args.CellrangerLocalCores
    cellrangerLocalMem = args.CellrangerLocalMem

    # Partial parameter of Drop-seq, or inDrop v1, v2, v3.
    DROPTAG_p = args.DropTag_p
    starIndex = args.StarIndex
    starRunThreadN = args.StarRunThreadN
    dropReport_m = args.DropReport_m

    # Pass parameters. -- Upstream analysis
    bulkSingleCellReferenceGenomeIndex = BulkSingleCellReferenceGenomeIndex(indexProjectPath, hisat2ThreadNum, rsemThreadNum, referenceGenomeFasta, referenceGenomeGtf, starPath)
    bulkSmartSeqMatrix = BulkSmartSeqMatrix(buildLibraryType, designatedAll, sampleList, sequencingType, readType, rawData, hisat2Index, rsemIndex, bedFile, fasterqDumpThread, FASTP_q, FASTP_u, FASTP_l, FASTP_W, FASTP_M, FASTP_w, HISAT2_p, samtoolsThread, rsemThread, starPath, generateMatrix)
    tenXMatrix = TenXMatrix(designatedAll, readType, rawData, cellrangerIndex, cellrangerLocalCores, cellrangerLocalMem, sampleList)
    dropSeqinDropMatrix = DropSeqinDropMatrix(buildLibraryType, designatedAll, sampleList, readType, rawData, DROPTAG_p, starIndex, starRunThreadN, referenceGenomeGtf, dropReport_m)

    # Record the start time of the index-building process.
    start_time = datetime.datetime.now()
    print("\n")
    print("******************************************************")
    print("Start time: " + str(start_time) + ".")
    print("******************************************************")
    print("\n")

    # Bulk RNA-seq
    if buildLibraryType == "Bulk":
        script_path = os.getcwd()
        if indexBuild == "index_build":
            bulkSingleCellReferenceGenomeIndex.Hisat2Index()
            bulkSingleCellReferenceGenomeIndex.RSEMIndex()
        os.chdir(script_path)
        bulkSmartSeqMatrix.BulkSmartSeq()

    # Smart-seq2 scRNA-seq
    if buildLibraryType == "Smart-seq2":
        script_path = os.getcwd()
        if indexBuild == "index_build":
            bulkSingleCellReferenceGenomeIndex.Hisat2Index()
            bulkSingleCellReferenceGenomeIndex.RSEMIndex()
        os.chdir(script_path)
        bulkSmartSeqMatrix.BulkSmartSeq()

    # 10X Genomics scRNA-seq
    if buildLibraryType == "10X":
        script_path = os.getcwd()
        if indexBuild == "index_build":
            bulkSingleCellReferenceGenomeIndex.cellrangerIndex()
        os.chdir(script_path)
        tenXMatrix.TenX()

    # Drop-seq or inDrop(v1, v2, v3) scRNA-seq
    if buildLibraryType == "Drop-seq" or buildLibraryType == "inDrop_v1" or buildLibraryType == "inDrop_v2" or buildLibraryType == "inDrop_v3":
        script_path = os.getcwd()
        if indexBuild == "index_build":
            bulkSingleCellReferenceGenomeIndex.RSEMIndex()
        os.chdir(script_path)
        dropSeqinDropMatrix.DropSeqinDrop()
            
    # Record the end time of the index-building process.
    end_time = datetime.datetime.now()
    total_time = end_time - start_time
    print("\n")
    print("****************************************************")
    print("End time: " + str(end_time) + ".")
    print("Total time:" + str(total_time) + ".")
    print("****************************************************")
    print("\n")

if __name__ == '__main__':
    main()
