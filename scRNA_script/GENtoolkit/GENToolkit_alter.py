#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File  : GENToolkit_alter.py
@Author: MING CHEN & Zhu TT (altered)
Contact: chenm@big.ac.cn
'''

import os
import argparse
import datetime
import glob

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DEFAULT_RESOURCES_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, '..', '..', '..'))

"""
        Note: The directory layout for your environment is as follows:

        /hpcdisk1/zhaowm_group/gaoxiaojing/CellLine
                |--- resources
                |       |--- ref_genome
                |       |       |--- <GenomeName> (e.g., CriGri-PICRH-1.0, hg38_Ensemble)
                |       |               |--- genome fasta (.fna/.fa/.fasta[.gz])
                |       |               |--- annotation (.gtf/.gtf.gz, optional .bed)
                |       |               |--- rsem.index (optional)
                |       |               |--- star.index (optional)
                |       |               |--- hisat2_index/ (genome.*.ht2)
                |       |               |--- cellranger_index/ (cellranger mkref output)
                |       |--- src/scRNA_script/GENtoolkit
                |               |--- GENToolkit_alter.py (this script)
                |               |--- BulkSmartSeqMatrix.sh
                |               |--- 10XMatrix.sh
                |               |--- DropSeqinDropMatrix.sh
                |               |--- config/
                |
                |--- scRNA_projects  (by BioProject ID)
                                |--- <BioProjectID>
                                                |--- 1_raw/        # raw sra/fastq
                                                |--- 2_output/     # QC, bam, intermediate
                                                |--- 3_expression_result/ # output matrices
                                                |--- 4_jobs/       # generated slurm scripts
                                                |--- 5_logs/       # job logs

        Adapted for directory structure:
            <resources_root>/ref_genome/<GenomeName>/...
            <resources_root>/src/scRNA_script/GENtoolkit/ (this script and *.sh)
"""

def _first_match(directory, patterns):
    for pattern in patterns:
        matches = glob.glob(os.path.join(directory, pattern))
        if matches:
            matches.sort()
            return matches[0]
    return None


def _infer_hisat2_index_prefix(genome_dir):
    preferred_prefix = os.path.join(genome_dir, 'hisat2_index', 'genome')
    if glob.glob(preferred_prefix + '.1.ht2'):
        return preferred_prefix

    ht2 = glob.glob(os.path.join(genome_dir, '*.1.ht2'))
    if not ht2:
        ht2 = glob.glob(os.path.join(genome_dir, '*/*.1.ht2'))
    if not ht2:
        return None
    ht2.sort()
    first = ht2[0]
    return first.replace('.1.ht2', '')


def _infer_reference_paths(genome_dir):
    reference = {
        'fasta': _first_match(genome_dir, ['*.fna', '*.fa', '*.fasta', '*.fa.gz', '*.fna.gz']),
        'gtf': _first_match(genome_dir, ['*.gtf', '*.gtf.gz']),
        'bed': _first_match(genome_dir, ['*.bed']),
        'rsem': None,
        'star': None,
        'hisat2': None,
        'cellranger': None,
    }

    rsem_candidate = os.path.join(genome_dir, 'rsem.index')
    if os.path.exists(rsem_candidate):
        reference['rsem'] = rsem_candidate

    star_candidate = os.path.join(genome_dir, 'star.index')
    if os.path.exists(star_candidate):
        reference['star'] = star_candidate

    hisat2_prefix = _infer_hisat2_index_prefix(genome_dir)
    if hisat2_prefix:
        reference['hisat2'] = hisat2_prefix

    cellranger_dir = os.path.join(genome_dir, 'cellranger_index')
    if os.path.isdir(cellranger_dir):
        reference['cellranger'] = cellranger_dir
    else:
        cellranger_candidates = glob.glob(os.path.join(genome_dir, '*.genome'))
        if cellranger_candidates:
            cellranger_candidates.sort()
            reference['cellranger'] = cellranger_candidates[0]

    return reference


class BulkSingleCellReferenceGenomeIndex():
    """Before you start generating expressions, you need to build reference genome index first."""

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
        hisat2_output_path = os.path.join(os.getcwd(), 'hisat2_index')
        os.makedirs(hisat2_output_path, exist_ok=True)
        os.system('hisat2-build -p %s %s %s/%s' % (self.hisat2_thread_num, self.reference_genome_fasta, hisat2_output_path, "genome"))

    def RSEMIndex(self):
        """
        Build bulk RSEM reference genome index (bulk, Smart-seq2, Drop-seq and inDrop).
        """
        os.chdir(self.project_path)
        species_name = os.path.basename(os.getcwd())
        RSEM_output_path = os.path.join(os.getcwd(), species_name + '_RSEM')

        os.system('mkdir %s' % (RSEM_output_path))
        os.system('rsem-prepare-reference --gtf %s -p %s --star --star-path %s %s %s/%s' \
            % (self.reference_genome_gtf, self.RSEM_thread_num, self.star_path, self.reference_genome_fasta, RSEM_output_path, species_name))

    def cellrangerIndex(self):
        """
        Build cellranger index reference genome (10X).
        """
        os.chdir(self.project_path)
        cellranger_output_path = 'cellranger_index'
        os.system('cellranger mkref --genome=%s --genes=%s --fasta=%s' % (cellranger_output_path, self.reference_genome_gtf, self.reference_genome_fasta))


class BulkSmartSeqMatrix():
    """Generate Bulk or Smart-seq2 matrixes."""

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
        script_path = os.path.join(SCRIPT_DIR, 'BulkSmartSeqMatrix.sh')
        os.system('bash %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s' \
            % (script_path, self.bulk_smart_seq, self.designated_all, self.sample_list, self.sequencing_type, self.read_type, self.raw_data_path, \
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
        script_path = os.path.join(SCRIPT_DIR, '10XMatrix.sh')
        os.system('bash %s %s %s %s %s %s %s %s' \
            % (script_path, self.designated_all, self.read_type, self.raw_data_path, self.cellranger_index, self.cellranger_localcores, self.cellranger_localcores, self.sample_list))


class DropSeqinDropMatrix():
    """Generate Drop-seq or inDrop matrix."""
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
        script_path = os.path.join(SCRIPT_DIR, 'DropSeqinDropMatrix.sh')
        os.system('bash %s %s %s %s %s %s %s %s %s %s %s' \
            % (script_path, self.DropinDrop, self.designated_all, self.sample_list, self.read_type, self.raw_data_path, self.dropTag_p, self.star_index, \
                self.star_runThreadN, self.dropEst_g, self.dropReport_m))


def main():
    """Pass parameters and enter the matrix generating process."""

    # Set parameters
    parser = argparse.ArgumentParser(description = "Please input the right path and choose suitable parameters.")

    # New parameters for resource layout
    parser.add_argument('--ResourcesRoot', '-rr', type = str, default = DEFAULT_RESOURCES_ROOT, help = "Root path that contains ref_genome/ and src/. Default is auto-detected.")
    parser.add_argument('--GenomeName', '-gn', type = str, help = "Genome folder name under <ResourcesRoot>/ref_genome/.")

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
    parser.add_argument('--Hisat2Index', '-hi', type = str, help = "The absolute path of Hisat2 index, for example, ../hisat2_index/genome")
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
    parser.add_argument('--CellrangerIndex', '-ci', type = str, help = "The absolute path of cellranger index, for example, ../cellranger_index")
    parser.add_argument('--CellrangerLocalCores', '-clc', type = int, default = 12, help = "Caution! Caution! Caution! The default localcores may not be appropriate all the time, you can adjust the localcores according to https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/count.")
    parser.add_argument('--CellrangerLocalMem', '-clm', type = int, default = 64, help = "Caution! Caution! Caution! The default localmem may not be enough all the time, you can adjust the localmem according to https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/count.")

    # Partial parameters of Drop-seq or inDrop.
    parser.add_argument('--DropTag_p', '-dp', type = int, default = 12, help = "The thread number of dropTag process.")
    parser.add_argument('--StarIndex', '-si', type = str, help = "The absolute path of STAR index.")
    parser.add_argument('--StarRunThreadN', '-srtn', type = int, default = 12, help = "The thread number of STAR alignment.")
    parser.add_argument('--DropReport_m', '-dm', type = str, help = "The reference organelle gene's rds file.")

    args = parser.parse_args()

    resources_root = os.path.abspath(args.ResourcesRoot)
    genome_name = args.GenomeName

    # Auto-infer reference paths from ref_genome/<GenomeName>
    if genome_name:
        genome_dir = os.path.join(resources_root, 'ref_genome', genome_name)
        if os.path.isdir(genome_dir):
            inferred = _infer_reference_paths(genome_dir)
            if not args.IndexProjectPath:
                args.IndexProjectPath = genome_dir
            if not args.ReferenceGenomeFasta:
                args.ReferenceGenomeFasta = inferred['fasta']
            if not args.ReferenceGenomeGtf:
                args.ReferenceGenomeGtf = inferred['gtf']
            if not args.BedFile:
                args.BedFile = inferred['bed']
            if not args.Hisat2Index:
                args.Hisat2Index = inferred['hisat2']
            if not args.RSEMIndex:
                args.RSEMIndex = inferred['rsem']
            if not args.StarIndex:
                args.StarIndex = inferred['star']
            if not args.CellrangerIndex:
                args.CellrangerIndex = inferred['cellranger']

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
        script_path = SCRIPT_DIR
        if indexBuild == "index_build":
            bulkSingleCellReferenceGenomeIndex.Hisat2Index()
            bulkSingleCellReferenceGenomeIndex.RSEMIndex()
        os.chdir(script_path)
        bulkSmartSeqMatrix.BulkSmartSeq()

    # Smart-seq2 scRNA-seq
    if buildLibraryType == "Smart-seq2":
        script_path = SCRIPT_DIR
        if indexBuild == "index_build":
            bulkSingleCellReferenceGenomeIndex.Hisat2Index()
            bulkSingleCellReferenceGenomeIndex.RSEMIndex()
        os.chdir(script_path)
        bulkSmartSeqMatrix.BulkSmartSeq()

    # 10X Genomics scRNA-seq
    if buildLibraryType == "10X":
        script_path = SCRIPT_DIR
        if indexBuild == "index_build":
            bulkSingleCellReferenceGenomeIndex.cellrangerIndex()
        os.chdir(script_path)
        tenXMatrix.TenX()

    # Drop-seq or inDrop(v1, v2, v3) scRNA-seq
    if buildLibraryType == "Drop-seq" or buildLibraryType == "inDrop_v1" or buildLibraryType == "inDrop_v2" or buildLibraryType == "inDrop_v3":
        script_path = SCRIPT_DIR
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
