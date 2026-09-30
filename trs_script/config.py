import os

# -----------------------------------------------------------------
# 配置区域: 所有可调整的参数都定义在这里
# -----------------------------------------------------------------

REF_BASE_DIR = "/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/resources/ref_genome"

CONFIG = {
    # 1. 路径设置
    "paths": {
        "sra_data_root": '/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/transcriptome_projects/rna_rawdata',
        "project_results_root": '/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/transcriptome_projects/rna_count_result',
        "main_script_path": '/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/resources/src/trs_script/bulkRNA-seq_E4_v2.sh',
        "custom_bin_path": '/hpcdisk1/zhaowm_group/gaoxiaojing/softwares/miniforge3/bin:/hpcdisk1/zhaowm_group/gaoxiaojing/softwares/pixi_0.59.0/trs_env/.pixi/envs/default/bin'
    },

    # 2. SLURM作业调度系统设置
    "slurm_settings": {
        "partition": 'corexd192',
        "time": '7-00:00:00',  # D-HH:MM:SS 格式 
        "log_dir": '/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/transcriptome_projects/logs',
    },

    # 3. 工具和资源参数
    "tool_params": {
        "threads": 4,
        "memory_gb": 60
    },

    # 4. 参考基因组设置
    "reference_genomes": {
        # 该基因组已弃用
#        "CriGri-PICRH-1.0": { 
#            "star_index": os.path.join(REF_BASE_DIR,"CriGri-PICRH-1.0","star.index"),
#            "kallisto_gene_idx": os.path.join(REF_BASE_DIR,"CriGri-PICRH-1.0","kallisto.index","GCF_003668045.3_CriGri-PICRH-1.0_genomic.gene.fa.idx"),
#            "kallisto_transcript_idx": os.path.join(REF_BASE_DIR,"CriGri-PICRH-1.0","kallisto.index","GCF_003668045.3_CriGri-PICRH-1.0_genomic.transcript.fa.idx"),
#            "rsem_ref_prefix": os.path.join(REF_BASE_DIR,"CriGri-PICRH-1.0","rsem.index","reference")
#        },
        "CH_Ensembl": { 
            "star_index": os.path.join(REF_BASE_DIR,"CriGri-PICRH-1.0_Ensembl","star.index"),
            "kallisto_gene_idx": os.path.join(REF_BASE_DIR,"CriGri-PICRH-1.0_Ensembl","kallisto.index","CriGri-PICRH-1.0.115.gene.idx"),
            "kallisto_transcript_idx": os.path.join(REF_BASE_DIR,"CriGri-PICRH-1.0_Ensembl","kallisto.index","CriGri-PICRH-1.0.115.transcript.idx"),
            "rsem_ref_prefix": os.path.join(REF_BASE_DIR,"CriGri-PICRH-1.0_Ensembl","rsem.index","reference")
        },
        "hg38_Ensembl": { 
            "star_index": os.path.join(REF_BASE_DIR,"hg38_Ensembl","star.index"),
            "kallisto_gene_idx": os.path.join(REF_BASE_DIR,"hg38_Ensembl","kallisto.index","Homo_sapiens.GRCh38.dna_sm.primary_assembly.gene.fa.idx"),
            "kallisto_transcript_idx": os.path.join(REF_BASE_DIR,"hg38_Ensembl","kallisto.index","Homo_sapiens.GRCh38.dna_sm.primary_assembly.transcript.fa.idx"),
            "rsem_ref_prefix": os.path.join(REF_BASE_DIR,"hg38_Ensembl","rsem.index","reference")
        },
        "Cattle_E_ARSUCD2": { 
            "star_index": os.path.join(REF_BASE_DIR, "Cattle_E_ARSUCD2", "star.index"),
            "kallisto_gene_idx": os.path.join(REF_BASE_DIR, "Cattle_E_ARSUCD2", "kallisto.index", "ARS-UCD2.0.115.gene.idx"),
            "kallisto_transcript_idx": os.path.join(REF_BASE_DIR, "Cattle_E_ARSUCD2", "kallisto.index", "ARS-UCD2.0.115.transcript.idx"),
            "rsem_ref_prefix": os.path.join(REF_BASE_DIR, "Cattle_E_ARSUCD2", "rsem.index", "reference")
        },
        "Chicken_E_GRCg7b": { 
            "star_index": os.path.join(REF_BASE_DIR, "Chicken_E_GRCg7b", "star.index"),
            "kallisto_gene_idx": os.path.join(REF_BASE_DIR, "Chicken_E_GRCg7b", "kallisto.index", "bGalGal1.mat.broiler.GRCg7b.115.gene.idx"),
            "kallisto_transcript_idx": os.path.join(REF_BASE_DIR, "Chicken_E_GRCg7b", "kallisto.index", "bGalGal1.mat.broiler.GRCg7b.115.transcript.idx"),
            "rsem_ref_prefix": os.path.join(REF_BASE_DIR, "Chicken_E_GRCg7b", "rsem.index", "reference")
        },
        "Dog_E_UUGSD": { 
            "star_index": os.path.join(REF_BASE_DIR, "Dog_E_UUGSD", "star.index"),
            "kallisto_gene_idx": os.path.join(REF_BASE_DIR, "Dog_E_UUGSD", "kallisto.index", "UU_Cfam_GSD_1.0.115.gene.idx"),
            "kallisto_transcript_idx": os.path.join(REF_BASE_DIR, "Dog_E_UUGSD", "kallisto.index", "UU_Cfam_GSD_1.0.115.transcript.idx"),
            "rsem_ref_prefix": os.path.join(REF_BASE_DIR, "Dog_E_UUGSD", "rsem.index", "reference")
        },
        "GreenMonkey_E_ChlSab1.1": { 
            "star_index": os.path.join(REF_BASE_DIR, "GreenMonkey_E_ChlSab1.1", "star.index"),
            "kallisto_gene_idx": os.path.join(REF_BASE_DIR, "GreenMonkey_E_ChlSab1.1", "kallisto.index", "ChlSab1.1.115.gene.idx"),
            "kallisto_transcript_idx": os.path.join(REF_BASE_DIR, "GreenMonkey_E_ChlSab1.1", "kallisto.index", "ChlSab1.1.115.transcript.idx"),
            "rsem_ref_prefix": os.path.join(REF_BASE_DIR, "GreenMonkey_E_ChlSab1.1", "rsem.index", "reference")
        },
        "Mouse_E_GRCm39": { 
            "star_index": os.path.join(REF_BASE_DIR, "Mouse_E_GRCm39", "star.index"),
            "kallisto_gene_idx": os.path.join(REF_BASE_DIR, "Mouse_E_GRCm39", "kallisto.index", "GRCm39.115.gene.idx"),
            "kallisto_transcript_idx": os.path.join(REF_BASE_DIR, "Mouse_E_GRCm39", "kallisto.index", "GRCm39.115.transcript.idx"),
            "rsem_ref_prefix": os.path.join(REF_BASE_DIR, "Mouse_E_GRCm39", "rsem.index", "reference")
        },
        "Pig_E_Sscrofa11.1": { 
            "star_index": os.path.join(REF_BASE_DIR, "Pig_E_Sscrofa11.1", "star.index"),
            "kallisto_gene_idx": os.path.join(REF_BASE_DIR, "Pig_E_Sscrofa11.1", "kallisto.index", "Sscrofa11.1.115.gene.idx"),
            "kallisto_transcript_idx": os.path.join(REF_BASE_DIR, "Pig_E_Sscrofa11.1", "kallisto.index", "Sscrofa11.1.115.transcript.idx"),
            "rsem_ref_prefix": os.path.join(REF_BASE_DIR, "Pig_E_Sscrofa11.1", "rsem.index", "reference")
        },
    }
}
