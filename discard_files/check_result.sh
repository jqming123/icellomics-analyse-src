# 统计目录/文件数，并检查 rsem.genes.results
expected_dirs=8
expected_files=25
# cd /gpfs/zhaowm_group/gaoxiaojing/CellLine/projects_results/rnseq
cd /hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/transcriptome_projects/rna_count_result



for d in SRR15559{102..121}; do
  if [ ! -d "$d" ]; then
    echo "$d : MISSING DIR"
    continue
  fi
  n_dirs=$(find "$d" -type d | wc -l)
  n_files=$(find "$d" -type f | wc -l)

  # 检查关键文件是否存在
  if [ -f "$d/rsem/${d}_rsem.genes.results" ]; then
    key="rsem.genes.results OK"
  else
    key="MISSING rsem.genes.results"
  fi

  # 打印并标注是否达标
  if [ "$n_dirs" -eq "$expected_dirs" ] && [ "$n_files" -eq "$expected_files" ]; then
    echo "$d : $n_dirs directories, $n_files files | $key | PASS"
  else
    echo "$d : $n_dirs directories, $n_files files | $key | FAIL (expected $expected_dirs dirs, $expected_files files)"
  fi
done

