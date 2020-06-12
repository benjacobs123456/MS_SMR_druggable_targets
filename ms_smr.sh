#! /bin/bash
#$ -pe smp 4
#$ -l h_vmem=4G
#$ -l h_rt=1:0:0
#$ -cwd
#$ -j y
#$ -t 1:22
#$ -o /data/scratch/hmy117



cd /data/Wolfson-UKBB-Dobson/SMR


./smr_Linux --bfile ../1kg_reference/filtered_chr${SGE_TASK_ID} \
--gwas-summary ./gwas/ms_gwas.ma \
--beqtl-summary /data/Wolfson-UKBB-Dobson/SMR/datasets/eqtlgen/cis-eQTLs-full_eQTLGen_AF_incl_nr_formatted_20191212.new.txt_besd-dense \
--out ms_smr_chr${SGE_TASK_ID} \
--thread-num $NSLOTS \
--diff-freq-prop 1

cat ms_smr_chr${SGE_TASK_ID}\.smr >> overall_ms_smr_results.txt
rm ms_smr_chr${SGE_TASK_ID}\.smr*


