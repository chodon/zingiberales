awk '$3>0' CURC/CURC.cov | awk 'BEGIN{FS="\t"}{print $1"\t"$2-1"\t"$2}' | mergeBed > sep_exons_full_cds/curc_covered.bed
awk '$3>0' COST/COST_full_cds_novo540/COST.cov | awk 'BEGIN{FS="\t"}{print $1"\t"$2-1"\t"$2}' | mergeBed > sep_exons_full_cds/cost_covered.bed
awk '$3>0' KNKV/KNKV_full_cds_novo540/KNKV.cov | awk 'BEGIN{FS="\t"}{print $1"\t"$2-1"\t"$2}' | mergeBed > sep_exons_full_cds/knkv_covered.bed
awk '$3>0' TZNS/TZNS_full_cds_novo540/TZNS.cov | awk 'BEGIN{FS="\t"}{print $1"\t"$2-1"\t"$2}' | mergeBed > sep_exons_full_cds/tzns_covered.bed
awk '$3>0' JNUB/JNUB_full_cds_novo540/JNUB.cov | awk 'BEGIN{FS="\t"}{print $1"\t"$2-1"\t"$2}' | mergeBed > sep_exons_full_cds/jnub_covered.bed
awk '$3>0' LSKK/LSKK_full_cds_novo540/LSKK.cov | awk 'BEGIN{FS="\t"}{print $1"\t"$2-1"\t"$2}' | mergeBed > sep_exons_full_cds/lskk_covered.bed
awk '$3>0' UOEL/UOEL_full_cds_novo540/UOEL.cov | awk 'BEGIN{FS="\t"}{print $1"\t"$2-1"\t"$2}' | mergeBed > sep_exons_full_cds/uoel_covered.bed
awk '$3>0' BDJQ/BDJQ_full_cds_novo540/BDJQ.cov | awk 'BEGIN{FS="\t"}{print $1"\t"$2-1"\t"$2}' | mergeBed > sep_exons_full_cds/bdjq_covered.bed
exit
