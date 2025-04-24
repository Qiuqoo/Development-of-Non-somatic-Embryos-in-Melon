###########################################################################################################################
#                                                              Step1
###########################################################################################################################
#!/bin/bash

input_dir="/data3/users/Qiuqoo/CM_ZJU/WGBS/CX_report"
promoter_bed="/data3/users/Qiuqoo/CM_ZJU/hintatac/genome/gene_promoter.bed"
genebody_bed="/data3/users/Qiuqoo/CM_ZJU/hintatac/genome/gene_genebody.bed"

mkdir -p "${input_dir}/meth_output"
cd "$input_dir"

process_sample() {
  file="$1"
  context="$2"
  name=$(basename "$file" | cut -d'_' -f1)
  echo "Processing $name [$context]..."

  awk -v type=$context 'BEGIN{OFS="\t"} $6==type {
    meth=$4; unmeth=$5; tot=meth+unmeth;
    ratio=(tot > 0 ? meth / tot : 0);
    print $1, $2-1, $2, ratio, tot, meth;
  }' "$file" | sort -k1,1 -k2,2n | bgzip -c > ${name}.${context}.bed.gz


  tabix -p bed ${name}.${context}.bed.gz

  bedtools map -a "$promoter_bed" -b ${name}.${context}.bed.gz -c 4,5,6 -o mean,sum,sum > "${name}_promoter_meth_${context}.txt"
  bedtools map -a "$genebody_bed" -b ${name}.${context}.bed.gz -c 4,5,6 -o mean,sum,sum > "${name}_genebody_meth_${context}.txt"
}

export -f process_sample
export promoter_bed genebody_bed

parallel --jobs 7 process_sample {} {1} ::: "$input_dir"/*.deduplicated.CX_report.txt ::: CG CHG CHH
###########################################################################################################################
#                                                              Step2
###########################################################################################################################
#!/bin/bash

input_dir="/data3/users/Qiuqoo/CM_ZJU/WGBS/CX_report"
cd "$input_dir"

contexts=(CG CHG CHH)
regions=(promoter genebody)

for region in "${regions[@]}"; do
  for ctx in "${contexts[@]}"; do
    echo "Merging: $region | $ctx"

    files=(*_${region}_meth_${ctx}.txt)

    sample_names=()
    for f in "${files[@]}"; do
      sample=$(basename "$f" | cut -d'_' -f1)
      sample_names+=("$sample")
    done

    cut -f1-6 "${files[0]}" > merged_${region}_${ctx}.txt


    for i in "${!files[@]}"; do
      cut -f7 "${files[i]}" > temp_${sample_names[$i]}.meth
      paste merged_${region}_${ctx}.txt temp_${sample_names[$i]}.meth > temp && mv temp merged_${region}_${ctx}.txt
    done

    header="chr\tstart\tend\tgene_id\tregion\tstrand"
    for s in "${sample_names[@]}"; do
      header="${header}\t${s}"
    done
    sed -i "1s/.*/$header/" merged_${region}_${ctx}.txt

    rm temp_*.meth
  done
done
###########################################################################################################################
#                                                              Step3
###########################################################################################################################
#!/bin/bash

input="merged_promoter_CG.txt" 
output="merged_promoter_CG_avg.txt"


head -n 1 "$input" | awk 'BEGIN{OFS="\t"} {
  printf "%s\t%s\t%s\t%s\t%s\t%s", $1, $2, $3, $4, $5, $6;
  print "\tGE\tHE\tLCE\tLTE\tNSE\tPCE\tPYE";
}' > "$output"


tail -n +2 "$input" | awk 'BEGIN{OFS="\t"} {
  printf "%s\t%s\t%s\t%s\t%s\t%s", $1, $2, $3, $4, $5, $6;

  # 样品组：每组3列，依次是 GE1-3、HE1-3...
  GE=($7 + $8 + $9) / 3;
  HE=($10 + $11 + $12) / 3;
  LCE=($13 + $14 + $15) / 3;
  LTE=($16 + $17 + $18) / 3;
  NSE=($19 + $20 + $21) / 3;
  PCE=($22 + $23 + $24) / 3;
  PTE=($25 + $26 + $27) / 3;

  printf "\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n", GE, HE, LCE, LTE, NSE, PCE, PTE;
}' >> "$output"
###########################################################################################################################
#                                                              Step3
###########################################################################################################################
#!/bin/bash


gene_file="merged_genebody_CHH_avg.txt"
prom_file="merged_promoter_CHH_avg.txt"
output="merged_gene_promoter_CHH_avg.txt"

head -n 1 "$gene_file" | awk 'BEGIN{OFS="\t"} {
  printf "%s\t%s\t%s\t%s\t%s", $1, $2, $3, $4, $6;
  print "\tGE_avg\tHE_avg\tLCE_avg\tLTE_avg\tNSE_avg\tPCE_avg\tPYE_avg\tGE_sum\tHE_sum\tLCE_sum\tLTE_sum\tNSE_sum\tPCE_sum\tPYE_sum";
}' > "$output"


tail -n +2 "$gene_file" | while read -r g_line; do
  gene_id=$(echo "$g_line" | cut -f4)
  p_line=$(grep -P "\t$gene_id\t" "$prom_file")

  if [[ -n "$p_line" ]]; then

    gb_vals=$(echo "$g_line" | cut -f7-)
    pr_vals=$(echo "$p_line" | cut -f7-)


    echo -e "$gb_vals\n$pr_vals" | awk -v meta="$g_line" 'BEGIN{OFS="\t"} {
      for (i=1; i<=NF; i++) {
        sum[i] += $i;
      }
    }
    END {
      split(meta, a, "\t");
      printf "%s\t%s\t%s\t%s\t%s", a[1], a[2], a[3], a[4], a[6];
      for (i=1; i<=NF; i++) printf "\t%.6f", sum[i]/2;
      for (i=1; i<=NF; i++) printf "\t%.6f", sum[i];
      printf "\n";
    }' >> "$output"
  fi
done
















