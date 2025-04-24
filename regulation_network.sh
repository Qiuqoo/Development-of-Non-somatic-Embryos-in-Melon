bedtools flank -i genes.bed -g genome_file.txt -l 2000 -r 0 -s > promoters.bed
cat promoters.bed genes.bed > gene_regions.bed

input_dir="/data3/users/Qiuqoo/CM_ZJU/hintatac/output/motif_output"

gene_bed="/data3/users/Qiuqoo/CM_ZJU/hintatac/genome/gene_regions.bed"

output_dir="/data3/users/Qiuqoo/CM_ZJU/hintatac/GRN/footprint_gene_mapping/"

bedtools intersect -wa -wb \
        -a PTE_cleaned.bed \
        -b "$gene_bed" \
        > "${output_dir}/${prefix}_motif_gene_overlap.bed"
		
######################################################################################################		
#-                                              PYTHON
######################################################################################################		

# -*- coding: utf-8 -*-
import pandas as pd
import glob
import os

input_dir = "/data3/users/Qiuqoo/CM_ZJU/hintatac/GRN/footprint_gene_mapping"  # 存放DAIx_motif_gene.bed的目录
output_file = "/data3/users/Qiuqoo/CM_ZJU/hintatac/GRN/footprint_gene_mapping/motif_gene_mapping_all.tsv"

rows = []

for file in glob.glob(f"{input_dir}/*_motif_gene_overlap.bed"):
    timepoint = os.path.basename(file).split("_")[0]  
    df = pd.read_csv(file, sep="\t", header=None)
    

    for _, row in df.iterrows():
        motif_id = row[3]
        gene_id = row[9]  
        rows.append([timepoint, motif_id, gene_id])


output_df = pd.DataFrame(rows, columns=["Timepoint", "Motif", "Gene"])
output_df.drop_duplicates().to_csv(output_file, sep="\t", index=False)
######################################################################################################		
#-                                              PYTHON
######################################################################################################		

# -*- coding: utf-8 -*-
import pandas as pd


motif_tf_df = pd.read_csv("/data3/users/Qiuqoo/CM_ZJU/hintatac/genome/motif-TF.csv")

motif_tf_df["motifid_clean"] = motif_tf_df["motifid"].str.replace(r"_\d+\.\d+$", "", regex=True)


mapping_df = pd.read_csv("/data3/users/Qiuqoo/CM_ZJU/hintatac/GRN/footprint_gene_mapping/motif_gene_mapping_all.tsv", sep="\t")


mapping_df["Motif_clean"] = mapping_df["Motif"].str.replace(r"_\d+\.\d+$", "", regex=True)


merged_df = mapping_df.merge(motif_tf_df, left_on="Motif_clean", right_on="motifid_clean", how="inner")


network_df = merged_df[["geneid", "Gene", "Timepoint", "Motif"]].drop_duplicates()
network_df.columns = ["TF", "TargetGene", "Timepoint", "Motif"]


network_df.to_csv("/data3/users/Qiuqoo/CM_ZJU/hintatac/GRN/TF_Target_Network.tsv", sep="\t", index=False)
