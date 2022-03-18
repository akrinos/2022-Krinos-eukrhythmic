import os
import pandas as pd
import numpy as np
import plotnine
from datetime import date
import scipy
import seaborn as sns
from Bio import SeqIO
import matplotlib.pyplot as plt

output_dir="/vortexfs1/omics/alexander/akrinos/2021-testing-eukrhythmic/eukrhythmic_paper_tara"
salmon_dir="/vortexfs1/omics/alexander/akrinos/2021-testing-eukrhythmic/eukrhythmic_paper_tara/intermediate-files/04-compare/xx-individual-mapping/salmon"
salmon_dir="/vortexfs1/omics/alexander/akrinos/2021-testing-eukrhythmic/eukrhythmic_paper_tara/intermediate-files/04-compare/xx-CDS-mapping/salmon"

## Read in BLAST searches against MATOU database

matou_fold=os.path.join("/vortexfs1/omics/alexander","akrinos/2022-Krinos-eukrhythmic/MATOU/BLAST_MATOU",
                        "blast_out_matou/top_1")
output_dir="/vortexfs1/omics/alexander/akrinos/2021-testing-eukrhythmic/eukrhythmic_paper_tara"
salmon_dir="/vortexfs1/omics/alexander/akrinos/2021-testing-eukrhythmic/eukrhythmic_paper_tara/intermediate-files/04-compare/09-CAG-mapping/salmon"
salmon_dir="/vortexfs1/omics/alexander/akrinos/2021-testing-eukrhythmic/eukrhythmic_paper_tara/intermediate-files/04-compare/xx-CDS-mapping/salmon"
all_blast_file=pd.DataFrame()
for dir_curr in os.listdir(matou_fold):
    #if not os.path.isfile(os.path.join(matou_fold,dir_curr,"combined","blast_formatted.out")):
    #    continue
    
    curr_sample = dir_curr.split("_CAG")[0]
    subdir = "EUKulele_merged_CAG_cds_14Mar22" #"EUKulele_merged"
    
    eggnog_CAG = pd.read_csv(os.path.join(output_dir, "eggnog_CAG", 
                       curr_sample + ".emapper.annotations"),
                                  sep="\t",comment="#",
                            header=None,names=["query","seed_ortholog","evalue","score",
                                               "eggNOG_OGs","max_annot_lvl","COG_category",
                                               "Description","Preferred_name","GOs","EC","KEGG_ko",
                                               "KEGG_Pathway","KEGG_Module","KEGG_Reaction","KEGG_rclass",
                                               "BRITE","KEGG_TC","CAZy","BiGG_Reaction","PFAMs"])
    
    EUKulele_file = pd.read_csv(os.path.join(output_dir,#subdir,
                                             subdir,
                                             "taxonomy_estimation",dir_curr.split("_CAG")[0]+\
                                             "_CAG-estimated-taxonomy.out"),sep="\t")
    
    eggnog_CAG = pd.read_csv(os.path.join(output_dir, "eggnog_CAG", 
                       curr_sample + ".emapper.annotations"),
                                  sep="\t",comment="#",
                            header=None,names=["query","seed_ortholog","evalue","score",
                                               "eggNOG_OGs","max_annot_lvl","COG_category",
                                               "Description","Preferred_name","GOs","EC","KEGG_ko",
                                               "KEGG_Pathway","KEGG_Module","KEGG_Reaction","KEGG_rclass",
                                               "BRITE","KEGG_TC","CAZy","BiGG_Reaction","PFAMs"])
    eggnog_CAG = eggnog_CAG.groupby("query")[["GOs","KEGG_ko"]].\
                    agg(lambda x: "-".join(sorted(list(set(x))))).reset_index(drop=False)
    eggnog_CAG["KEGG_ko"] = [str(curr) if str(curr) != "-" else \
                         "NoAnnot" for curr in eggnog_CAG.KEGG_ko]
    eggnog_CAG["GOs"] = [str(curr) if str(curr) != "-" else \
                         "NoAnnot" for curr in eggnog_CAG.GOs]
    
    EUKulele_file["split_tscpt_name"] = EUKulele_file["transcript_name"]
    
    salmon_file = pd.read_csv(os.path.join(salmon_dir,
                                           dir_curr.split("_CAG")[0] + "_CAG_quant",
                                           "quant.sf"),sep="\t")
    blast_file=pd.read_csv(os.path.join(matou_fold,dir_curr),
                                        #,"combined","blast_formatted.out"),
                                        #header=None,
                         #names=["qseqid","sseqid","pident","length","mismatch",
                         #       "gapopen","qstart","qend","sstart","send","evalue","bitscore"],
                           sep=",")
    blast_file["split_tscpt"] = [curr.split(".p")[0] for curr in blast_file.qseqid]
    blast_file=pd.merge(blast_file,salmon_file,how="right",left_on="qseqid",right_on="Name")
    blast_file=pd.merge(blast_file,EUKulele_file,how="left",left_on="Name",right_on="transcript_name")
    
    eggnog_CAG["split_tscpt_name_egg"] = [curr.split(".p")[0] for curr in eggnog_CAG["query"]]
    merged_kos = pd.merge(eggnog_CAG,blast_file, how="right",left_on="query",right_on="Name").drop("query",axis=1)
    
    merged_kos["KEGG_ko"] = ["NoAnnot" if curr!=curr else curr for curr in merged_kos.KEGG_ko]
    merged_kos["GOs"] = ["NoAnnot" if curr!=curr else curr for curr in merged_kos.GOs]
    merged_kos["Sample"] = dir_curr.split("_CAG")[0]
    
    all_blast_file = pd.concat([all_blast_file,merged_kos])
    
    
all_blast_file["HasMATOUMatch"] = ["Yes" if (curr == curr) else "No"\
                                   for curr in all_blast_file.split_tscpt]
all_blast_file["HasEUKuleleAnnot"] = ["Yes" if (curr == curr) else "No"\
                                   for curr in all_blast_file.classification]
all_blast_file["HasGOAnnot"] = ["Yes" if (curr != "NoAnnot") & (curr == curr) else "No"\
                                for curr in all_blast_file.GOs]
all_blast_file["HasKEGGAnnot"] = ["Yes" if (curr != "NoAnnot") & (curr == curr) else "No"\
                                for curr in all_blast_file.KEGG_ko]
all_blast_file["HasSalmonAnnot"] = "Yes"

all_blast_file = all_blast_file.drop_duplicates(subset=["HasMATOUMatch",
               "HasEUKuleleAnnot","HasGOAnnot","Length","HasSalmonAnnot",
               "HasKEGGAnnot","Sample","Name"])
print(all_blast_file[(all_blast_file.HasGOAnnot=="Yes") & \
                       (all_blast_file.HasEUKuleleAnnot=="Yes") &
                       (all_blast_file.HasMATOUMatch=="No")], flush=True)
all_blast_file.drop_duplicates(subset=["HasMATOUMatch",
               "HasEUKuleleAnnot","HasGOAnnot","HasSalmonAnnot",
               "HasKEGGAnnot","Sample","Name"]).loc[(all_blast_file.Length>=150),["HasMATOUMatch",
               "HasEUKuleleAnnot","HasGOAnnot","HasSalmonAnnot",
               "HasKEGGAnnot","Sample","Name"]].groupby(["HasMATOUMatch",
               "HasEUKuleleAnnot","HasGOAnnot",
               "HasKEGGAnnot","Sample"])["Name"].count().reset_index().to_csv(os.path.join("..","data-output","counted_blast_file_length_cutoff_150bp_revised.csv"))
all_blast_file.loc[(all_blast_file["Length"]>=150),["HasMATOUMatch",
               "HasEUKuleleAnnot","HasGOAnnot","Length","HasSalmonAnnot",
               "HasKEGGAnnot","Sample","Name"]].groupby(["HasMATOUMatch",
               "HasEUKuleleAnnot","HasGOAnnot",
               "HasKEGGAnnot","Sample"])["Length"].mean().reset_index().to_csv(os.path.join("..","data-output","avg_len_blast_file_length_cutoff_150bp.csv"))

length_plot=(plotnine.ggplot(all_blast_file) + plotnine.geom_boxplot(plotnine.aes(x = "HasMATOUMatch",
                                                                     y = "Length")) + plotnine.scale_y_log10()+
                plotnine.theme_bw(base_size=16))
try:
    length_plot.save(filename = "lengths_matou_vs_not.pdf",width=5,height=5,units="in")
except:
    print("no plot")

all_blast_file.loc[(all_blast_file.Length>=150),["HasMATOUMatch",\
               "HasEUKuleleAnnot","HasGOAnnot","HasSalmonAnnot",\
               "HasKEGGAnnot","Sample","Name"]].to_csv(os.path.join("..","data-output","all_blast_file_length_cutoff_150bp.csv"))

matou_fold=os.path.join("/vortexfs1/omics/alexander","akrinos/2022-Krinos-eukrhythmic/MATOU/BLAST_MATOU",
                        "blast_out_matou/top_1")
output_dir="/vortexfs1/omics/alexander/akrinos/2021-testing-eukrhythmic/eukrhythmic_paper_tara"
salmon_dir="/vortexfs1/omics/alexander/akrinos/2021-testing-eukrhythmic/eukrhythmic_paper_tara/intermediate-files/04-compare/09-CAG-mapping/salmon"
salmon_dir="/vortexfs1/omics/alexander/akrinos/2021-testing-eukrhythmic/eukrhythmic_paper_tara/intermediate-files/04-compare/xx-CDS-mapping/salmon"
all_blast_file=pd.DataFrame()
all_salmon_sum=pd.DataFrame()
for dir_curr in os.listdir(matou_fold):
    #if not os.path.isfile(os.path.join(matou_fold,dir_curr,"combined","blast_formatted.out")):
    #    continue
    
    curr_sample = dir_curr.split("_CAG")[0]
    subdir = "EUKulele_merged_CAG_cds_14Mar22" #"EUKulele_merged"
    
    eggnog_CAG = pd.read_csv(os.path.join(output_dir, "eggnog_CAG", 
                       curr_sample + ".emapper.annotations"),
                                  sep="\t",comment="#",
                            header=None,names=["query","seed_ortholog","evalue","score",
                                               "eggNOG_OGs","max_annot_lvl","COG_category",
                                               "Description","Preferred_name","GOs","EC","KEGG_ko",
                                               "KEGG_Pathway","KEGG_Module","KEGG_Reaction","KEGG_rclass",
                                               "BRITE","KEGG_TC","CAZy","BiGG_Reaction","PFAMs"])
    
    EUKulele_file = pd.read_csv(os.path.join(output_dir,#subdir,
                                             subdir,
                                             "taxonomy_estimation",dir_curr.split("_CAG")[0]+\
                                             "_CAG-estimated-taxonomy.out"),sep="\t")
    
    eggnog_CAG = pd.read_csv(os.path.join(output_dir, "eggnog_CAG", 
                       curr_sample + ".emapper.annotations"),
                                  sep="\t",comment="#",
                            header=None,names=["query","seed_ortholog","evalue","score",
                                               "eggNOG_OGs","max_annot_lvl","COG_category",
                                               "Description","Preferred_name","GOs","EC","KEGG_ko",
                                               "KEGG_Pathway","KEGG_Module","KEGG_Reaction","KEGG_rclass",
                                               "BRITE","KEGG_TC","CAZy","BiGG_Reaction","PFAMs"])
    eggnog_CAG = eggnog_CAG.groupby("query")[["GOs","KEGG_ko"]].\
                    agg(lambda x: "-".join(sorted(list(set(x))))).reset_index(drop=False)
    eggnog_CAG["KEGG_ko"] = [str(curr) if str(curr) != "-" else \
                         "NoAnnot" for curr in eggnog_CAG.KEGG_ko]
    eggnog_CAG["GOs"] = [str(curr) if str(curr) != "-" else \
                         "NoAnnot" for curr in eggnog_CAG.GOs]
    
    EUKulele_file["split_tscpt_name"] = EUKulele_file["transcript_name"]
    
    salmon_file = pd.read_csv(os.path.join(salmon_dir,
                                           dir_curr.split("_CAG")[0] + "_CAG_quant",
                                           "quant.sf"),sep="\t")
    blast_file=pd.read_csv(os.path.join(matou_fold,dir_curr),
                                        #,"combined","blast_formatted.out"),
                                        #header=None,
                         #names=["qseqid","sseqid","pident","length","mismatch",
                         #       "gapopen","qstart","qend","sstart","send","evalue","bitscore"],
                           sep=",")
    blast_file["split_tscpt"] = [curr.split(".p")[0] for curr in blast_file.qseqid]
    blast_file=pd.merge(blast_file,salmon_file,how="right",left_on="qseqid",right_on="Name")
    blast_file=pd.merge(blast_file,EUKulele_file,how="left",left_on="Name",right_on="transcript_name")
    
    eggnog_CAG["split_tscpt_name_egg"] = [curr.split(".p")[0] for curr in eggnog_CAG["query"]]
    merged_kos = pd.merge(eggnog_CAG,blast_file, how="right",left_on="query",right_on="Name").drop("query",axis=1)
    
    merged_kos["KEGG_ko"] = ["NoAnnot" if curr!=curr else curr for curr in merged_kos.KEGG_ko]
    merged_kos["GOs"] = ["NoAnnot" if curr!=curr else curr for curr in merged_kos.GOs]
    merged_kos["Sample"] = dir_curr.split("_CAG")[0]

    merged_kos["HasMATOUMatch"] = ["Yes" if (curr == curr) else "No"\
                                   for curr in merged_kos.split_tscpt]
    merged_kos["HasEUKuleleAnnot"] = ["Yes" if (curr == curr) else "No"\
                                       for curr in merged_kos.classification]
    merged_kos["HasGOAnnot"] = ["Yes" if (curr != "NoAnnot") & (curr == curr) else "No"\
                                    for curr in merged_kos.GOs]
    merged_kos["HasKEGGAnnot"] = ["Yes" if (curr != "NoAnnot") & (curr == curr) else "No"\
                                    for curr in merged_kos.KEGG_ko]
    merged_kos["HasSalmonAnnot"] = "Yes"
    
    merged_kos = merged_kos[merged_kos.Length >= 150]
    merged_kos["Sample"] = dir_curr.split("_CAG")[0]
    salmon_sum = merged_kos.groupby(["HasSalmonAnnot","HasKEGGAnnot","Sample",
                                     "HasGOAnnot","HasEUKuleleAnnot",
                                     "HasMATOUMatch"])[["TPM","NumReads"]].sum()
    
    all_salmon_sum = pd.concat([all_salmon_sum,salmon_sum.reset_index()])
    all_blast_file = pd.concat([all_blast_file,merged_kos[merged_kos.HasMATOUMatch == "No"]])
    
all_blast_file.to_csv(os.path.join("..","data-output","eukulele_only.csv"))
all_salmon_sum.to_csv(os.path.join("..","data-output","salmon_summed_reads.csv"))