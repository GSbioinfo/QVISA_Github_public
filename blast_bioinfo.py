#from Bio import SeqIO
import numpy as np
import subprocess, sys
import pipes
import pandas as pd
#from joblib import Parallel, delayed
#import multiprocessing
import sqlite3
import datetime as dt
import itertools
from cruzdb import Genome

#def fun_Bioinfo(alignDict):
#----------------------------------------------------------------------------------
#Bioinformatics section of the VIS analysis
#----------------------------------------------------------------------------------


#take in visAllTSVPath
fpath=sys.argv[1]
final_vis_data=pd.read_csv(fpath, sep='\t')

finalResultsDirectory=sys.argv[1].replace(".tsv","")
genomeLabel=sys.argv[2]
genome = Genome(db=genomeLabel)

coln=[str(cc) for cc in final_vis_data.columns]
#final_vis_data=final_vis_data[coln[:-1]].copy()
site_table= final_vis_data[["CLONE_ID","CHROMO","STRAND","SITE"]].copy()
site_table["TEMPSITE"]=list(itertools.repeat(0.0,len(site_table.index)))
site_table.loc[site_table["STRAND"]=='+',"TEMPSITE"]= site_table.loc[site_table["STRAND"]=='+',"SITE"]+5.0
site_table.loc[site_table["STRAND"]=='-',"TEMPSITE"]= site_table.loc[site_table["STRAND"]=='-',"SITE"]-5.0
site_table["SITE1"]=site_table[["SITE","TEMPSITE"]].min(axis=1).astype(int)
site_table["SITE2"]=site_table[["SITE","TEMPSITE"]].max(axis=1).astype(int)
site_table["DOT1"]=list(itertools.repeat('.',len(site_table.index)))
site_table["DOT2"]=list(itertools.repeat('.',len(site_table.index)))
site_table.sort_values(['CHROMO', 'SITE1'], ascending=[True, True],inplace=True)
site_table["KEY"]=site_table.CHROMO.map(str)+'|'+site_table.SITE1.map(str)+'|'+site_table.SITE2.map(str)

vis_bioinfo_all=pd.DataFrame(columns=["CLONE_ID","CHR","GENE_ID","GENE_NAME","GENE_STRAND","DISTANCE","FEATURE"])
vis_bioinfo_closest=pd.DataFrame(columns=["CLONE_ID","CHR","GENE_ID","GENE_NAME","GENE_STRAND","DISTANCE","FEATURE"])
ensemid_list=[]
for ro in site_table.iterrows():
    loc1=ro[1].tolist()
    #loc1=site_table.ix[0].tolist()
    tupgens=genome.knearest("refGene",str(loc1[1]),loc1[5],loc1[6],k=10)
    temp=pd.DataFrame()
    gen_name_list=[]
    for ygen in tupgens:
        if ygen.name2 in gen_name_list:
            continue
        gen_name_list.append(ygen.name2)
        geneid=ygen.name
        genesym=ygen.name2
        dist,feat=ygen.distance(loc1[5],loc1[6],True)
        temp=temp.append(pd.Series([loc1[0],ygen.chrom,geneid,genesym,ygen.strand,dist,feat]).to_frame().transpose())
    temp.columns=["CLONE_ID","CHR","GENE_ID","GENE_NAME","GENE_STRAND","DISTANCE","FEATURE"]
    temp1=temp.loc[temp["DISTANCE"]==min(temp.DISTANCE),].copy()
    vis_bioinfo_all=vis_bioinfo_all.append(temp)
    vis_bioinfo_closest=vis_bioinfo_closest.append(temp1)
vis_dict=vis_bioinfo_closest.set_index("CLONE_ID").to_dict('index')
for clone in vis_dict.keys():
    vis_dict[clone]['FEATURE']=vis_bioinfo_closest.loc[vis_bioinfo_closest["CLONE_ID"]==clone,"FEATURE"].tolist()
    vis_dict[clone]['GENE_ID']=vis_bioinfo_closest.loc[vis_bioinfo_closest["CLONE_ID"]==clone,"GENE_ID"].tolist()
    vis_dict[clone]['GENE_NAME']=vis_bioinfo_closest.loc[vis_bioinfo_closest["CLONE_ID"]==clone,"GENE_NAME"].tolist()
    vis_dict[clone]['GENE_STRAND']=vis_bioinfo_closest.loc[vis_bioinfo_closest["CLONE_ID"]==clone,"GENE_STRAND"].tolist()
vis_bioinfo_closest=pd.DataFrame(vis_dict).transpose()
vis_dict=vis_bioinfo_all.set_index("CLONE_ID").to_dict('index')
for clone in vis_dict.keys():
    vis_dict[clone]['DISTANCE']=vis_bioinfo_all.loc[vis_bioinfo_all["CLONE_ID"]==clone,"DISTANCE"].tolist()
    vis_dict[clone]['FEATURE']=vis_bioinfo_all.loc[vis_bioinfo_all["CLONE_ID"]==clone,"FEATURE"].tolist()
    vis_dict[clone]['GENE_ID']=vis_bioinfo_all.loc[vis_bioinfo_all["CLONE_ID"]==clone,"GENE_ID"].tolist()
    vis_dict[clone]['GENE_NAME']=vis_bioinfo_all.loc[vis_bioinfo_all["CLONE_ID"]==clone,"GENE_NAME"].tolist()
    vis_dict[clone]['GENE_STRAND']=vis_bioinfo_all.loc[vis_bioinfo_all["CLONE_ID"]==clone,"GENE_STRAND"].tolist() 
vis_bioinfo_all=pd.DataFrame(vis_dict).transpose()

vis_bioinfo_all["CLONE_ID"]=vis_bioinfo_all.index.astype(str).tolist()
vis_bioinfo_closest["CLONE_ID"]=vis_bioinfo_closest.index.astype(str).tolist()
final_vis_bio_all=final_vis_data.merge(vis_bioinfo_all,on='CLONE_ID',how='outer')
final_vis_bio_closest=final_vis_data.merge(vis_bioinfo_closest,on='CLONE_ID',how='outer')

visBioAllTsvPath=finalResultsDirectory+"_bioinformatics_all_genes.tsv"
xlwriterbioinfo = pd.ExcelWriter(finalResultsDirectory+"_bioinformatics_genes.xlsx")
final_vis_bio_all.to_csv(visBioAllTsvPath,sep="\t")
visBioClosestTsvPath=finalResultsDirectory+"_bioinformatics_closest_genes.tsv"
final_vis_bio_closest.to_csv(visBioClosestTsvPath,sep="\t")
final_vis_bio_closest.to_excel(xlwriterbioinfo,sheet_name='closest_genes',index=False)
final_vis_bio_all.to_excel(xlwriterbioinfo,sheet_name='all_genes',index=False)
xlwriterbioinfo.save()
