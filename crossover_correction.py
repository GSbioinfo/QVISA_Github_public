#This code is to correct signal crossovers and errors in signature mutations 
import sqlite3
import pandas as pd 
import itertools
from itertools import repeat
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.SeqIO import FastaIO
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
import subprocess, sys, os
import pipes
import datetime as dt
def cross_correct_fun(ani_detail_dic,run_dir):
    #*******Creating temp folder for the BLAST analysis**********
    cctemp=run_dir+'/'+'cc_temp/'
    if not os.path.exists(cctemp):
        os.makedirs(cctemp)
    #*****Creating SQL database************
    cross_crr_db = sqlite3.connect(':memory:')
    vis_cursor=cross_crr_db.cursor()
    
    records1=[]
    records1_shrt=[]
    records1_long=[]
    
    animalID=list(ani_detail_dic.keys())
    
    #print(ani_detail_dic[animalID[0]]["sampleIDList"])
    colnam_animals={}
    #sys.exit(0)
    print("Hello")
    for animal in animalID:
        alignDict_cc=ani_detail_dic[animal]
        record_animal=[]
        final_file=alignDict_cc["finalFile"]
        cc_dir=alignDict_cc['cross_corrDirectory']
        cc_file=alignDict_cc['cross_corrFile']
        vis_raw_data=pd.read_csv(final_file,sep='\t') 
        samples_IDs=ani_detail_dic[animal]["sampleIDList"]
        table_name="TABLE_"+animal
        vis_cursor.execute("DROP TABLE IF EXISTS {tn}".format(tn=table_name))
        colname_list=[]
        colnam_animals[animal]=vis_raw_data.columns.tolist()
        for cln in vis_raw_data.columns.tolist():
            if "_C" not in cln:
                colname_list.append(cln)
        vis_freq= vis_raw_data[colname_list].copy()
        sam_total_list=[]
        sam_count_list=[]
        #*****Determining valid and invalid VIS within a same animal *************
        #*****Here only the conflicting VIS with different signature mutation are resolved**********
        for samp in samples_IDs:
            vis_freq[samp+"_TOTAL_FREQ"]=vis_freq[samp+"_TOTAL"]/sum(vis_freq[samp+"_TOTAL"])
            sam_total_list.append(samp+"_TOTAL_FREQ")
            sam_count_list.append(samp+"_TOTAL")
        vis_freq["MAX_FREQ"]=vis_freq[sam_total_list].max(1)
        vis_freq["COUNT_MAX"]=vis_freq[sam_count_list].max(1)
        vis_freq["COUNT_SAMP"]=vis_freq[sam_count_list].idxmax(1)
        colname_list_freq=[xx.replace("_TOTAL_FREQ","") for xx in vis_freq.columns.tolist()]
        vis_freq.columns=colname_list_freq
        vis_freq["MAX_SAMP"]=vis_freq[samples_IDs].idxmax(1)
        vis_freq["STATUS"]=0*vis_freq["STATUS"]
        vis_freq["NEO_VIS"]=vis_freq["CHROMO"].map(str)+"!"+vis_freq["STRAND"].map(str)+"!"+vis_freq["SITE"].map(str)
        vis_freq["VALIDITY"]=list(itertools.repeat("VALID",len(vis_freq.index)))
        groups=vis_freq.groupby("NEO_VIS")
        list_group=list(groups.groups.keys())
        for glist in list_group:
            temp_group_df=groups.get_group(glist).copy()
            if len(temp_group_df)>1:
                temp_group_df=temp_group_df.sort_values("MAX_FREQ",ascending=False) 
                temp_group_df=temp_group_df.reset_index()
                max_fre = temp_group_df["MAX_FREQ"].max()
                max_cln = temp_group_df.loc[temp_group_df["MAX_FREQ"]==max_fre,"CLONE_ID"].tolist()[0]
                clone_list= temp_group_df.loc[temp_group_df["CLONE_ID"]!=max_cln,"CLONE_ID"].tolist()
                for clid in clone_list:
                    #***********10x cutoff is done here**********
                    if max_fre/temp_group_df.loc[temp_group_df["CLONE_ID"]==clid,"MAX_FREQ"].tolist()[0]>10.0:
                        vis_freq.loc[vis_freq["CLONE_ID"]==clid,"VALIDITY"]="INVALID" 
                    else:
                        if 10*vis_freq.loc[vis_freq["CLONE_ID"]==clid,"TOTAL_ALL"].tolist()[0] < vis_freq.loc[vis_freq["CLONE_ID"]==max_cln,"TOTAL_ALL"].tolist()[0]:
                            vis_freq.loc[vis_freq["CLONE_ID"]==clid,"VALIDITY"]="INVALID"
                        else:
                            vis_freq.loc[vis_freq["CLONE_ID"]==clid,"VALIDITY"]="XVALID" 
                            vis_freq.loc[vis_freq["CLONE_ID"]==max_cln,"VALIDITY"]="XVALID" 
        vis_raw_data=vis_raw_data.merge(vis_freq[["CLONE_ID","NEO_VIS","MAX_FREQ","MAX_SAMP","VALIDITY","COUNT_MAX","COUNT_SAMP"]],on='CLONE_ID',how='outer')
        vis_raw_data.to_sql(table_name, cross_crr_db)

        for visrow in vis_raw_data.loc[vis_raw_data["VALIDITY"]!="INVALID",["NEO_VIS","MAX_FREQ","COUNT_SAMP","TYPE","GSEQ","GLEN","COUNT_MAX"]].itertuples():
            vis_row=list(visrow)[1:]
            seqid=vis_row[0]+"|"+animal+"|"+vis_row[2].replace("_TOTAL","")+'|'+'{:.2e}'.format(vis_row[1])+'|'+'{:d}'.format(int(vis_row[6]))+'|'+vis_row[3]+'|'+str(vis_row[5])
            records1.append(SeqRecord(Seq(vis_row[4],IUPAC.ambiguous_dna),id=seqid,description=""))
            if vis_row[5]<70:
                records1_shrt.append(SeqRecord(Seq(vis_row[4],IUPAC.ambiguous_dna),id=seqid,description=""))
            else:
                records1_long.append(SeqRecord(Seq(vis_row[4],IUPAC.ambiguous_dna),id=seqid,description=""))

    #****Creating VIS seqeunces data for all the animals in one run     
    vis_fasta_path=cctemp+"all_vis.fasta"
    handle_records1=open(vis_fasta_path, 'w')
    fasta_out = FastaIO.FastaWriter(handle_records1, wrap=None)
    fasta_out.write_file(records1)
    handle_records1.close()
    
    vis_short_path=cctemp+"vis_short.fasta"
    handle_records_short=open(vis_short_path, 'w')
    fasta_out_short = FastaIO.FastaWriter(handle_records_short, wrap=None)
    fasta_out_short.write_file(records1_shrt)
    handle_records_short.close()

    vis_long_path=cctemp+"vis_long.fasta"
    handle_records_long=open(vis_long_path, 'w')
    fasta_out_long = FastaIO.FastaWriter(handle_records_long, wrap=None)
    fasta_out_long.write_file(records1_long)
    handle_records_long.close()
    optblast="/opt/latesblast/ncbi-blast-2.7.1+/bin/"
    make_blastdb_cmd=optblast+'./makeblastdb -in '+vis_fasta_path+' -input_type fasta -dbtype nucl -parse_seqids'
    pdb=subprocess.Popen(make_blastdb_cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE, shell=True)
    pdb.communicate()
    
    #*******************Running BLAST *********************
    blast_cmd=optblast+'./blastn -query %s -db %s -perc_identity 95 -max_target_seqs 500 -num_threads 8 -outfmt "6 qseqid sseqid pident evalue qstart qend sstart send qlen" '%(vis_long_path,vis_fasta_path) 
    blast_cmd=blast_cmd+"|awk "+pipes.quote('{OFS="\t"; gsub(/[|]/,"\t",$1);gsub(/[|]/,"\t",$2);if($1 != $2 && $5 == 1 && $7 == 1 && ($6-($5-1))/$9 > 0.9) print $0 }')+'| gzip -c >%sout_long_cross_crrect.txt.gz'%cctemp
    p1=subprocess.Popen(blast_cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE, shell=True)
    p1.communicate()
    blast_cmd=optblast+'./blastn -query %s -db %s -task blastn-short -perc_identity 95 -max_target_seqs 500 -num_threads 8 -outfmt "6 qseqid sseqid pident evalue qstart qend sstart send qlen" '%(vis_short_path,vis_fasta_path) 
    blast_cmd=blast_cmd+"|awk "+pipes.quote('{OFS="\t"; gsub(/[|]/,"\t",$1);gsub(/[|]/,"\t",$2);if($1 != $2 && $5 == 1 && $7 == 1 && ($6-($5-1))/$9 > 0.9) print $0 }')+'| gzip -c >%sout_short_cross_crrect.txt.gz'%cctemp
    p1=subprocess.Popen(blast_cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE, shell=True)
    p1.communicate()
    index_start = 1
    j = 0
    colname_df=["R_NEO_VIS","R_ANIMAL","R_SAMPLEID","R_MAX_FREQ","R_MAX_COUNT","R_TYPE","R_GLEN","Q_NEO_VIS","Q_ANIMAL","Q_SAMPLEID","Q_MAX_FREQ","Q_MAX_COUNT","Q_TYPE","Q_GLEN","PIDENT","EVALUE","R_START","R_END","Q_START","Q_END"]
    blast_df=pd.DataFrame()
    chunksize=20000
    start = dt.datetime.now()
    vis_cursor.execute("DROP TABLE IF EXISTS blast_out")
    #**********Reading BLAST output**************************
    for df in pd.read_table(cctemp+"out_long_cross_crrect.txt.gz",sep='\t', header=None,names=colname_df,compression='gzip',engine='c',chunksize=20000, iterator=True,index_col=False):
        if df.dropna().empty:
            continue
        df.index += index_start
        print ('{} seconds: completed {} rows'.format((dt.datetime.now() - start).seconds, j*chunksize))
        df.to_sql("blast_out", cross_crr_db, if_exists='append')
        
        index_start = df.index[-1] + 1
        j+=1
    for df in pd.read_table(cctemp+"out_short_cross_crrect.txt.gz",sep='\t', header=None,names=colname_df,compression='gzip',engine='c',chunksize=20000, iterator=True,index_col=False):
        if df.dropna().empty:
            continue
        df.index += index_start
        print ('{} seconds: completed {} rows'.format((dt.datetime.now() - start).seconds, j*chunksize))
        df.to_sql("blast_out", cross_crr_db, if_exists='append')
        
        index_start = df.index[-1] + 1
        j+=1
    in_cr=vis_cursor.execute("CREATE INDEX qind1 ON blast_out (R_NEO_VIS)")
    in_cr.fetchone()

    #*********Writing crossover corrected data for each animal***************************
    for tname in animalID:
        alignDict_cc=ani_detail_dic[tname]
        dx = pd.read_sql_query('SELECT * FROM blast_out WHERE R_ANIMAL== "{an}"'.format(an=tname), cross_crr_db)
        tname_sql="TABLE_"+tname
        vis_animal = pd.read_sql_query('SELECT * FROM {tn}'.format(tn=tname_sql), cross_crr_db)
        
        vis_animal["DROP_VIS"]=list(itertools.repeat(False,len(vis_animal.index)))
        vis_animal["MAX_CROSS_COUNT"]=list(itertools.repeat(0.0,len(vis_animal.index)))
        vis_animal["MAX_CROSS_FREQ"]=list(itertools.repeat(0.0,len(vis_animal.index)))
        for nvis in dx["R_NEO_VIS"].unique().tolist():
            tempdx=dx.loc[dx["R_NEO_VIS"]==nvis]
            vis_animal.loc[vis_animal["NEO_VIS"]==nvis,"DROP_VIS"]=tempdx["R_MAX_COUNT"].max()/tempdx["Q_MAX_COUNT"].max()<=10
            vis_animal.loc[vis_animal["NEO_VIS"]==nvis,"MAX_CROSS_COUNT"]=tempdx["Q_MAX_COUNT"].max()
            vis_animal.loc[vis_animal["NEO_VIS"]==nvis,"MAX_CROSS_FREQ"]=tempdx["Q_MAX_FREQ"].max()
        cc_out_file=alignDict_cc["cross_corrFile"]
        cc_out_file_extra=cc_out_file.replace('.tsv','_errors_included.tsv')
        temp_vis_animal=vis_animal[colnam_animals[tname]+["NEO_VIS","VALIDITY","DROP_VIS","MAX_CROSS_COUNT","MAX_CROSS_FREQ"]].copy()
        temp_vis_animal.loc[temp_vis_animal["DROP_VIS"]!=True].to_csv(cc_out_file,sep='\t',index=False)
        temp_vis_animal.to_csv(cc_out_file_extra,sep='\t',index=False)
        xlwriter = pd.ExcelWriter(cc_out_file.replace('.tsv','.xlsx'))
        temp_vis_animal.loc[temp_vis_animal["DROP_VIS"]!=True].to_excel(xlwriter,index=False)
        xlwriter.save()
