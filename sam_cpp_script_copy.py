import pandas as pd
import sqlite3
import datetime as dt
import re
import numpy as np
from joblib import Parallel, delayed
import multiprocessing
from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
import gzip
import math
from types import *
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import subprocess, sys
from pyfaidx import Fasta
import gzip
from Bio.SeqIO import FastaIO
import itertools
from itertools import repeat
import platform


#******** make fxn, pass df
def readSAMwriteDB(alignDict):
	os_name = platform.system()
	#********Create SQL database ***********
	vis_db = sqlite3.connect(':memory:')
	vis_cursor=vis_db.cursor()
	animal_name=alignDict["AnimalID"]
	samDirectory=alignDict["samDirectory"]
	animalDirectory = alignDict["animalDirectory"]
	blastDirectory = alignDict["blastDirectory"]
	genomeFastaPath=alignDict["genomeFastaPath"]
	sampleIDList=alignDict["sampleIDList"]
	Restrict_site=alignDict["Restrict_site"]
	
	
	samFileList=[samDirectory+"noheader_%s.sam.gz"%(smid) for smid in sampleIDList]
	print(samDirectory)
	sam_count=0
	shrtFileList=[]
	print(samFileList)
	for sam_file in samFileList:
		j = 0
		index_start = 1
		print("Processing : ",sam_file)
		sampleID=(sam_file.replace('.sam.gz','')).split('_')[-1]
		if (sampleID not in sampleIDList):
			print("Error: sam filename mismatch")
			sys.exit(0)
		table_name=sampleID
		shrtFileList.append(table_name)
		vis_cursor.execute("DROP TABLE IF EXISTS {tn}".format(tn=table_name))
		pos_col_names=['QNAME','STRAND','CHRID','VIS','FPOS','RPOS','TLEN','FID','RID','READ_NO','TYPE','MULTI_HIT','TOP_MAP','FSEQ','RSEQ']
		chunksize = 20000
		start = dt.datetime.now()
		#****Running the cpp program to process sam file*********
		if (os_name=="Darwin"):
			zcatCmd='gzcat '+sam_file+' | ./sam_reader.out'
		else:
			zcatCmd='zcat '+sam_file+' | ./sam_reader.out'
		p1=subprocess.Popen(zcatCmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE, shell=True)
		#****Reading output of sam_reader.out c++ script *********
		for df in pd.read_table(p1.stdout,sep='\t', header=None,names=pos_col_names,engine='c',chunksize=chunksize, iterator=True,index_col=False):
		    df.index += index_start
		    print ('{} seconds: completed {} rows'.format((dt.datetime.now() - start).seconds, j*chunksize))
		    df.loc[:,pos_col_names].to_sql(table_name, vis_db, if_exists='append')
		    index_start = df.index[-1] + 1
		    j+=1
		    
		print(table_name)
		print("CREATE INDEX {} ON {} (QNAME)".format(table_name,table_name))

		in_cr=vis_cursor.execute("CREATE INDEX {} ON {} (QNAME)".format(table_name+"qin",table_name))
		in_cr.fetchone()
		vis_char=pd.read_sql_query('SELECT CHRID FROM {tn}'.format(tn=table_name), vis_db)
		uni_char=[str(jj) for jj in vis_char.loc[:,"CHRID"].unique()]
		quali_total_VIS=pd.DataFrame(columns=['VIS','Count'])
		unqualified_total=pd.DataFrame(columns=['VIS','Count'])
		total_qname_dict={}
		uni_char.sort()
		if '*' in uni_char:
		    uni_char.remove('*')
		unqualified_seq={}
		for cr in uni_char:
		    for strnd in ['-','+']:
		        qname_dict={}
		        vis_data=pd.read_sql_query('SELECT QNAME,STRAND,CHRID,VIS,FPOS,RPOS,TLEN,FID,RID,READ_NO,TYPE,MULTI_HIT,TOP_MAP FROM {tn} WHERE STRAND==? and CHRID==? and TLEN != 0 and MULTI_HIT != 1'.format(tn=table_name), vis_db,params=[strnd,cr])
		        #vis_data=pd.read_sql_query('SELECT QNAME,STRAND,CHRID,VIS,FPOS,RPOS,TLEN,FID,RID,READ_NO,TYPE,MULTI_HIT,TOP_MAP FROM {tn} WHERE STRAND==? and CHRID==? '.format(tn=table_name), vis_db,params=[strnd,cr])
		        if len(vis_data.index)<1:
		            print ('No VIS found on %s strand of %s'%(strnd,cr))
		            continue
		        vis_data=vis_data.sort_values(by="VIS")
		        vis_data.index=vis_data.QNAME.map(str)+'|'+vis_data.READ_NO.map(str)+'|'+vis_data.TOP_MAP.map(str)
		        uni_vis=vis_data["VIS"].unique()
		        vis_group=vis_data.groupby(['CHRID','TYPE','VIS']).groups
		        vis_num=[[jj,str(str(jj[0])+'|'+strnd+'|'+str(jj[1])),int(jj[2]),int(len(vis_group[jj]))] for jj in vis_group.keys()]
		        my_vis=pd.DataFrame(vis_num,columns=["KEY","TYPE","site","Count"])
		        my_vis.index=my_vis.KEY
		        
		        my_vis_sorted=my_vis.iloc[:,1:].sort_values(by='site',ascending=True)
		        
		        my_vis_sorted["STATUS"]=list(itertools.repeat(0,len(my_vis_sorted.index)))
		        ind1=my_vis_sorted.index.tolist()[0]
		        merged_vis=my_vis_sorted.groupby(["TYPE","site"]).sum() 
		        merg_group=my_vis_sorted.groupby(["TYPE","site"]).groups
		        unqualified_vis=merged_vis.loc[merged_vis["Count"]<=1]
		        qualified_vis=merged_vis.loc[merged_vis["Count"]>1]
		        if len(qualified_vis.index)<1:
		            continue
		        qualified_vis.loc[:,"VIS"]=pd.Series(qualified_vis.index.tolist(),index=qualified_vis.index)
		        for vis1 in qualified_vis.loc[:,"VIS"].tolist():
		            qlist=[]                
		            for key1 in merg_group[vis1].tolist():
		                qlist=qlist+(vis_group[key1].tolist())
		            qname_dict[vis1]=qlist
		        total_qname_dict.update(qname_dict)
		        quali_total_VIS=quali_total_VIS.append(qualified_vis[["VIS","Count"]])           
		        if len(unqualified_vis.index)>0:
		            unqualified_vis.loc[:,"VIS"]=pd.Series(qualified_vis.index.tolist(),index=qualified_vis.index)
		            unqualified_total=unqualified_total.append(unqualified_vis)
		            for ky in unqualified_vis.index:
		                unqualified_v=merg_group[ky]
		                for seq_key in unqualified_v:
		                    seq_list=vis_group[seq_key]
		                    #print seq_key, len(seq_list)
		                    unqualified_seq[seq_key]=seq_list
		                    #unqualified_total=unqualified_total.append(unqulified_vis)
		#cr=0
		seq_name_list=[]
		seq_name_list1=[]
		seq_name_list2=[]
		for seqn in [xval[1] for xval in total_qname_dict.items()]:
		    seq_name_list=seq_name_list+seqn
		
		records1=[]
		records1_shrt=[]
		records2=[]
		records2_shrt=[]
		#******Collecting VIS seqeunces for BLAST remapping******
		for seqname in seq_name_list:
		    #rec1,rec2=select_reads(table_name,seqname)
		    seqid,type1,read_no,top_m=seqname.split('|')
		    qury_seqid=seqid+'|'+type1
		    dx = pd.read_sql_query('SELECT FSEQ,RSEQ FROM {tn} WHERE QNAME==?'.format(tn=table_name), vis_db,params=[qury_seqid])
		    
		    fastaseq1=str(dx['FSEQ'].tolist()[0])
		    fastaseq2=str(dx['RSEQ'].tolist()[0])
		   # print(fastaseq1)
		   # print(fastaseq2)
		    rec1=SeqRecord(Seq(fastaseq1,IUPAC.ambiguous_dna),id=seqid,description="")
		    rec1.id=rec1.id+'|'+read_no+top_m+'|'+type1
		    
		    #****Separating long and short seqeunces********
		    rec2=SeqRecord(Seq(fastaseq2,IUPAC.ambiguous_dna),id=seqid,description="")
		    rec2.id=rec2.id+'|'+read_no+top_m+'|'+type1
		    if len(rec1.seq)<50:
		        records1_shrt.append(rec1)
		    else:
		        records1.append(rec1)
		    if len(rec2.seq)<50:
		        records2_shrt.append(rec2)
		    else:
		        records2.append(rec2) 

#modified
		pathPrefix=blastDirectory+sampleID
		R1_path=pathPrefix+"_R1.fasta"
		R2_path=pathPrefix+"_R2.fasta"
		R1_shrt_path=pathPrefix+"_shrt_R1.fasta"
		R2_shrt_path=pathPrefix+"_shrt_R2.fasta"

		handle_records1=open(R1_path, 'w')
		handle_records2=open(R2_path, 'w')
		handle_records1_shrt=open(R1_shrt_path, 'w')
		handle_records2_shrt=open(R2_shrt_path, 'w')
		#*****Writing the VIS seqeunces to fasta files **********
		fasta_out = FastaIO.FastaWriter(handle_records1, wrap=None)
		fasta_out.write_file(records1)
		handle_records1.close()
		fasta_out = FastaIO.FastaWriter(handle_records2, wrap=None)
		fasta_out.write_file(records2)
		handle_records2.close()
		fasta_out = FastaIO.FastaWriter(handle_records1_shrt, wrap=None)
		fasta_out.write_file(records1_shrt)
		handle_records1_shrt.close()
		fasta_out = FastaIO.FastaWriter(handle_records2_shrt, wrap=None)
		fasta_out.write_file(records2_shrt)
		handle_records2_shrt.close()
		unmapped_seq=pd.read_sql_query("SELECT QNAME,CHRID FROM {tn} WHERE CHRID=='*' ".format(tn=table_name),vis_db)
		vis_data_multi=pd.read_sql_query('SELECT QNAME,MULTI_HIT FROM {tn} WHERE MULTI_HIT == 1 '.format(tn=table_name), vis_db)
		vis_data_read1=pd.read_sql_query('SELECT QNAME,READ_NO FROM {tn} WHERE READ_NO == 1 '.format(tn=table_name), vis_db)

		f = open(sam_file.replace('.sam.gz','.txt'), 'w')
		
#not writing bug

		f.write("Total number of sequences %f\n"%len(vis_char))
		f.write("Total unmapped of sequences %f\n"%len(unmapped_seq.index))
		f.write("Total Multi-mapped of sequences %f\n"%len(vis_data_multi.index))
		f.write("Total Forward read sequences %f\n"%len(vis_data_read1.index))
		f.write("Percentage of unmapped sequences %f\n"%(100.0*len(unmapped_seq.index)/len(vis_char)))
		print("Total number of sequences %f\n"%len(vis_char))
		print("Total unmapped of sequences %f\n"%len(unmapped_seq.index))
		print("Total Multi-mapped of sequences %f\n"%len(vis_data_multi.index))
		print("Total Forward read sequences %f\n"%len(vis_data_read1.index))
		f.write("Percentage of unmapped sequences %f\n"%(100.0*len(unmapped_seq.index)/len(vis_char)))
		print ("Percentage of unmapped sequences %f"%(100.0*len(unmapped_seq.index)/len(vis_char)))
		f.write("Percentage of mapped reads in VIS %f\n" %(100.0*sum(quali_total_VIS.Count)/len(vis_char)))
		print ("Percentage of mapped reads in VIS %f" %(100.0*sum(quali_total_VIS.Count)/len(vis_char)))

		quali_total_VIS.columns=['VIS',table_name]
		quali_file_name=sam_file.replace('.sam.gz','.tsv') #probably correct
		quali_total_VIS.to_csv(quali_file_name,sep='\t',index=False)
		if sam_count==0:  
		    all_quali_VIS=quali_total_VIS.copy()
		else:
		    all_quali_VIS = all_quali_VIS.merge(quali_total_VIS,on='VIS',how='outer')
		sam_count=sam_count+1 

	all_quali_VIS.fillna(value=0,inplace=True)
	tsvFilePath=samDirectory+"vis_all_samples.tsv"
	all_quali_VIS.to_csv(tsvFilePath,sep='\t')
	all_quali_VIS["TOTAL"]=all_quali_VIS[shrtFileList].sum(axis=1)
	all_quali_VIS[["CHRID","SITE"]]=all_quali_VIS['VIS'].astype(str).apply(eval).apply(pd.Series)
	chrid_vis=all_quali_VIS["CHRID"].unique().tolist()
	sel_col=shrtFileList+['TOTAL']
	improved_quali_VIS=pd.DataFrame(columns=['CHRID','SITE','GLEN','GSEQ']+sel_col)
	
	#*****Loading host genome to data base************
	genome_index=Fasta(genomeFastaPath,sequence_always_upper=True,one_based_attributes=False)
	nStep = 2000

	chr1_list=[chr1 for chr1 in chrid_vis]
	#******Extracting host DNA sequence from VIS to nearest restriction site***********
	for chrid in chr1_list:
	    chromo,strand,typ=chrid.split('|')
	    temp_vis=all_quali_VIS.loc[all_quali_VIS["CHRID"]==chrid,].copy()
	    temp_vis["STATUS"]=list(itertools.repeat(0,len(temp_vis.index)))
	    vis_range=5
	    while len(temp_vis["STATUS"].loc[temp_vis["STATUS"]==0].index)!=0:
	        first_site=temp_vis.loc[temp_vis["STATUS"]==0,"SITE"].astype('int').tolist()[0]
	        first_count=temp_vis.loc[temp_vis["STATUS"]==0,"TOTAL"].astype('int').tolist()[0]
	        seter=1
	        while seter:
	            temp = temp_vis.loc[(temp_vis["STATUS"]==0)&((temp_vis["SITE"]>=(first_site-vis_range))&(temp_vis["SITE"]<=(first_site+vis_range))),].copy()
	            if len(temp.index)==1:
	                temp_vis.loc[(temp_vis["SITE"]==first_site)&(temp_vis["TOTAL"]==first_count),"STATUS"]=1
	                seter=0
	                continue
	            else:
	                temp_second_count=max(temp.TOTAL)
	                tmp2=temp.loc[(temp["STATUS"]==0)&(temp["TOTAL"]==temp_second_count)].sort_values(by=['SITE'],ascending=False,axis=0)
	                second_site,second_count=tmp2[["SITE","TOTAL"]].values.tolist()[0] 
	                if first_count == max(temp.TOTAL) and first_site==second_site:
	                            tmp_second_count=max(temp.TOTAL)
	                            temp_vis.loc[(temp_vis["STATUS"]==0)&((temp_vis["SITE"]>=first_site-vis_range)&(temp_vis["SITE"]<=first_site+vis_range)),["SITE","STATUS"]]=(first_site,1) 
	                            seter=0
	                            continue
	                else:
	                    first_site,first_count=second_site,second_count
	    uniq_sites=temp_vis["SITE"].astype('int').unique().tolist()
	    for usite in uniq_sites:
	        if strand == '+':
	            tstart=int(usite)-1
	            seqstr1,seqstr2=str(genome_index[chromo][tstart:tstart+5]),''
	            for iccurent in range(tstart+5,genome_index[chromo][:].end,nStep):
	                seqstr2 +=str(genome_index[chromo][iccurent:iccurent+nStep])
	                if Restrict_site in seqstr2: 
	                    if seqstr2.find(Restrict_site,20)<20: 
	                        continue 
	                    else: 
	                        break
	            cut_loc=seqstr2.find(Restrict_site,20)
	            t_seq=seqstr1+seqstr2[:cut_loc+4]
	            cut_end=tstart+len(t_seq)
	        if strand=='-':
	             tstart=int(usite)-1
	             seqstr1,seqstr2=str(-genome_index[chromo][tstart-5:tstart]),''
	             for iccurent in range(tstart-5,0,-nStep):
	                 seqstr2 +=str(-genome_index[chromo][iccurent-nStep:iccurent])
	                 if Restrict_site in seqstr2: 
	                    if seqstr2.find(Restrict_site,20)<20: 
	                        continue
	                    else:
	                        break
	             cut_loc=seqstr2.find(Restrict_site,20)
	             t_seq=seqstr1+seqstr2[:cut_loc+4]
	             cut_end=tstart-cut_loc+len(t_seq)-4
	        tmps=temp_vis.loc[(temp_vis["SITE"]==usite),sel_col].copy()
	        g_seq=t_seq
	        
	        append_row=[chrid,usite,len(t_seq),g_seq]+tmps.sum(axis=0).tolist()
	        df=pd.Series(append_row).to_frame().transpose()
	        df.columns=improved_quali_VIS.columns
	        improved_quali_VIS=improved_quali_VIS.append(df,ignore_index=True)
	#*******Writing select VIS table data for each sample*****
	bwaFilePath=samDirectory+animal_name+"_bwa_vis_all_samples.tsv"
	improved_quali_VIS.loc[(improved_quali_VIS["TOTAL"]>10)].to_csv(bwaFilePath,sep='\t')
	Ref_VIS_seqrecords=[]
	for ro in improved_quali_VIS.itertuples():
	    rox=[str(x) for x in ro]
	    recid=str(rox[1])+'|'+str(rox[2])+'|'+str(rox[3])
	    record=SeqRecord(Seq(rox[4],IUPAC.ambiguous_dna),id=recid,description='')
	    Ref_VIS_seqrecords.append(record)
	#*******Writing the VIS seqeunce data for one animal*****
	handle_records=open(blastDirectory+'Ref_VIS_seqs.fasta', 'w')
	fasta_out = FastaIO.FastaWriter(handle_records, wrap=None)
	fasta_out.write_file(Ref_VIS_seqrecords)    
	handle_records.close()

	vis_db.close()

#***************
