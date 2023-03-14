#This code is desinged to read the sample sheet and generate input files for 
# vector and linker trimming c++ program. This script also runs the c++ program

import pandas as pd
import subprocess, sys, os
from pathlib import Path
import itertools
from Bio import SeqIO
from Bio.Seq import Seq
import gzip
from sam_cpp_script_copy import readSAMwriteDB
from blast_runner import blastRunner
from crossover_correction import cross_correct_fun
import datetime
import platform
import shutil
import pickle

#take in sample sheet file path
print ('Number of arguments:', len(sys.argv), 'arguments.')
print ('Argument List:', str(sys.argv))
if len(sys.argv)<3:
    print("$ipython3 VIS_runner_copy_GS Sample_sheet.xlsx run_dir")
    sys.exit(0)
samplesheetname=sys.argv[1]
run_dir=sys.argv[2]#+"_"+datetime.datetime.now()strftime("%d/%m/%Y")
os_name = platform.system()

#*************************##***********************#
#parameters 
Max_Spacer_Length=9

#************************##************************#
#delete old result_out and runinput_dir folders, make new ones
file_list=[]
Complete_sample_sheet=pd.read_excel(samplesheetname)


print("*******************************************************************")
print("Warning: Running the script may delete output from the previous run.\n")
print("Make sure that your output directory is not same as previous run.\n")
print("Do you want to continue y/n:  ")
deci_in=input()
#**************Creating run directory for the analysis************
if deci_in.lower()!='y':
    print("Aborting run")
    sys.exit(0)

if Complete_sample_sheet.empty:
    print("Error:Empty sample information sheet.\n")
    sys.exit(0)

if run_dir[-1:]=="/":
    run_dir=run_dir[:-1]

if os.path.exists(run_dir):
    shutil.rmtree(run_dir)
    os.makedirs(run_dir)
else:
    os.makedirs(run_dir)
#****************************************************************

#*******************User input options***************************

print("To run the script select one of the options:\n")
print("a To run for select AnimalID\n")
print("b To run for select SampleID\n")
print("c To run for select SampleNames\n")
print("d To run for all samples \n")
print("e To exit\n")
first_option=input()
if first_option.lower() =='e':
    print('Exited the program without running\n')
    sys.exit(0)
if first_option.lower() not in ['a','b','c','d']:
    print("Rerun the program and select correct option\n")
    sys.exit(0)
if first_option.lower() =='a':
    print("Enter Animal IDs with , separating each ID. When done press enter")
    animalID_list=input().split(',')
    Select_sample_sheet=Complete_sample_sheet.loc[Complete_sample_sheet['AnimalID'].isin(animalID_list),:].copy()
    if Select_sample_sheet.empty:
        print('Incorrect Animal IDs')
        sys.exit(0)
if first_option.lower() =='b':
    print("Enter Sample IDs with , separating each ID. When done press enter")
    sampleID_list=input().split(',')
    Select_sample_sheet=Complete_sample_sheet.loc[Complete_sample_sheet['SampleID'].isin(sampleID_list),:].copy()
    if Select_sample_sheet.empty:
        print('Incorrect Sample IDs')
        sys.exit(0)
if first_option.lower() =='c':
    print("Enter SampleNames with , separating each ID. When done press enter")
    samplename_list=input().split(',')
    Select_sample_sheet=Complete_sample_sheet.loc[Complete_sample_sheet['SampleName'].isin(samplename_list),:].copy()
    if Select_sample_sheet.empty:
        print('Incorrect Sample names')
        sys.exit(0)
if first_option.lower() =='d':
    print("All samples in the list will be run")
    Select_sample_sheet=Complete_sample_sheet.copy()
    if Select_sample_sheet.empty:
        print('Excel file is emplty or problem with reading file')
        sys.exit(0)
#*****************************************************************

#***********Create lists of Animal Id and Samples IDs********************
Select_animalID=Select_sample_sheet['AnimalID'].unique().tolist()
Select_SampleID=Select_sample_sheet['SampleID'].unique().tolist()
Select_SampleName=Select_sample_sheet['SampleName'].unique().tolist()

#***********Bioinfomatics analysis selections******************
print("Do you want to perform bioinformatics analysis? (y/n)")
bioinfo_check=input()
if bioinfo_check.upper()=="Y":
    bioinfo_run=True
    print("Bioinfomatics analysis will be done")
else:
    print("Bioinfomatics analysis will not be done")
    bioinfo_run=False

#check for empty columns and check if files are there
#**************Writing the sample file in the run directory *********************
Select_sample_sheet.to_csv(run_dir+'/samples_process_this_run.tsv',sep='\t',index=False,doublequote=False)

#*******Choice to provide input file for the VL trimmer ***************
print("Do you want to provide parameter file? (y/n)")
parafile_check=input()
if parafile_check.upper()=='Y':
    print("Provide location of the parameter file and name \t (/home/user/myparameterfile.txt):\n")
    parafile_location=input()
    if not parafile_location:
        print("Location of the parameter file is not provided:")
        sys.exit(0)
    cmd_check_para='ls '+parafile_location
    p=subprocess.Popen(cmd_check_para, stdout=subprocess.PIPE,stderr=subprocess.PIPE, shell =True)
    get_para_list=p.stdout.read().decode('utf8').split('\n')
    get_para_list.remove('')
    if len(get_para_list)!=1:
        print("Error: Missing parameter file")
        sys.exit(0)
    file_list.append('Parameter:'+ parafile_location)   
else:
    print("Trimming will be done using default parameter values\n")

#****************Write date and time analysis started ********************
run_date_pkl=str(datetime.date.today()).replace("-","_")
date_ti= open(run_dir+'/run_date_time.txt', 'w')
date_ti.write(str(datetime.datetime.now()))
date_ti.close()
#animalID case
cc_animal_dir={}


for animal in Select_animalID:
    scstart_time = datetime.datetime.now()
    if first_option.lower() =='d':
        Select_sample_sheet=Complete_sample_sheet.loc[Complete_sample_sheet['AnimalID']==animal,:].copy()
    #*************Checking inputs for genome files ********************** 
    genomeLabelList=Select_sample_sheet['UCSC_Ref_Genome'].unique().tolist()

    if (len(genomeLabelList)>1):
        print("Error: more than one UCSC_Ref_Genome listed")
        sys.exit(0)
    elif (len(genomeLabelList)==0):
        print("Error: missing UCSC_Ref_Genome")
        sys.exit(0)
    genomeLabel=genomeLabelList[0]
    #**********Checking Restriction_Site information ************************
    Restriction_Site_list=Select_sample_sheet['Restriction_Site'].unique().tolist()
    if (len(Restriction_Site_list)>1):
        print("Error: more than one Restriction Site listed for one animal")
        sys.exit(0)
    elif (len(Restriction_Site_list)==0):
        print("Error: missing Restriction_Site")
        sys.exit(0)
    Restriction_Site=Restriction_Site_list[0]
    #***********make output directories, trimmed seq output dir for each animal ***************
    outdir_list=Select_sample_sheet['Output_dir'].unique().tolist()
    if(len(outdir_list)==1):
        outdir=outdir_list[0]
    else:
        #*************Checking output directory details ********************** 
        if (len(outdir_list)>1):
            print("Error: more than one output directories for one animalID not permited")
            sys.exit(0)
        elif (len(genomeLabelList)==0):
            print("Using default dirctory dir_AnimalID")
            outdir='dir_'+animal
    out_put_dir=run_dir+'/'+outdir
    if not os.path.exists(out_put_dir):
        os.makedirs(out_put_dir)
    result_dir=out_put_dir+'/result_out'
    if not os.path.exists(result_dir):
        os.makedirs(result_dir)
    
    #*************Create input directory to store vl trimmer input files for each animal********************** 
    runinput_dir=out_put_dir+'/runinput_dir' #--'alt'
    if not os.path.exists(runinput_dir):
        os.makedirs(runinput_dir)

    #*************Create output dictory for each animal ********************** 
    animal_outdir=result_dir+'/'+animal #str(samplerow.loc[sindex,'AnimalID'])
    check_make_trimdir='mkdir -p '+animal_outdir
    chkmake=subprocess.Popen(check_make_trimdir, stdout=subprocess.PIPE,stderr=subprocess.PIPE, shell =True)
    chkmake.communicate()
    proc = subprocess.Popen('ls -d ./'+animal_outdir, stdout=subprocess.PIPE, shell =True)
    direct_list=proc.stdout.read().decode('UTF-8').split('\n')
    
    #*************Create sub output dictory for each animal ********************** 
    myDirectoryList=['sam_out/','vl_trimmed/','dump/','blast/','final_results/','cross_corrected/']
    for di in myDirectoryList:
        if di not in direct_list:
            cmd1='mkdir -p %s' %animal_outdir+'/'+di
            xxx=subprocess.Popen(cmd1, stdout=subprocess.PIPE, shell =True)
            xxx.communicate()

    sample_sheet=Select_sample_sheet.loc[Select_sample_sheet['AnimalID']==animal].copy().reset_index()
    checked_sampleID_list = []
    for sindex, samplerow in sample_sheet.iterrows():
        samplerow=samplerow.to_frame().transpose()
        checked_sampleID = str(samplerow.loc[sindex,'SampleID'])
        checked_sampleID_list.append(checked_sampleID)
        

        #***create sample_input_dir for each sample in one animal*************8
        sample_input_dir=runinput_dir+'/'+checked_sampleID
        check_make_dir='mkdir -p '+sample_input_dir
        chkmake=subprocess.Popen(check_make_dir, stdout=subprocess.PIPE,stderr=subprocess.PIPE, shell =True)
        chkmake.communicate()
        
        #*************check if the fastq files are present or not ***********
        if pd.isnull(samplerow.loc[sindex,'Location_of_fastq_files']):
            print("Error: Please provide location of fastq files")
            sys.exit(0)
        if pd.isnull(samplerow.loc[sindex,'Location_of_genome_fasta_file']):
            print("Error: Please provide location of reference genome files")
            sys.exit(0)
        cmd_check_fastq='ls '+str(samplerow.loc[sindex,'Location_of_fastq_files']).strip()+'/'+str(samplerow.loc[sindex,'Fastq_File_name_starts']).strip()+'*'
        p=subprocess.Popen(cmd_check_fastq, stdout=subprocess.PIPE,stderr=subprocess.PIPE, shell =True)
        get_fastq_list=p.stdout.read().decode('utf8').split('\n')
        get_fastq_list.remove('')
        if len(get_fastq_list)!=2:
            print("Error: Missing fastq files")
            sys.exit(0)

        #************check for reference genome file *********************
        cmd_check_fa='ls '+str(samplerow.loc[sindex,'Location_of_genome_fasta_file']).strip()+'/'+str(samplerow.loc[sindex,'UCSC_Ref_Genome']).strip()+'*.fa'
        print(cmd_check_fa)
        p=subprocess.Popen(cmd_check_fa, stdout=subprocess.PIPE,stderr=subprocess.PIPE, shell =True)
        get_fa_list=p.stdout.read().decode('utf8').split('\n')
        print(get_fa_list)
        get_fa_list.remove('')
        if len(get_fa_list)!=1 or not get_fa_list:
            print("Error: Missing reference genome fasta files or incorrect format of .fa")
            sys.exit(0)

        #*******************check if probe, vector and linker sequences are provided
        if pd.isnull(samplerow.loc[sindex,'Starting_Vector_sequence']):
            print("Error:No vector sequence.\nPlease provide the vector sequence.\n")
            sys.exit(0)
        if pd.isnull(samplerow.loc[sindex,'Vector_min_%identity']):
            print('Warning: Vector_min_%%identity not provided, using default value of 90.0')
            samplerow.loc[1,'Vector_min_%identity']=90.0 
        if pd.isnull(samplerow.loc[sindex,'Linker_sequence']):
            print("Error:No vector sequence.\nPlease provide the vector sequence.\n")
            sys.exit()
        if pd.isnull(samplerow.loc[sindex,'Linker_min_%identity']):
            print('Warning: Linker_min_%%identity not provided, using default value of 90.0')
            samplerow.loc[sindex,'Linker_min_%identity']=90.0 
        if pd.isnull(samplerow.loc[sindex,'ForwardPrimer_sequence']):
            samplerow.loc[sindex,'ForwardPrimer_sequence']=samplerow.loc[sindex,'Starting_Vector_sequence'].upper()[:18]
            if pd.isnull(samplerow.loc[sindex,'Primer_%identity']):
                samplerow.loc[sindex,'Primer_%identity']=90.0

        Fprimer_seq='>ForwardPrimer|'+str(samplerow.loc[sindex,'Primer_%identity'])+'\n'+samplerow.loc[sindex,'ForwardPrimer_sequence'].upper()+'\n'
        vector_seq='>Vector|'+str(samplerow.loc[sindex,'Vector_min_%identity'])+'\n'+samplerow.loc[sindex,'Starting_Vector_sequence'].upper()+'\n'
        linker_seq='>Linker|'+str(samplerow.loc[sindex,'Linker_min_%identity'])+'\n'+samplerow.loc[sindex,'Linker_sequence'].upper()

        #******************write vector and linker sequences in fasta file*********************
        vec_link_file=open(sample_input_dir+'/'+samplerow.loc[sindex,'Fastq_File_name_starts']+'_vector_linker.fasta','w')
        vec_link_file.write(Fprimer_seq+vector_seq+linker_seq)
        vec_link_file.flush()
        vec_link_file.close()
        file_list.append('Vector_Linker:'+sample_input_dir+'/'+samplerow.loc[sindex,'Fastq_File_name_starts']+'_vector_linker.fasta')

        #***************check for signature mutations************************
        
        if samplerow.loc[sindex,'Number_of_Signature_Mutations']>0 or not pd.isnull(samplerow.loc[sindex,'Number_of_Signature_Mutations']):
            mut_seq_dict={}
            
            mut_seqs=samplerow.loc[sindex,'Sig_Mut_sequences'].split(';')
            mut_seqs = list(filter(None, mut_seqs))
            if pd.isnull(samplerow.loc[sindex,'Sig_Mut_%Identity']):
                for mseq in mut_seqs:
                    if not mseq.split(':')[0] or not mseq.split(':')[1]:
                        print("Error: problem with the signature mutation sequence information")
                        sys.exit(0)
                    else:
                        mut_seq_dict[mseq.split(':')[0]]=[mseq.split(':')[0],mseq.split(':')[1].upper(),'85.0']
            else:
                mut_idt=samplerow.loc[sindex,'Sig_Mut_%Identity'].split(';')
                mut_idt = list(filter(None, mut_idt))
                if len(mut_seqs) != len(mut_idt):
                    print("Error: Elements in column Sig_Mut_sequences and Sig_Mut_%Identity are not equal")
                    sys.exit(0)
                
                uni_keys=[]
                for mseq in mut_seqs:
                    if not mseq.split(':')[0] or not mseq.split(':')[1]:
                        print("Error: problem with the signature mutation sequence information")
                        sys.exit(0)
                    else:
                        mut_seq_dict[mseq.split(':')[0]]=[mseq.split(':')[0],mseq.split(':')[1]]
                if len(mut_seq_dict.keys()) != len(mut_seqs):
                    print("Error: names of the signature mutation are not unique")
                    sys.exit(0)
                for mseq in mut_idt:
                    if ':' not in mseq:
                        if not mseq:
                            print("Error: problem with the signature mutation identity information")
                            sys.exit(0)
                        else:
                            mut_seq_dict[mseq].append('85.0')
                    else:
                        if not mseq.split(':')[1]:
                            mut_seq_dict[mseq].append('85.0')
                        else:
                            mut_seq_dict[mseq.split(':')[0]].append(mseq.split(':')[1])
                if len(mut_seq_dict.keys()) != len(mut_seqs):
                    print("Error: Elements in column Sig_Mut_sequences and Sig_Mut_%%Identity do not match")
                    sys.exit(0)
            
            #*******************write signature mutations in fasta file******************************
            sig_mut_fasta=open(sample_input_dir+'/'+samplerow.loc[sindex,'Fastq_File_name_starts']+"_sig_sequence.fasta",'w')
            for skey in mut_seq_dict.keys():
                sig_mut_fasta.write('>'+mut_seq_dict[skey][0]+'|'+mut_seq_dict[skey][2]+'\n')
                sig_mut_fasta.write(mut_seq_dict[skey][1]+'\n')
            sig_mut_fasta.flush()
            sig_mut_fasta.close()
            file_list.append('Signature_seq:'+sample_input_dir+'/'+samplerow.loc[sindex,'Fastq_File_name_starts']+'_sig_sequence.fasta')
        else:
            print("No signature mutations in the samples")

        #******************check for internal vector sequences***********************
        if not pd.isnull(samplerow.loc[sindex,'Internal_Vector_sequences']):
            inter_seq_dict={}
            inter_seq=samplerow.loc[sindex,'Internal_Vector_sequences'].split(';')
            if pd.isnull(samplerow.loc[sindex,'Internal_Vector_%identity']):
                cnt=0
                for intseq in inter_seq:
                    inter_seq_dict['internal_seq_'+str(cnt)]=[intseq.upper(),'90.0']
                    cnt +=1
            else:
                inter_idt=samplerow.loc[sindex,'Internal_Vector_%identity'].split(';')
                if len(inter_idt)<len(inter_seq):
                    difflen=len(inter_seq)-len(inter_idt)
                    for i in range(sindex,difflen):
                        inter_idt.append('90.0')
                cnt=0
                for intseq in inter_seq:
                    inter_seq_dict['internal_seq_'+str(cnt)]=[intseq.upper(),str(inter_idt[cnt])]
                    cnt +=1

            #****************write internal sequences in fasta file
            interseq_fasta=open(sample_input_dir+'/'+samplerow.loc[sindex,'Fastq_File_name_starts']+"_internal_sequence.fasta",'w')
            for skey in inter_seq_dict.keys():
                interseq_fasta.write('>'+skey+'|'+inter_seq_dict[skey][1]+'\n')
                interseq_fasta.write(inter_seq_dict[skey][0]+'\n')
            interseq_fasta.flush()
            interseq_fasta.close()
            file_list.append('Internal_seq:'+sample_input_dir+'/'+samplerow.loc[sindex,'Fastq_File_name_starts']+'_internal_sequence.fasta')
        
        #****************Check spacer sequences to remove corss contamination 
        if (samplerow.loc[sindex,'Check_spacer_seqs']=='Y'):
            spacer_seq_dict={}
            if pd.isnull(samplerow.loc[sindex,'Forward_spacer_seq']):
                for_spacer_seq = (samplerow.loc[sindex,'Starting_Vector_sequence'].upper()[:Max_Spacer_Length])
                for_spacer_len = 0.0
            else:
                for_spacer_seq = (samplerow.loc[sindex,'Forward_spacer_seq'].upper()+samplerow.loc[sindex,'Starting_Vector_sequence'].upper()[:Max_Spacer_Length-len(samplerow.loc[sindex,'Forward_spacer_seq'])])
                for_spacer_len = 1.0*(len(samplerow.loc[sindex,'Forward_spacer_seq']))
            spacer_seq_dict["Forward_spacer"]=[for_spacer_seq,str(for_spacer_len)]
            
            link_seq = Seq(samplerow.loc[sindex,'Linker_sequence'].upper())
            RC_link_seq = str(link_seq.reverse_complement())
            if pd.isnull(samplerow.loc[sindex,'Reverse_spacer_seq']):
                rev_spacer_seq = (RC_link_seq.upper()[:Max_Spacer_Length])
                rev_spacer_len = 0.0
            else:
                rev_spacer_seq = (samplerow.loc[sindex,'Reverse_spacer_seq'].upper()+RC_link_seq.upper()[:Max_Spacer_Length-len(samplerow.loc[sindex,'Reverse_spacer_seq'])])
                rev_spacer_len = 1.0*len(samplerow.loc[sindex,'Reverse_spacer_seq'])
            spacer_seq_dict["Reverse_spacer"]=[rev_spacer_seq,str(rev_spacer_len)]
            #**********write spacer sequences in fasta file
            spacerseq_fasta=open(sample_input_dir+'/'+samplerow.loc[sindex,'Fastq_File_name_starts']+"_spacer_sequence.fasta",'w')
            for skey in spacer_seq_dict.keys():
                spacerseq_fasta.write('>'+skey+'|'+str(spacer_seq_dict[skey][1])+'\n')
                spacerseq_fasta.write(spacer_seq_dict[skey][0]+'\n')
            spacerseq_fasta.flush()
            spacerseq_fasta.close()
            file_list.append('Spacer_seq:'+sample_input_dir+'/'+samplerow.loc[sindex,'Fastq_File_name_starts']+'_spacer_sequence.fasta')
        #******************identify read1 fastq file
        if not pd.isnull(samplerow.loc[sindex,'Fastq_identifier']):
            read_file_dict={}
            read_identifier=samplerow.loc[sindex,'Fastq_identifier'].split(':')
            read_file_dict["R1"]=[read_identifier[0],'1']
            read_file_dict["R2"]=[read_identifier[1],'2']
            print("Using %s (forward) and %s (reverse) as identifier:\n" %(read_identifier[0],read_identifier[1]))
            for file_n in get_fastq_list:
                
                #R1 read
                if read_file_dict["R1"][0] in file_n:
                    handle1=gzip.open(file_n,"rt")
                    jj=0
                    num_of_chk_read=10
                    file1_seqids={}

                    for record1 in SeqIO.parse(handle1,"fastq"):
                        read_no = record1.description.split(' ')[1].split(':')[0]
                        if(read_no == '1' or read_no == '2'):
                            file1_seqids[read_no]=record1.description 
                            jj=jj+1
                            if jj>(num_of_chk_read-1):
                                break
                        else:
                            print("Read information is missing\n")
                            sys.exit(0)
                    if len(file1_seqids.keys()) > 1:
                        print("Unable to identify file reads\n")
                        sys.exit(0)
                    if read_file_dict["R1"][1] in file1_seqids.keys():
                        file_list.append("R1:"+file_n)
                        forward_fastq_file=file_n.split('/')[len(file_n.split('/'))-1]
                    if read_file_dict["R2"][1] in file1_seqids.keys():
                        file_list.append("R2:"+file_n) 
                        reverse_fastq_file=file_n.split('/')[len(file_n.split('/'))-1]
                
                #R2 read
                if read_file_dict["R2"][0] in file_n:
                    handle1=gzip.open(file_n,"rt")
                    jj=0
                    num_of_chk_read=10
                    file1_seqids={}

                    for record1 in SeqIO.parse(handle1,"fastq"):
                        read_no = record1.description.split(' ')[1].split(':')[0]
                        if(read_no == '1' or read_no == '2'):
                            file1_seqids[read_no]=record1.description 
                            jj=jj+1
                            if jj>(num_of_chk_read-1):
                                break
                        else:
                            print("Read information is missing\n")
                            sys.exit(0)
                    if len(file1_seqids.keys()) > 1:
                        print("Unable to identify file reads\n")
                        sys.exit(0)
                    if read_file_dict["R1"][1] in file1_seqids.keys():
                        file_list.append("R1:"+file_n)
                        forward_fastq_file=file_n.split('/')[len(file_n.split('/'))-1]
                    if read_file_dict["R2"][1] in file1_seqids.keys():
                        file_list.append("R2:"+file_n)
                        reverse_fastq_file=file_n.split('/')[len(file_n.split('/'))-1]
        else:
            print("Trying to identify read1 and read2 fastq files\n")
            handle1=gzip.open(get_fastq_list[0],"rt")
            handle2=gzip.open(get_fastq_list[1],"rt")
            jj=0
            num_of_chk_read=10

            file1_seqids={}

            for record1 in SeqIO.parse(handle1,"fastq"):
                read_no = record1.description.split(' ')[1].split(':')[0]
                if(read_no == '1' or read_no == '2'):
                    file1_seqids[read_no]=record1.description 
                    jj=jj+1
                    if jj>(num_of_chk_read-1):
                        break
                else:
                    print("Read information is missing\n")
                    sys.exit(0)
            if len(file1_seqids.keys()) != 1 :
                    print("Unable to identify file reads\n")
                    sys.exit(0)
            if '1' in file1_seqids.keys():
                    file_list.append("R1:"+get_fastq_list[0])
                    file_n = get_fastq_list[0]
                    forward_fastq_file=file_n.split('/')[len(file_n.split('/'))-1]
            if '2' in file1_seqids.keys():
                    file_list.append("R2:"+get_fastq_list[0])
                    file_n = get_fastq_list[0]
                    reverse_fastq_file=file_n.split('/')[len(file_n.split('/'))-1]
            jj=0
            num_of_chk_read=10
            file2_seqids={}
            for record1 in SeqIO.parse(handle2,"fastq"):
                read_no = record1.description.split(' ')[1].split(':')[0]
                if(read_no == '1' or read_no == '2'):
                    file2_seqids[read_no]=record1.description 
                    jj=jj+1
                    if jj>(num_of_chk_read-1):
                        break
                else:
                    print("Read information is missing\n")
                    sys.exit(0)
            if len(file2_seqids.keys()) != 1:
                    print("Unable to identify file reads\n")
                    sys.exit(0)
            if '1' in file2_seqids.keys():
                    file_list.append("R1:"+get_fastq_list[1])
                    file_n = get_fastq_list[1]
                    forward_fastq_file=file_n.split('/')[len(file_n.split('/'))-1]
            if '2' in file2_seqids.keys():
                    file_list.append("R2:"+get_fastq_list[1])
                    file_n = get_fastq_list[1]
                    reverse_fastq_file=file_n.split('/')[len(file_n.split('/'))-1]

        print("Forward read fastq= %s \n"%forward_fastq_file)
        print("Reverse read fastq= %s \n"%reverse_fastq_file)
        if forward_fastq_file == reverse_fastq_file:
            print("Error in fastq files: contain incorrect read information\n")
            sys.exit(0)
        trim_forward_fastq=animal_outdir+'/'+'vl_trimmed/'+forward_fastq_file.replace('.fastq.gz','_VL_trimmed.fastq')
        trim_reverse_fastq=animal_outdir+'/'+'vl_trimmed/'+reverse_fastq_file.replace('.fastq.gz','_VL_trimmed.fastq')
        file_list.append("outR1:"+trim_forward_fastq)
        file_list.append("outR2:"+trim_reverse_fastq)

        dump_forward_fastq=animal_outdir+'/'+'dump/'+forward_fastq_file.replace('.fastq.gz','_dump.fastq')
        dump_reverse_fastq=animal_outdir+'/'+'dump/'+reverse_fastq_file.replace('.fastq.gz','_dump.fastq')

        file_list.append("dumpR1:"+dump_forward_fastq)
        file_list.append("dumpR2:"+dump_reverse_fastq)
        #**************Wirite vl trimmer input file ***************
        cppinputfile=open(sample_input_dir+'/'+'cppinput_filelist.txt','w')
        cppinputfile.write('\n'.join(file_list))
        cppinputfile.flush() 
        cppinputfile.close()

        run_trim_cmd='./vltrimmer.out '+sample_input_dir+'/'+'cppinput_filelist.txt'
        

        trim_seq=subprocess.Popen(run_trim_cmd, stdout=subprocess.PIPE,stderr=subprocess.PIPE, shell =True)
        stdout_data, stderr_data = trim_seq.communicate()
        print(stdout_data)
        if trim_seq.returncode != 0:
            print("%r failed, stderr %r" % (run_trim_cmd, stderr_data))
            sys.exit(0)
        zipfile_list=[trim_forward_fastq,trim_reverse_fastq,dump_forward_fastq,dump_reverse_fastq]
        for tozip in zipfile_list:
            zipcmd='gzip -f9 '+tozip
            zippro=subprocess.Popen(zipcmd, stdout=subprocess.PIPE,stderr=subprocess.PIPE, shell =True)
            stdout_data, stderr_data = zippro.communicate()
            if zippro.returncode != 0:
                print("%r failed, stderr %r" % (run_trim_cmd, stderr_data))
                sys.exit(0)

        file_list[:]=[]

        #*************run BWA on trim_forward and trim_reverse, insert created SAM file into sam_out folder
        current_working_dir=animal_outdir+'/'+'vl_trimmed/'
        R1=trim_forward_fastq+'.gz'
        R2=trim_reverse_fastq+'.gz'
        sam_file=checked_sampleID+'.sam'
        
        ref_genome_path=get_fa_list[0].replace("//","/")

        aln_command1='bwa aln %s %s> %stemp1.sai' %(ref_genome_path,R1,current_working_dir)
        aln1=subprocess.Popen(aln_command1, stdout=subprocess.PIPE, shell =True) 
        print ('\n *****  aligning %s *****\n' %R1)
        aln1.communicate()

        aln_command2='bwa aln %s %s> %stemp2.sai' %(ref_genome_path,R2,current_working_dir)
        aln2=subprocess.Popen(aln_command2, stdout=subprocess.PIPE, shell =True) 
        print ('\n *****  aligning %s *****\n' %R2)
        aln2.communicate()

        make_sam="bwa sampe -a 1500 -o 1000000 -N 1000 -P %s %stemp1.sai %stemp2.sai %s %s |grep -v '^@' |cut -f11 --complement |sed 's/[A,N,S,X][B-Z0-9]:[A,B,H,Z,i,f]:[a-zA-Z0-9]*\t//g' > %s/sam_out/noheader_%s" %(ref_genome_path,current_working_dir,current_working_dir,R1,R2,animal_outdir,sam_file)
        sam1=subprocess.Popen(make_sam, stdout=subprocess.PIPE, shell =True)
        print ('\n ***** generating sam output for %s and %s ******\n' %(R1,R2))
        sam1.communicate()

        zipsam_cmd='gzip -f9 %s/sam_out/noheader_%s'%(animal_outdir,sam_file)
        zipsam=subprocess.Popen(zipsam_cmd,stdout=subprocess.PIPE, shell =True)
        print('\n Zipping sam file %s'%sam_file)
        zipsam.communicate()
    
    alignDict={}
    samDirectory=result_dir+'/'+animal+'/sam_out/'
    animalDirectory=samDirectory.replace("sam_out/","")

    alignDict["AnimalID"]=animal
    alignDict["samDirectory"]=samDirectory
    alignDict["animalDirectory"]=animalDirectory
    alignDict["blastDirectory"]=animalDirectory+"blast/"
    alignDict["finalDirectory"]=animalDirectory+"final_results/"
    alignDict["cross_corrDirectory"]=animalDirectory+"cross_corrected/"
    alignDict["genomeFastaPath"]=ref_genome_path
    alignDict["sampleIDList"]=checked_sampleID_list
    alignDict["genomeLabel"]=genomeLabel
    alignDict["Restrict_site"]=Restriction_Site
    with open(run_dir+'/'+animal + '.pkl', 'wb') as fani:
        pickle.dump(alignDict, fani, pickle.HIGHEST_PROTOCOL)
    
    #**********Processing of sam files outputed by BWA *************
    readSAMwriteDB(alignDict)
    
    #**********Remapping using BLAST****************
    file_4_cc=blastRunner(alignDict)
    
    alignDict["finalFile"]= file_4_cc
    alignDict["cross_corrFile"]=alignDict["cross_corrDirectory"]+animal+'_CC_Final_vis_data.tsv'
    cc_animal_dir[animal]=alignDict
print("Starting crossover correction")

#*******Write pickel file for the run

with open(run_dir+'/'+run_date_pkl+ '.pkl', 'wb') as f:
    pickle.dump(cc_animal_dir, f, pickle.HIGHEST_PROTOCOL)

#*********Runing crossover correction code ***********
cross_correct_fun(cc_animal_dir,run_dir)
    #print("Time for processing animal %s is %f hrs"%(animal,(datetime.datetime.now() - scstart_time).seconds/3600))
#*********Running bioinformatics analysis ************
if(bioinfo_run):
    print("Now starting bioinformatics analysis")
    bio_animalID=list(cc_animal_dir.keys())
    for bio_ani in bio_animalID:
        bio_alignDict_cc=cc_animal_dir[bio_ani]
        file_4_bioinfo=bio_alignDict_cc["cross_corrFile"]
        genomeLabel=bio_alignDict_cc["genomeLabel"]
        bio_info_cmd="python3 blast_bioinfo.py "+file_4_bioinfo+" "+genomeLabel 
        run_bioinfo=subprocess.Popen(bio_info_cmd, stdout=subprocess.PIPE,stderr=subprocess.PIPE, shell =True)
        run_bioinfo.communicate()
    print("Completed bioinformatics analysis")    
