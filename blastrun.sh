#!/bin/bash
echo $1 sample mapping start
/opt/latesblast/ncbi-blast-2.7.1+/bin/./blastn -query $1_R1.fasta -db $2 -perc_identity 95 -max_target_seqs 500 -num_threads 8 -outfmt "6 qseqid sseqid pident evalue qstart qend sstart send qlen" |awk '{split($1,a,"|"); split($2,b,"|"); if($5 == 1 && $7 <= 2 && a[3]==b[3] && ($6-($5-1))/$9 > 0.9) print $0 }' |gzip -c >$1_output1.gz &
/opt/latesblast/ncbi-blast-2.7.1+/bin/./blastn -query $1_shrt_R1.fasta -db $2 -task blastn-short -perc_identity 95 -max_target_seqs 500 -num_threads 8 -outfmt "6 qseqid sseqid pident evalue qstart qend sstart send qlen" |awk '{split($1,a,"|"); split($2,b,"|"); if($5 == 1 && $7 <= 2 && a[3]==b[3] && ($6-($5-1))/$9 > 0.9) print $0 }' |gzip -c >$1_shrt_output1.gz &
/opt/latesblast/ncbi-blast-2.7.1+/bin/./blastn -query $1_R2.fasta -db $2 -perc_identity 95 -max_target_seqs 500 -num_threads 8 -outfmt "6 qseqid sseqid pident evalue qstart qend sstart send qlen" |awk '{split($1,a,"|"); split($2,b,"|"); if(a[3]==b[3] && ($6-$5)/$9 > 0.9) print $0 }' |gzip -c >$1_output2.gz &
/opt/latesblast/ncbi-blast-2.7.1+/bin/./blastn -query $1_shrt_R2.fasta -db $2 -task blastn-short -perc_identity 95 -max_target_seqs 500 -num_threads 8 -outfmt "6 qseqid sseqid pident evalue qstart qend sstart send qlen" |awk '{split($1,a,"|"); split($2,b,"|"); if(a[3]==b[3] && $5 == 1 && ($6-$5)/$9 > 0.9) print $0 }' |gzip -c >$1_shrt_output2.gz &
wait
echo $1 sample mapping done
