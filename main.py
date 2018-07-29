import sys,os,time,commands
import second_step
config={'trimmomatic_path':'','mod_path':'','virus_db':'','raw2clean':0}
fin=open('virus_analysis.conf')
for line in fin:
	line=line.strip()
	name=line.split('=')[0]
	value=line.split('=')[1].replace("'","")
	config[name]=value
fin.close()

fin=open('file.list')
total_sample_num=0
for line in fin:
	line=line.strip()
	total_sample_num+=1
	fq1=line
	fq2=fq1.replace('_1.clean.fq.gz','_2.clean.fq.gz')
	sample_name=fq1.split('/')[-1].split('_')[0]
	sh_name=sample_name+'_first_step.sh'
	fout=open(sh_name,'w')
	fout.write("#!/bin/sh\n")
	cmd='python getvirusgene.py '+fq1+' '+fq2+' '+config['trimmomatic_path']+' '+config['mod_path']+' '+config['raw2clean']
	fout.write("%s\n" %cmd)
	fout.close()
	os.system('yhbatch '+sh_name)
fin.close()	
os.system('mv *_first_step.sh 0tmp/sh_file')

###second step for total unigene 
while(num<total_sample_num):
	num=int(commands.getstatusoutput('ls 0tmp/done_file/*predict.done | wc -l')[1])
	time.sleep(15)
os.system('cat 3prediction/*nucl.fa > 3prediction/total_nucl.fasta')

ss=second_step.second_step('3prediction/total_nucl.fasta')
ss.changename()
ss.cdhit()
#start mysql service for species annotation
os.system("cd /HOME/anjiemed_1/mysql")
os.system("mysql_install_db --defaults-file=/HOME/anjiemed_1/.my.cnf --basedir=/HOME/anjiemed_1/mysql")
os.system("mysqld_safe --defaults-file=/HOME/anjiemed_1/.my.cnf --ledir=/HOME/anjiemed_1/mysql/bin &")
os.system("sleep 15")
os.system("cd -")
ss.species_annotation(config['virus_db'])

####last step for LCA species annotation and generate results(abundance and counts) for plot
fin=open('file.list')
#total_unigene='3prediction/totalgene/unigene.fa'
for line in fin:
	line=line.strip()
	fq1=line
	fq2=fq1.replace('_1.clean.fq.gz','_2.clean.fq.gz')
	sample_name=fq1.split('/')[-1].split('_')[0]
	sh_name=sample_name+'_LCA.sh'
	fout=open(sh_name,'w')
	fout.write("#!/bin/sh\n")
	cmd='python last_step.py '+fq1+' '+fq2+' 3prediction/totalgene/unigene.fa'
	fout.write("%s\n" %cmd)
	fout.close()
	os.system('yhabtch '+sh_name)
fin.close()
os.system('mv *LCA.sh 0tmp/sh_file')
