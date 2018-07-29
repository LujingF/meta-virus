import sys,os
from multiprocessing import cpu_count

class second_step():
	#function: this step including get unigene(collect all predicted genes from all samples into a total predicted genes file, rename the each sequence and remove redundant sequences)
	#mapping: align sample reads to unigenes and calculate unigenes expression value(unigene abundance and unigene reads count)
	#annotation: annotate unigene, and get species classification of each unigene
	def __init__(self,total_gene):
		self.total_gene=total_gene
	def changename(self):
	#function: there may be duplicated sequence name since all predicted results were  combined into a single file
		unigene='3prediction/totalgene/totalgene.fa'
		fin=open(self.total_gene)
		fout=open(unigene,'w')
		num=0
		for line in fin:
			line=line.strip()
			if len(line)!=0:
				if line[0]=='>':
					num+=1
					length=line.split('=')[-1]
					name='>gene_'+str(num)+'|'+length
					fout.write("%s\n"%name)
				else:
					fout.write("%s\n"%line)
		fin.close()
		fout.close()
		self.total_gene=unigene
		return self.total_gene
	def cdhit(self):
	#function: use cdhit to remove redundant sequence
		cdhit_cmd='cd-hit-est -i 3prediction/totalgene/totalgene.fa -o 3prediction/totalgene/unigene.fa -T 0 -M 0'
		os.system(cdhit_cmd)
		build_cmd='bowtie2-build 3prediction/totalgene/unigene.fa 3prediction/totalgene/unigene.fa'
		os.system(build_cmd)
		self.total_gene='3prediction/totalgene/unigene.fa'
		return self.total_gene


	def species_annotation(self,virus_db):
		diamond_cmd='diamond blastx -q '+self.total_gene+' -d '+virus_db+' -o 5annotation/total_virus_annot.tab  -f 6 qseqid evalue sseqid qlen'
		os.system(diamond_cmd)
		fin=open('5annotation/total_virus_annot.tab')
		fout=open('5annotation/total_virus_annot_filter.tab','w')
		query_name=''
		for line in fin:
			line=line.strip()
			query_id=line.split()[0]
			evalue=line.split()[1]
			if line.split()[2].count('|')!=0:
				if line.split()[2].split('|')[1].count('.')==0:
					acc=line.split()[2].split('|')[1]+'.1'
				else:
					acc=line.split()[2].split('|')[1]
			else:
				if line.split()[2].split('|')[1].count('.')==0:
					acc=line+'.1'
				else:
					acc=line
			query_len=line.split()[3]
			if query_id!=query_name:
				query_names=query_id
				threshold=evalue*10
				fout.write("%s\t%s\t%s\t%s\n" %(query_id,evalue,acc,query_len))
			else:
				if evalue<=threshold:
					fout.write("%s\t%s\t%s\t%s\n" %(query_id,evalue,acc,query_len))
		fin.close()
		fout.close()
		os.system("python scripts/species_mysql_annot.py -i 5annotation/total_virus_annot_filter.tab -m 1 -o 5annotation/total_species_annot.tab")
