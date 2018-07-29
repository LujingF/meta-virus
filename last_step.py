class last_step():
	def __init__(self,total_gene,fq1,fq2):
		self.total_gene=total_gene
		self.fq1=fq1
		self.fq2=fq2

	def mapping(self):
		sample_name=fq1.split('/')[-1].split('_')[0]
		align_cmd='bowtie2 -x '+self.total_gene+' -1 '+self.fq1+' -2 '+self.fq2+' -p '+str(cpu_count())+' -S 3prediction/totalgene/'+sample_name+'.sam'
		os.system(align_cmd)
		return '3prediction/totalgene/'+sample_name+'.sam'
	@staticmethod
	def lca(align_result):
		#generate gene annotation result according to LCA algorithm 
		cmd='python scripts/lca_sp.py -s 5annotation/total_species_annot.tab -a '+align_result
		os.system(cmd)
		sample_name=align_result.split('/')[-1].split('.sam')[0]
		os.system("mv "+sample_name+"*.abundance 4abundance/")
		os.system("mv "+sample_name+"*.lefse 6result/lefse")
		os.system("mv "+sample_name+"*krona 6result/krona")
		os.system("mv "+sample_name+"*count 6result/counts")


fq1=sys.argv[1]
fq2=sys.argv[2]
total_gene=sys.argv[3]
ls=last_step(total_gene,fq1,fq2)
align_result=ls.mapping()
ls.lca(align_result)
