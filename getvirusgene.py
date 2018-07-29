import sys,os,commands
class virus():
	def __init__(self,fq1,fq2,output):
		self.fq1=fq1
		self.fq2=fq2
		self.output=output
	@staticmethod
	def check():
		#this function is check the enviroment
		#fastqc,megahit,gmhmmp,bowtie2,diamond,cd-hit-est
		tools=['fastqc','megahit','gmhmmp','bowtie2','diamond','cd-hit-est']
		num=0
		for i in tools:
			check_result=commands.getstatusoutput('which '+i)[1]
			if check_result.startswith('which'):
				print i+' was not installed on your computer or not added into your enviroment path'
				exit(1)
			else:
				num+=1
		if num==len(tools):
			return 1
		else:
			return 0
	@staticmethod
	def check_done(done_file):
		if os.path.exists(done_file):
			return 0
		else:
			return 1

	@staticmethod
	def mkdir():
		os.system("mkdir -p 0tmp/done_file 0tmp/sh_file 1quality_control/cleandata 1quality_control/fastqc_result 2assembly/megahit_result 2assembly/quast_result 3prediction/mappingsample 3prediction/totalgene 3prediction/index 4abundance 5annotation 6result/lefse 6result/krona 6result/counts")
	def raw2clean(self,trimmomatic_path):
		if self.fq1.count('/')!=0:
			tmp=self.fq1.split('/')
			fq1_name=tmp[-1]
			basename=fq1_name.split('_1')[0]
			del(tmp[-1])
			path='/'.join(tmp)+'/'
			fq2_name=self.fq2.split('/')[-1]
		else:
			path='./'
			fq1_name=self.fq1
			basename=fq1_name.split('_1')[0]
			fq2_name=self.fq2
		cmd='java -jar '+trimmomatic_path+'trimmomatic-0.36.jar PE -threads 24 '+self.fq1+' '+self.fq2+' -baseout 1quality_control/cleandata/'+basename+' ILLUMINACLIP:'+trimmomatic_path+'adapters/TruSeq2-PE.fa:2:30:10:8:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:120'
		os.system(cmd)
		trimmomatic1='1quality_control/cleandata/'+basename+'_1P'
		cleanfq1=trimmomatic1.replace('_1P','_1.fq')
		trimmomatic2='1quality_control/cleandata/'+basename+'_2P'
		cleanfq2=trimmomatic2.replace('_2P','_2.fq')
		os.system('mv '+trimmomatic1+' '+cleanfq1)
		os.system('mv '+trimmomatic2+' '+cleanfq2)
		self.fq1=cleanfq1
		self.fq2=cleanfq2
		return self.fq1,self.fq2

	def fastqc(self):
		cmd='fastqc -t 24 -o 1quality_control/fastqc_result '+self.fq1
		os.system(cmd)
		cmd='fastqc -t 24 -o 1quality_control/fastqc_result '+self.fq2
		os.system(cmd)
	def assembly(self,status):
		if status:
			sample_name=self.fq1.split('/')[-1].split("_1")[0]
			cmd='megahit -t 24 -m 0.95 --out-dir 2assembly/megahit_result/'+sample_name+'  --out-prefix '+sample_name+' -1 '+self.fq1+' -2 '+self.fq2
			os.system(cmd)
			assembly_result='2assembly/megahit_result/'+sample_name+'/'+sample_name+'.contigs.fa'
			quast_cmd='quast.py -t 24 -o 2assembly/quast_result/'+sample_name+' '+assembly_result
	#		print quast_cmd
			os.system(quast_cmd)
			self.output=assembly_result
			return self.output
	def predict(self,mod_path,status):
		if status:
			sample_name=self.output.split('/')[2]
			genemark_cmd='gmhmmp '+self.output+' -f G -A '+'3prediction/'+sample_name+'.prot.fa -D 3prediction/'+sample_name+'.nucl.fa -m '+mod_path+' -o 3prediction/'+sample_name+'.gff'
			os.system(genemark_cmd)
			os.system('touch 0tmp/done_file/'+sample_name+'_predict.done')
	def mapping(self,status):
		if status:
			sample_name=self.output.split('/')[2]
			build_cmd='bowtie2-build 3prediction/'+sample_name+'.nucl.fa 3prediction/index/'+sample_name
			os.system(build_cmd)
			align_cmd='bowtie2 -x 3prediction/index/'+sample_name+' -1 '+self.fq1+' -2 '+self.fq2+' -p 24 -S 3prediction/mappingsample/'+sample_name+'.sam 2> 3prediction/mappingsample/'+sample_name+'_align.result'
#			print align_cmd
			os.system(align_cmd)
			os.system('rm 3prediction/mappingsample/'+sample_name+'.sam')

fq1=sys.argv[1]
fq2=sys.argv[2]
trimmomatic_path=sys.argv[3]
mod_path=sys.argv[4]
flag=sys.argv[5]
Virus=virus(fq1,fq2,'')
if(Virus.check()):
	Virus.mkdir()
	if flag==1:
		Virus.raw2clean(trimmomatic_path)
	Virus.fastqc()
	status=Virus.check_done(Virus.assembly(1))
	status=Virus.check_done(Virus.predict(mod_path,status))
	Virus.mapping(status)
