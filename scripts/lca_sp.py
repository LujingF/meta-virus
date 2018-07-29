import sys,os,commands
from optparse import OptionParser
#this script for virus analysis
def cal_abundance(alignment):
	#calculate gene abundance,alignment could be sam or bam file
	sample_name=alignment.split('/')[-1].split('.')[0]
	sample_readscount=sample_name+'.counts'
	sample_abundance=sample_name+'.abundance'
	cmd='samtools view -@ 24 -SF 260 '+alignment+ '| cut -f 3| sort | uniq -c > '+sample_readscount
#	os.system(cmd)
	cmd='samtools view -@ 24 -SF 260 '+alignment+'| cut -f 3| sort | uniq -c | awk '+"'{split($2,a,"+r'"|");printf("%s\t%s\t%e\n", $2,$1,$1/a[2])}'+"' > "  +sample_abundance
#	os.system(cmd)

	return sample_abundance

def LCA(annotation_result,sample_abundance):
	fin=open(sample_abundance)
	gene_counts={}
	gene_relative_abundance={}
	for line in fin:
		line=line.strip()
		gene=line.split()[0]
		counts=int(line.split()[1])
		abundance=float(line.split()[2])
		if gene not in gene_counts:
			gene_counts[gene]=counts
		else:
			gene_counts[gene]+=counts
		if gene not in gene_relative_abundance:
			gene_relative_abundance[gene]=abundance
		else:
			gene_relative_abundance[gene]+=abundance
	fin.close()
	fin=open(annotation_result)
	gene=''
	results={}
	annot_count={}
	annot_abundance={}
	num=0
	names=['Root','Kingdom','Phylum','Class','Order','Family','Genus','Species']
	for line in fin:
		line=line.strip()
		gene_name=line.split()[0]
		annot=line.split()[-1]
		Root=annot.split(';')[0]
		Kingdom=annot.split(';')[1]
		Phylum=annot.split(';')[2]
		Class=annot.split(';')[3]
		Order=annot.split(';')[4]
		Family=annot.split(';')[5]
		Genus=annot.split(';')[6]
		Species=annot.split(';')[7]
		if gene!=gene_name:
			if len(results)!=0:
				Result=''
				for i in xrange(0,len(names)):
					element_num=len(list(set(results[names[i]])))#count non-redundant element in list, if the element_num is 1, it means all annotation result in this level is same
					annot_result=list(set(results[names[i]]))[0]
					if element_num == 1:
						Result+=(annot_result+'\t')
					else:
						break
				Result=Result[:-1]
				num+=1
				if gene in gene_counts:
					counts=int(gene_counts[gene])
					abundance=float(gene_relative_abundance[gene])
				else:
					counts=0
					abundance=0.0
				if Result not in annot_count:
					annot_count[Result]=counts
				else:
					annot_count[Result]+=counts
				if Result not in annot_abundance:
					annot_abundance[Result]=abundance
				else:
					annot_abundance[Result]+=abundance
			gene=gene_name
			results={'Root':[],'Kingdom':[],'Phylum':[],'Class':[],'Order':[],'Family':[],'Genus':[],'Species':[]}
		results['Root'].append(Root)
		results['Kingdom'].append(Kingdom)
		results['Phylum'].append(Phylum)
		results['Class'].append(Class)
		results['Order'].append(Order)
		results['Family'].append(Family)
		results['Genus'].append(Genus)
		results['Species'].append(Species)
	fin.close()
	if len(results)	!=0:
		Result=''
		for i in xrange(0,len(names)):
			element_num=len(list(set(results[names[i]])))
			annot_result=list(set(results[names[i]]))[0]
			if element_num == 1:
				Result+=(annot_result+'\t')
			else:
				break
		Result=Result[:-1]
		num+=1
		if gene in gene_counts:
			counts=int(gene_counts[gene])
		else:
			counts=0
		if Result not in annot_count:
			annot_count[Result]=counts
		else:
			annot_count[Result]+=counts
		if gene in gene_relative_abundance:
			abundance=float(gene_relative_abundance[gene])
		else:
			abundance=0.0
		if Result not in annot_abundance:
			annot_abundance[Result]=abundance
		else:
			annot_abundance[Result]+=abundance
	total_abundance=0
	for key in annot_abundance:
		total_abundance+=float(annot_abundance[key])
	species_relative_abundance={}
	for key in annot_abundance:
		species_relative_abundance[key]=annot_abundance[key]/total_abundance
#		print key,annot_abundance[key]/total_abundance
	return annot_count,species_relative_abundance

def dict_count(key,value,dict):
	if key not in dict:
		dict[key]=value
	else:
		dict[key]+=value
		
def generate_result(species_counts,species_abundance,prefix_name):
	krona=prefix_name+'.krona'
	lefse_species=prefix_name+'_species.lefse'
	lefse_genus=prefix_name+'_genus.lefse'
	lefse_family=prefix_name+'_family.lefse'
	species_plot=prefix_name+'_species.abundance'
	genus_plot=prefix_name+'_genus.abundance'
	family_plot=prefix_name+'_family.abundance'
	total_sp_count=prefix_name+'_total_species.count'
	sp_count=prefix_name+'_species.count'
#	fout=open(krona,'w')
	dicts={'Roots':{},'Kingdoms':{},'Phylums':{},'Classs':{},'Orders':{},'Familys':{},'Genuss':{},'Speciess':{}}
	dict_order=['Roots','Kingdoms','Phylums','Classs','Orders','Familys','Genuss','Speciess']
	Species_ab={}
	Genus_ab={}
	Family_ab={}
	Species_count={}
	krona_dict={}
	count_dict={}
#	fout1=open(total_sp_count,'w')
	for key in species_abundance:
#		fout.write("%s\t%s\n" %(species_abundance[key],key))
		if len(key.split())==8:
			name=key
		elif len(key.split())==7:
			name=key+'\tS__undef'
		elif len(key.split())==6:
			name=key+'\tG__undef\tS__undef'
		elif len(key.split())==5:
			name=key+'\tF__undef\tG__undef\tS__undef'
		elif len(key.split())==4:
			name=key+'\tO__undef\tF__undef\tG__undef\tS__undef'
		if name not in krona_dict:
			krona_dict[name]=species_abundance[key]
		else:
			krona_dict[name]+=species_abundance[key]
		if name.replace('\t',';') not in count_dict:
			count_dict[name.replace('\t',';')]=species_counts[key]
		else:
			count_dict[name.replace('\t',';')]+=species_counts[key]
#		fout.write("%s\t%s\n" %(species_abundance[key],name))
#		fout1.write("%s\t%s\n" %(name.replace('\t',';'),species_counts[key]))
		species_name=name.split()[-1]
		genus_name=name.split()[-2]
		family_name=name.split()[-3]
		order_name=name.split()[-4]
		class_name=name.split()[-5]
		phylum_name=name.split()[-6]
		kingdom_name=name.split()[-7]
		root_name=name.split()[-8]
		Root=root_name
		Kingdom=Root+'|'+kingdom_name
		Phylum=Kingdom+'|'+phylum_name
		Class=Phylum+'|'+class_name
		Order=Class+'|'+order_name
		Family=Order+'|'+family_name
		Genus=Family+'|'+genus_name
		Species=Genus+'|'+species_name
		value=float(species_abundance[key])
		dict_count(species_name,value,Species_ab)
		dict_count(genus_name,value,Genus_ab)
		dict_count(family_name,value,Family_ab)
		dict_count(species_name,int(species_counts[key]),Species_count)
		for k in dicts:
			Key=k[:-1]
			dict_count(eval(Key),value,dicts[k])
#	fout.close()
#	fout1.close()
	fout=open(krona,'w')
	for key in krona_dict:
		fout.write("%s\t%s\n" %(krona_dict[key],key))
	fout.close()
	fout=open(total_sp_count,'w')
	for key in count_dict:
		fout.write("%s\t%s\n" %(key,count_dict[key]))
	fout.close()
	fout=open(species_plot,'w')
	for key in Species_ab:
		fout.write("%s\t%s\n" %(key,Species_ab[key]))	
	fout.close()
	fout=open(genus_plot,'w')
	for key in Genus_ab:
		fout.write("%s\t%s\n" %(key,Genus_ab[key]))
	fout.close()
	fout=open(family_plot,'w')
	for key in Family_ab:
		fout.write("%s\t%s\n" %(key,Family_ab[key]))
	fout.close()
	fout=open(lefse_species,'w')
	for i in xrange(0,len(dict_order)):
		for k in dicts[dict_order[i]]:
			fout.write("%s\t%s\n" %(k,dicts[dict_order[i]][k]))

	fout.close()
	fout=open(lefse_genus,'w')
	for i in xrange(0,len(dict_order)-1):
		for k in dicts[dict_order[i]]:
			fout.write("%s\t%s\n" %(k,dicts[dict_order[i]][k]))
	fout.close()
	fout=open(lefse_family,'w')
	for i in xrange(0,len(dict_order)-2):
		for k in dicts[dict_order[i]]:
			fout.write("%s\t%s\n" %(k,dicts[dict_order[i]][k]))
	fout.close()
	fout=open(sp_count,'w')
	for key in Species_count:
		fout.write("%s\t%i\n" %(key,Species_count[key]))
	fout.close()
def main():
	usage='usage: %prog [options] arg'
	p=OptionParser()
	p.add_option('-s','--species',type='string',dest='species',action='store',help="input species annotation result,the first column is gene name and the second column is species annotation result")
	p.add_option('-a','--alignment',type='string',dest='alignment',action='store',help='input reads mapped gene set bam or sam alignment result, samtools should be installed in your computer ')
#	p.add_option('-o','--outdir',type='string',dest='alignment',action='store',help='output dir')
	(opt,args)=p.parse_args()
	species=opt.species
	alignment=opt.alignment
#	out=opt.outdir
	sample_abundance=cal_abundance(alignment)
	species_counts,species_abundance=LCA(species,sample_abundance)
#	for key in species_counts:
#		print key,species_abundance[key]
	prefix_name=alignment.split('/')[-1].split('.')[0].split('_')[0]
	generate_result(species_counts,species_abundance,prefix_name)

if __name__ == "__main__":
	main()
