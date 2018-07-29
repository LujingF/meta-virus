import MySQLdb,sys
from optparse import OptionParser
def mysql_annot(blast_out,annot_out):
	#this function for meta analysis annotation
	con=MySQLdb.connect(host="localhost",user="username",passwd="password",db="database",unix_socket="socket_path")
	cur=con.cursor()
	fin=open(blast_out)
	fout=open(annot_out,'w')
	gene=''
	ID=[]
	ids=[]
	category=['species','genus','family','order','class','phylum','kingdom','root']
	Category=['S__','G__','F__','O__','C__','P__','K__','R__']
	for line in fin:
		line=line.strip()
		gene=line.split('\t')[0]
#		accession_num=str(line.split('\t')[1].split('|')[3])
		accession_num=str(line.split('\t')[2])
#		print accession_num
#		accession_num=str(line.split('\t')[13])
		cur.execute("select * from acc2taxid where prot_acc_id = '{0}' limit 1;".format(accession_num))
		if int(cur.rowcount)>0:
			tax_id=cur.fetchone()[1]
			if tax_id==0:
				pass
				# genes which unassigned
#				print gene+'\tunassign' 
			else:
				num=1
				results=[]
				while(num<40):
					num+=1
					cur.execute("select * from taxid2name where tax_id = '{0}' limit 1;".format(tax_id))
#					tax_name=cur.fetchone()[1].replace(' ','_')
					if int(cur.rowcount)>0:
						tax_name=cur.fetchone()[1].replace(' ','_')
						cur.execute("select * from nodes where tax_id= '{0}' limit 1;".format(tax_id))
						if int(cur.rowcount)>0:
							content=cur.fetchone()
							former_id=content[1]
							rank=content[2]
							tax_id=former_id
							if rank=='species':
								result='S__'+tax_name
								results.append(result)
							elif rank=='genus':
								result='G__'+tax_name
								results.append(result)
							elif rank=='family':
								result='F__'+tax_name
								results.append(result)
							elif rank=='order':
								result='O__'+tax_name
								results.append(result)
							elif rank=='class':
								result='C__'+tax_name
								results.append(result)
							elif rank=='phylum':
								result='P__'+tax_name
								results.append(result)
							elif rank=='kingdom' or rank=='superkingdom' or rank=='subkingdom':
								result='K__'+tax_name
								results.append(result)
							
							if tax_name=='root':
								result='R__root'
								results.append(result)
								break
						else:
							break
					
				if len(results)<8:
					Tmp_dict={}
					for i in results:
						category_name=i.split('__')[0]+'__'
						name=i.split('__')[1]
						Tmp_dict[category_name]=name
#					print Tmp_dict
					tmp=[]
					for k in Category:
						if k in Tmp_dict:
							tmp.append(k+Tmp_dict[k])
						else:
							tmp.append(k+'undef')
					results=tmp		
				Result=';'.join(results[::-1])
				fout.write("%s\t%s\n" %(gene,Result))
#				print gene+'\t'+Result
				
	fin.close()
	fout.close()

def species_annot(blast_out,annot_out):
	#function: to annotation transcriptome
	con=MySQLdb.connect(host="localhost",user="username",passwd="password",db="database",unix_socket="sock_path")
	cur=con.cursor()
	fin=open(blast_out)
	fout=open(annot_out,'w')
	sp_count={}
	for line in fin:
		line=line.strip()
		gene=line.split()[0]
		accession_num=str(line.split('\t')[1])
		cur.execute("select * from acc2taxid where prot_acc_id = '{0}' limit 1;".format(accession_num))
		if int(cur.rowcount)>0:
			tax_id=cur.fetchone()[1]
			cur.execute("select * from taxid2name where tax_id = '{0}' limit 1;".format(tax_id))
			if int(cur.rowcount)>0:
				tax_name=cur.fetchone()[1]
				if tax_name not in sp_count:
					sp_count[tax_name]=1
				else:
					sp_count[tax_name]+=1
#			fout.write("%s\t%s\n" %(gene,tax_name))
#			print gene+'\t'+tax_name
	for key in sp_count:
		fout.write("%s\t%i\n" %(key,sp_count[key]))
	fin.close()	
	fout.close()

def main():
	usage='usage: %prog [options] arg'
	p=OptionParser()
	p.add_option('-i','--input',type='string',dest='input',action='store',help='input file must be blast or diamond result and file format must be tab,the first column must be gene name and second column must be accession number')
	p.add_option('-m','--meta',type='int',dest='meta',action='store',help='this option for meta data annotation, if set 0, annotation result will only display species name, if set 0 annotation result will display root,kingdom,phylum,class,order,family,genus,species name')
	p.add_option('-o','--out',type='string',dest='output',action='store',help='output file')
	(opt,args)=p.parse_args()
	blast_out=opt.input
	annot_out=opt.output
	meta_flag=opt.meta
	if meta_flag==1:
		mysql_annot(blast_out,annot_out)
	else:
		species_annot(blast_out,annot_out)

if __name__ == "__main__":
	main()
