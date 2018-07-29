# Meta Virus analysis scripts
These scripts for Meta Virus analysis. It can be used to detected virus species from your NGS sequence data. These scripts was desined at slurm cluster.  
****
## Set enviroment for scripts  
Required tools: [Python2.7](https://www.python.org/download/releases/2.7/),  [GeneMark](http://opal.biology.gatech.edu/GeneMark/), [megahit](https://github.com/voutcn/megahit), [mysql](http://mirrors.sohu.com/mysql/MySQL-5.5/mysql-5.5.42-linux2.6-x86_64.tar.gz), [trimmoatic](https://github.com/timflutre/trimmomatic), [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/), [Quast](http://bioinf.spbau.ru/quast), [lefse](https://bitbucket.org/biobakery/biobakery/wiki/lefse), [samtools](http://samtools.sourceforge.net/), [krona](https://github.com/marbl/Krona/wiki)   
R packages required: [vegan](http://cc.oulu.fi/~jarioksa/softhelp/vegan.html), [reshape2](https://cran.r-project.org/web/packages/reshape2/index.html), [ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html), [RColorBrewer](https://cran.r-project.org/web/packages/RColorBrewer/index.html), [ggthemes](https://cran.r-project.org/web/packages/ggthemes/index.html), [ggsignif](https://cran.r-project.org/web/packages/ggsignif/index.html), [stringr](https://cran.r-project.org/web/packages/stringr/index.html), [pheatmap](https://cran.r-project.org/web/packages/pheatmap/index.html)  
## Usage  
### Files prepared  
#### file.list
Note: file list file was needed, and each line represent one of paired end files, example:  
sample1_1.clean.fq.gz  
sample2_1.clean.fq.gz  
...  
The file which contain sequence file name must be '__file.list__', and the name of sequence file must be contain __'_1.clean.fq.gz'__, if your data is raw data, please rename your file like that and set raw2clean as 1.   
#### virus_analysis.conf
- __trimmomatic_path__, adaptor sequence file should be under trimmomatic path.  
- __mode_path__, module file path for GeneMark.  
- __virus_db__, virus database path deposited virus protein database(download from ncbi),  you must use blast or diamond to create index for virus database
- __raw2clean__, if your data is raw data, then you can set this as  1, otherwise, set this value as 0.  
- __mysql_user__, this is for mysql database user name  
- __mysql_password__, mysql password  
- __mysql_database__, select sub database from mysql
- __unix_socket__, this would be used for python package MySQLdb

#### import ncbi species, accession number files to mysql  
Download species information, protein accession number etc from ncbi. Three tables should be in the database, the first table named 'acc2taxid', and this table store protein accession number(column name is 'prot_acc_id') and its related taxonomy id information(column name is tax_id). The second table named 'taxid2name' and this table also has two column: 'tax_id' and 'tax_name'. The last table named 'nodes' which contains three columns: 'tax_id','former_tax_id','rank'. 'tax_id' means taxonomy id in ncbi. 'former_tax_id' means the parent taxonomy id, 'rank' means species classification: kingdom, phylum, class, order, family, genus, species.   
## Virus Annotation
Make sure getvirusgene.py, second_step.py, last_step.py, main.py, main.sh and scripts catalogue in your work folder.  
Run command: __sbatch main.sh__  
## alpha diversity analysis
Put counts result into a directory, then you can use 'alpha_diversity.R' in scripts to do alpha diversity analysis.  
## lefse difference analysis  
Firstly, you should use the annotation result for lefse analysis to generate a matrix data for lefse analysis,
then you can use 'lefse_analysis.sh' to do lefse difference analysis.  
## Beta diversity analysis
'beta_diversity.R' is for beta diversity analysis. If you want to use this script, you shold create a directory named 'speciescounts' and put species counts result file in this directory. 
