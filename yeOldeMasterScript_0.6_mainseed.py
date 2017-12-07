#!/usr/bin/env python
import os, datetime, glob, time, argparse, subprocess





#################### RNAaseq pipeline for Rabt assembly of deep sequening reads #########
##      This python script will write a .pbs for every fastq file      #
##      in a directory.  it will call tophat, cufflinks, cuffmerge and HTseq.  #
##      These will be output to the same directory branch as FastQ
##      assumes single-end reads
##      n and N for aprun will always be 1.                 #
##	 Open source. By:Carlos PerezCervantes. Moskowitz lab # 
########################################################################################




parser = argparse.ArgumentParser(add_help=True)
	
def readable_dir(string):
	if not os.path.isdir(string):
		raise argparse.ArgumentTypeError("%s directory not found" % string)
	return os.path.abspath(string)


parser.add_argument("-d", "--fastqDir", type=readable_dir,
		help="directory of Fastq files")
				  
parser.add_argument("-o", "--out", type=str,
		help="specify output destination for all steps in the pipeline")
				  
parser.add_argument("-g", "--genome", type=str,default='mm10',
		help="specify genome: mm9, mm10, hg19, hg38")
				  
parser.add_argument("-p", "--processors", type=str, default='32',
		help="specify number of processors for computation node (a mutiple of 32) and tophat") 

parser.add_argument("-w", "--walltime", type=str, default='8', 
		help="specify running time in hours, i.e 02, 4, 16")

parser.add_argument("-l", "--libraryType", type=str, default='firststrand', 
		help="specify library type, options are 'firststrand','secondstrand', or 'unstranded' ")
		
args = parser.parse_args()
parent = os.path.abspath(args.out)

if not os.path.exists(args.out):
    os.makedirs(args.out)


###############################################################################
##
##  assign annotation from user input
##
###############################################################################


if args.genome=='mm9':
    genome_index =' /lustre/beagle2/xyang/Ivan/Carlos/Reference_Sequences/Mouse/mm9/Bowtie2Index/genome '
    genome =' /lustre/beagle2/xyang/Ivan/Carlos/Reference_Sequences/Mouse/mm9/Bowtie2Index/genome.fa '
    gtf= '/lustre/beagle2/xyang/Ivan/Carlos/Reference_Sequences/Mouse/mm9/genes.gtf '
    bed = '/lustre/beagle2/xyang/Ivan/Carlos/Reference_Sequences/Mouse/mm9/genes.mRNA.bed ' 

elif args.genome=='mm10':

    genome_index =' /lustre/beagle2/xyang/Ivan/Carlos/Reference_Sequences/Mouse/release_M9/Sequence/Bowtie2Index/GRCm38.p4.genome '
    genome =' /lustre/beagle2/xyang/Ivan/Carlos/Reference_Sequences/Mouse/release_M9/Sequence/Bowtie2Index/GRCm38.p4.genome.fa '
    gtf=' /lustre/beagle2/xyang/Ivan/Carlos/Reference_Sequences/Mouse/release_M9/Annotation/gencode.vM9.annotation.gtf'
    bed=' /lustre/beagle2/xyang/Ivan/Carlos/Reference_Sequences/Mouse/release_M9/Annotation/gencode.vM9.gene_protein_coding.gtf'

elif args.genome=='hg19':

    genome =' /lustre/beagle2/xyang/Ivan/Carlos/Reference_Sequences/Human/FASTA/GRCh37.v19/GRCh37.p13.genome.fa '
    gtf='/lustre/beagle2/xyang/Ivan/Carlos/Reference_Sequences/Human/Annotation/gencode.v19.annotation.gtf'
    bed='/lustre/beagle2/xyang/Ivan/Carlos/Reference_Sequences/Human/Annotation/gencode.v19.gene_protein_coding.gtf'

elif args.genome == 'hg38':

    genome =' /lustre/beagle2/xyang/Ivan/Carlos/Reference_Sequences/Human/FASTA/GRCh38.primary_assembly.genome.fa '
    gtf='/lustre/beagle2/xyang/Ivan/Carlos/Reference_Sequences/Human/Annotation/gencode.v25.primary_assembly.annotation.gtf'
    bed='/lustre/beagle2/xyang/Ivan/Carlos/Reference_Sequences/Human/Annotation/gencode.v25.gene_protein_coding.gtf'

else:
    print("specify mm9, mm10, hg19, or hg38")

##############################################################################
##  make directories for output
##
##############################################################################


#parent = os.path.dirname(os.path.abspath(args.out))


parentList= glob.glob(os.path.join(args.fastqDir, '[!Undetermined]*.fastq.gz'))


mainseed = list(set([x.strip().split('/')[-1].strip().split('_')[0] for x in parentList]))


tophatDir = parent+'/tophat_'+args.genome+'/'

if not os.path.exists(tophatDir):
    os.makedirs(tophatDir)


tophatDir = os.path.abspath(tophatDir)


all_bam = parent+'/all_bam_'+args.genome+'/'

if not os.path.exists(all_bam):
   	os.makedirs(all_bam)



all_bam = os.path.abspath(all_bam)



cleaned_bam = parent+'/cleaned_bam_'+args.genome+'/'

if not os.path.exists(cleaned_bam):
   	os.makedirs(cleaned_bam)

cleaned_bam = os.path.abspath(cleaned_bam)




cuffl = parent+'/cufflinks_'+args.genome+'/'

if not os.path.exists(cuffl):
   	os.makedirs(cuffl)

cuffl = os.path.abspath(cuffl)


if not os.path.exists(cuffl+'/cuffmerge/'):
	os.makedirs(cuffl+'/cuffmerge/')


cmerge = cuffl+'/cuffmerge/' 


htseq = parent+'/htseq_'+args.genome+'/'

if not os.path.exists(htseq):
   	os.makedirs(htseq)

htseq = os.path.abspath(htseq)



##############################################################################
##
##   make pbs and log directories
##
##############################################################################


if not os.path.exists(tophatDir+'/log/'):
    os.makedirs(tophatDir+'/log/')

thJob = tophatDir+'/log/'


if not os.path.exists(all_bam+'/log/'):
    os.makedirs(all_bam+'/log/')

mergeJob = all_bam+'/log/'


if not os.path.exists(cleaned_bam+'/log/'):
        os.makedirs(cleaned_bam+'/log/')

qcjob = cleaned_bam+'/log/'


if not os.path.exists(cuffl+'/log/'):
        os.makedirs(cuffl+'/log/')


cuffPBS = cuffl+'/log/'



if not os.path.exists(cmerge+'/log/'):
        os.makedirs(cmerge+'/log/')

cmPBS = cmerge+'log/'


if not os.path.exists(htseq+'/log/'):
        os.makedirs(htseq+'/log/')

htPBS = htseq+'/log/'



##############################################################################
##  define function that will be looped to write many 
##  PBS files
##############################################################################


def pbs_writer(pbsDir, seed, processors, walltime, shell):
    pbsf=open(pbsDir+seed+'.pbs','w')
    pbsf.write('#!/bin/bash' + ' \n')
    pbsf.write('#PBS -N ' +seed+ '\n')
    pbsf.write('#PBS -S /bin/bash' +'\n')
    pbsf.write('#PBS -l mppwidth='+processors+ '\n')
    pbsf.write('#PBS -l walltime='+walltime+':00:00' + '\n')
    pbsf.write('#PBS -j oe  \n')
    pbsf.write('#PBS -o ' + pbsDir +seed+'.log'+'\n')
    pbsf.write('#PBS -e ' + pbsDir +seed+'.err' + '\n')
    pbsf.write('./opt/modules/default/init/bash '+'\n')
    pbsf.write('cd $PBS_O_WORKDIR'+'\n')
    pbsf.write('module swap PrgEnv-cray PrgEnv-gnu'+'\n')
    pbsf.write('module load samtools/1.5 '+'\n')
    pbsf.write('aprun -n 1 -N 1 -d ' +processors+' ' +shell)
    pbsf.close()


##############################################################################
##
##  run tophat2 on FastQ file directory.  Also merge and rename
##
##############################################################################
   

for fs in parentList:
    f=fs.strip().split('/')[-1].strip().split('.')[0]
    pbs_writer(pbsDir = thJob, seed=f.replace("-",""), processors=str(args.processors), walltime=str(args.walltime), shell = '/lustre/beagle2/cperez5/soft/shellScripts/tophat2_runner.sh  -o ' +tophatDir+ ' -g '+genome_index+ '  -f ' +fs+  ' -c ' +str(args.processors)+ '  -s '+f)

samples = {}

for fastqzfile in parentList:
  f=fastqzfile.strip().split('/')[-1].strip().split('.')[0]
  name=f.strip().split('_')[0]
  if not ("Undetermined" in f):
        if not (name in samples):
            samples[name]=[]
        src = tophatDir +'/'+ f +'/accepted_hits.bam'
        samples[name].append(src)

for name in samples:
       dst = all_bam +'/'+ name +'.bam'
       if len(samples[name]) == 1:
           os.system('cp ' + samples[name][0] + ' ' + dst)
       if len(samples[name]) > 1:
           file_list = ' '.join(samples[name])
           pbs_writer(pbsDir = mergeJob , seed = name, walltime = str(args.walltime), processors = str(args.processors), shell = ' samtools merge -f ' + dst + ' ' + file_list)
    



##############################################################################
##	Index and clean bam files
##
##############################################################################


for fs in mainseed:
    pbs_writer(pbsDir = qcjob,  seed='qc_'+fs, walltime = str(args.walltime), processors=str(args.processors), shell = 'python /lustre/beagle2/cperez5/soft/shellScripts/samtoolsQC_runner.sh -o ' +cleaned_bam+ ' -f ' +all_bam+'/'+fs+ '.bam -s '+fs)


##############################################################################
##	Run Cuff links in Rabt mode
##
##############################################################################


for fs in mainseed:
    pbs_writer(pbsDir = cuffPBS, seed='cl_'+fs, walltime=str(args.walltime), processors=str(args.processors),  shell = '/lustre/beagle2/cperez5/soft/shellScripts/cufflinks2_runner.sh -o ' +cuffl+  '  -f ' +cleaned_bam+'/'+fs+  '.bam -c ' +str(args.processors)+ '  -s  '+fs+ ' -a ' +gtf )


##############################################################################
##	Run cuffmerge and bedops/bedtools to parse merged GTF
##
##############################################################################

for i in mainseed:
    if not os.path.exists(cuffl+'/'+i):
        os.makedirs(cuffl+'/'+i)

toappend = "/transcripts.gtf"

cuffList= glob.glob(os.path.join(cuffl, '*[!log][!cuffmerge]*'))

gtfList = [x + toappend for x in cuffList]


if not os.path.isfile(cuffl+'/assembly_GTF_list.txt'):
	with open(cuffl+'/assembly_GTF_list.txt', 'a') as outfile:
		outfile.write('\n'.join(gtfList) + '\n')
	mergeList = cuffl+'/assembly_GTF_list.txt'
else:
	mergeList = cuffl+'/assembly_GTF_list.txt'


pbs_writer(pbsDir= cmPBS, seed='cuffmerge', walltime='01', processors=str(args.processors), shell = ' /lustre/beagle2/cperez5/soft/shellScripts/cuffmerge_runner.sh  -o ' +cmerge+ ' -g '+genome+ '  -f ' +mergeList +  ' -c '+str(args.processors)+ ' -a ' +gtf+ ' ' )

pbs_writer(pbsDir= cmPBS, seed='bedopsBedtools', walltime='01', processors=str(args.processors), shell = '/lustre/beagle2/cperez5/soft/shellScripts/bedopsBedtools_runner.sh  -o ' +cmerge+ ' -g '+bed+ ' -f ' +cmerge+'merged.gtf ')



##############################################################################
##	Run HTseq on mergedGTF
##
##############################################################################


for fs in mainseed:
    pbs_writer(pbsDir= htPBS, seed='ht_'+fs, walltime=str(args.walltime), processors=str(args.processors), shell = '/lustre/beagle2/cperez5/soft/shellScripts/htseq_runner.sh -t bam -i gene_id -r reverse -a ' + cmerge+'merged.stranded.gtf ' + '-o ' +htseq+ ' -s ' +fs+ ' -f '+cleaned_bam+'/'+fs+'.bam')

##############################################################################
##
##  Run all qsub commands
##
##############################################################################



thPBSlist= glob.glob(os.path.join(thJob, '*.pbs'))
cflist=[]
thlist=[]

if os.path.isfile(args.out+"/submitALL.sh"):
  os.remove(args.out+"/submitALL.sh")
with open(args.out+'/submitALL.sh', "a") as runnit:
    runnit.write("#!/bin/bash\n \n \n")

for fs in thPBSlist:
  f=fs.strip().split('/')[-1].replace("_R1_001.pbs","")
  thlist.append('$th'+f)
  thcmd = 'th'+f+'=$(qsub '+fs+ ')\necho -e $th'+f+ ' \n'
  
  with open(args.out+'/submitALL.sh', "a") as runnit:
    runnit.write(thcmd+ '\nsleep 2\n\n\n\n')
    runnit.close()
#print("\n \n \n submitting all tophat jobs \n \n \n ")
for name in mainseed:
  bmcmd='mrg'+name.replace("-","")+'=$(qsub -W depend=afterok:'+ ":".join([str(i) for i in thlist])+' '+mergeJob+name+'.pbs) \n echo -e $mrg'+name.replace("-","")+ ' \n'
  qccm='qc'+name.replace("-","")+'=$(qsub -W depend=afterok:$mrg'+name.replace("-","")+' '+qcjob+'qc_'+name+'.pbs) \n echo -e $qc'+name.replace("-","")+ ' \n'
  clcm='cl'+name.replace("-","")+'=$(qsub -W depend=afterok:$qc'+name.replace("-","")+' '+cuffPBS+'cl_'+name+'.pbs) \n echo -e $cl'+name.replace("-","") + ' \n'
  cflist.append('$cl'+name.replace("-",""))
  cmrcmd='cmerge=$(qsub -W depend=afterok:'+":".join([str(i) for i in cflist])+' '+cmPBS+'cuffmerge.pbs) \n echo -e $cmerge' + ' \n'
  bbcmd='btools=$(qsub -W depend=afterok:$cmerge '+cmPBS+ 'bedopsBedtools.pbs) \n echo -e $btools ' + '\n'
  
  with open(args.out+'/submitALL.sh', "a") as runnit:
    runnit.write(bmcmd+ '\nsleep 2\n\n\n\n')
    runnit.write(qccm+ '\nsleep 2\n\n\n\n')
    runnit.write(clcm+ '\nsleep 2\n\n\n\n')
    runnit.close()

#print("\n \n \n submitting all BAM merge and rename jobs \n \n \n ")
#print("\n \n \n QC jobs \n \n \n ")
#print("\n \n \n submitting cufflinks \n \n \n ")


with open(args.out+'/submitALL.sh', "a") as runnit:
  runnit.write(cmrcmd+ '\nsleep 2\n\n\n\n')
  runnit.write(bbcmd+ '\nsleep 2\n\n\n\n')
  runnit.close()


#print("\n \n \n submit cuffmerge and bedtools/bedops \n \n \n ")

for name in mainseed:
  htscmd='hts'+name.replace("-","")+'=$(qsub -W depend=afterok:$btools ' +htPBS+'ht_'+name+ '.pbs) \n echo -e $hts'+name.replace("-","") + ' \n'
  with open(args.out+'/submitALL.sh', "a") as runnit:
    runnit.write(htscmd+ '\nsleep 2\n\n\n\n')
#print("\n \n \n submit HTSeq \n \n \n ")
os.system('chmod u+x ' +args.out+"/submitALL.sh")
#os.system( +args.out+"/submitALL.sh")





######################################
##
##  Get mapping stats and plot with R
##
######################################


pbs_writer(pbsDir = thJob, seed = 'tophSummry' , walltime= '01', processors= str(args.processors), shell='/lustre/beagle2/cperez5/soft/python2Scripts/get_tophat_summary.py -d ' +tophatDir)

for fs in mainseed:
  pbs_writer(pbsDir = qcjob,  seed='sm_'+fs, walltime ='01', processors=str(args.processors), shell = '/lustre/beagle2/cperez5/soft/shellScripts/samtools_bamsummary.sh -o ' +cleaned_bam+ '/bam_summary.txt -f ' +cleaned_bam+ '/' +fs+'.bam' )


pbs_writer(pbsDir = htPBS, seed = 'htsqSummry' , walltime= '01', processors= str(args.processors), shell=' /lustre/beagle2/cperez5/soft/python2Scripts/htseq_summary.py -d '+htseq)

