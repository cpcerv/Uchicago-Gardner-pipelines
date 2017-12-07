#!/usr/bin/env python
import os, datetime, glob, time, argparse, subprocess





####################     Pipeline for ATAC/ChIPseq     ################################
##      This python script will write a .pbs for every fastq file                     #
##      in a directory.  it will call bowtie,samtools and macs2.                      #
##      These will be output to the same directory branch as FastQ                    #
##      can be single end or paired end reads                                         #
##      n and N for aprun will always be 1                                            #
##      Open source. By:Carlos Perez-Cervantes. Moskowitz lab                         # 
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
		
parser.add_argument("-r", "--readType", type=str, default='single_end', 
    help="specify whether these reads are paired or single_end")

parser.add_argument("-e", "--experiment", type=str, required=True, 
    help="specify experiment type, ATACseq or ChIPseq ")


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
    gensize='mm'

elif args.genome=='mm10':

    genome_index =' /lustre/beagle2/xyang/Ivan/Carlos/Reference_Sequences/Mouse/release_M9/Sequence/Bowtie2Index/GRCm38.p4.genome '
    genome =' /lustre/beagle2/xyang/Ivan/Carlos/Reference_Sequences/Mouse/release_M9/Sequence/Bowtie2Index/GRCm38.p4.genome.fa '
    gensize='mm'

elif args.genome=='hg19':

    genome =' /lustre/beagle2/xyang/Ivan/Carlos/Reference_Sequences/Human/FASTA/GRCh37.v19/GRCh37.p13.genome.fa '
    genome_index =' /lustre/beagle2/xyang/Ivan/Carlos/Reference_Sequences/Human/FASTA/GRCh37.v19/GRCh37.p13.genome '
    gensize='hs'


elif args.genome == 'hg38':

    genome =' /lustre/beagle2/xyang/Ivan/Carlos/Reference_Sequences/Human/FASTA/GRCh38.primary_assembly.genome.fa '
    genome_index =' /lustre/beagle2/xyang/Ivan/Carlos/Reference_Sequences/Human/FASTA/GRCh38.primary_assembly.genome '
    gensize='hs'


else:
    print("specify mm9, mm10, hg19, or hg38")


##############################################################################
##  make directories for output
##
##############################################################################


#parent = os.path.dirname(os.path.abspath(args.out))


parentList= glob.glob(os.path.join(args.fastqDir, '[!Undetermined]*.fastq.gz'))


mainseed = list(set([x.strip().split('/')[-1].strip().split('_')[0] for x in parentList]))


bowtie2Dir = parent+'/bowtie2_'+args.genome+'/'

if not os.path.exists(bowtie2Dir):
    os.makedirs(bowtie2Dir)

bowtie2Dir = os.path.abspath(bowtie2Dir)


macs2 = parent+'/macs2_'+args.genome+'/'

if not os.path.exists(macs2):
   	os.makedirs(macs2)

macs2 = os.path.abspath(macs2)



##############################################################################
##
##   make pbs and log directories
##
##############################################################################


if not os.path.exists(bowtie2Dir+'/log/'):
    os.makedirs(bowtie2Dir+'/log/')

bwPBS = bowtie2Dir+'/log/'

if not os.path.exists(macs2+'/log/'):
        os.makedirs(macs2+'/log/')

mcPBS = macs2+'/log/'



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
##  write bowtie2 PBS scripts
##  runs with default.
##############################################################################
   

for fs in parentList:
  if args.readType == "paired":
    f=fs.strip().split('/')[-1].strip().split('.')[0]
    pbs_writer(pbsDir = bwPBS, seed=f, processors=str(args.processors),
    walltime=str(args.walltime), shell = '/lustre/beagle2/cperez5/soft/shellScripts/bowtie2_runner_paired.sh  -o ' 
    +bowtie2Dir+ ' -g '+genome_index+ ' -p ' +str(args.processors)+ '  -s '+f+ ' -q 30 -d '+args.fastqDir)
  else:
    f=fs.strip().split('/')[-1].strip().split('_')[0]
    pbs_writer(pbsDir = bwPBS, seed=f, processors=str(args.processors),
    walltime=str(args.walltime), shell = '/lustre/beagle2/cperez5/soft/shellScripts/bowtie2_runner.sh  -o ' 
    +bowtie2Dir+ ' -g '+genome_index+ ' -p ' +str(args.processors)+ '  -s '+f+ ' -q 30 -d '+args.fastqDir)


##############################################################################
##	Write macs2 pbs scripts.  add any additional parameters 
##
##############################################################################


for fs in mainseed:
  if args.experiment == 'ATACseq':
    pbs_writer(pbsDir= mcPBS, seed='mc_'+fs, walltime=str(args.walltime), 
    processors=str(args.processors), shell =
    "/lustre/beagle2/cperez5/soft/shellScripts/macs2_runner.sh -o macs2 -t BAM -f "+bowtie2Dir+fs+".bam -g "+gensize+" -s "+fs+" -x \'--keep-dup --call-summit -nomodel --shift -100 --extsize 200\'")
  
  else:
    pbs_writer(pbsDir= mcPBS, seed='mc_'+fs, walltime=str(args.walltime), 
    processors=str(args.processors), shell =
    "/lustre/beagle2/cperez5/soft/shellScripts/macs2_runner.sh -o macs2 -t BAM -f "+bowtie2Dir+fs+".bam -g "+gensize+" -s "+fs+" -x \'--keep-dup --call-summit\'" )




######################################
##
##  Get mapping stats and plot with R
##
######################################



#pbs_writer(pbsDir = qcjob,  seed='bamsummary', walltime ='01', processors=str(args.processors), shell = 'python /lustre/beagle2/cperez5/soft/python2Scripts/samtools_bamsummary_local.py -d ' +cleaned_bam+ )

#pbs_writer(pbsDir = mcPBS, seed = 'mcsqSummry' , walltime= '01', processors= str(args.processors), shell='python /lustre/beagle2/cperez5/soft/python2Scripts/macs2_summary.py -d '+macs2)

#pbs_writer(pbsDir = mcPBS, seed = 'plotRline' , walltime= '01', processors= str(args.processors), shell='R CMD BATCH /lustre/beagle2/cperez5/soft/R_scripts/plot_line.py '+bwPBS+'tophSummry.pbs' + qcjob+' bamsummary.pbs '+mcPBS+'mcsqSummry.pbs ' )



##############################################################################
##
##  Run all qsub commands
##
##############################################################################

if os.path.isfile(args.out+"/submitALL.sh"):
  os.remove(args.out+"/submitALL.sh")
with open(args.out+'/submitALL.sh', "a") as runnit:
    runnit.write("#!/bin/bash\n \n \n")

for name in mainseed:
  bwcmd = 'bw'+name+'=$(qsub '+bwPBS+name+ '.pbs)\necho -e $bw'+name+ ' \n'
  mcscmd='mcs'+name.replace("-","")+'=$(qsub -W depend=afterok:$bw'+name +' '+mcPBS+'mc_'+name+ '.pbs) \n echo -e $mcs'+name.replace("-","") + ' \n'

  
  with open(args.out+'/submitALL.sh', "a") as runnit:
    runnit.write(bwcmd+ '\nsleep 2\n\n\n\n')
    runnit.write(mcscmd+ '\nsleep 2\n\n\n\n')
    runnit.close()
#for name in mainseed:
#  qccm='qc'+name.replace("-","")+'=$(qsub -W depend=afterok:'+ ' $bw+'f + ' '+qcjob+'qc_'+name+'.pbs) \n echo -e $qc'+name.replace("-","")+ ' \n'
#  summ2='bamsummary=$(qsub -W depend=afterok:+ $qc'+name.replace("-","") +' '+'bamsummary.pbs) \n echo -e $bamsummary.pbs\n'

#  with open(args.out+'/submitALL.sh', "a") as runnit:
#    runnit.write(qccm+ '\nsleep 2\n\n\n\n')
#    runnit.close()


os.system('chmod u+x ' +args.out+"/submitALL.sh")
#os.system( +args.out+"/submitALL.sh")

