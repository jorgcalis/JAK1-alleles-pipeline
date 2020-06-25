#!/usr/local/bin/python2.7

# Jorg Calis
# April 19, 2018
# Last update: August 14, 2018

# This is a script to examine the JAK1 library (aligned) reads from a JAK1+polyA inDrop scRNAseq sequencing experiment. Using polyA-capturing beads 
# modified to also include JAK1 primers. The JAK1 RT product was prepped into a separate library that was analyzed by the inDrop pipeline 
# (https://github.com/indrops/indrops). Per sequencing run, a alignment (bam) file was obtained for the JAK1 libraries. 
# Now, this script aims to analyze and combine the alignments from different sequencing runs on the same sample, to give a per cell count 
# of the number of wild-type and mutant JAK1 UMIs.
# In addition, the script will output stats on the number of reads/etc. As well as on non-JAK1 alignments in the library. 
# - inDrop version3 beads were used. 

# Samples:
# - For each sample, a set of bam files needs to be present without specific naming:
#   . Specified in the command line per sample in the following way:
#     - "sample,sampleName,bam1path,bam2path,bamNpath" 
#     - A space separates different sample entries
#     - Each sample entry is comma-separated and starts with the text "sample", followed by the name of the sample
#     - The 3rd and subsequent items in the sample entry are paths to the bam-files for this sample
#     - No special signs in the sampleName, only a-z A-Z 0-9 _
# - Assume the alignments of individual read are in the same file, and that the alignments are name sorted

# goals of this script:
# - Determine the per cell number of wildtype and mutant JAK1 UMIs

# Methods:
# - Analyze each mono-aligned read, get cell-barcode, UMI, and alignment info
#   . alignment info can be JAK1 wt/mutant, or non-JAK alignment
# - Analyze per cell
#   . Number of JAK1 wt UMIs, number of JAK1 mutant UMIs
#   . Consider read-counts, UMI similarity, etc in detailed pipeline
#     - Keep/discard UMIs with too few reads or too high similarity to high-read UMI in the same cell
# - Give alignment/workflow stats(+figures) 
#   . Number of reads mono/multi-aligned. Per bam 
#   . Per (top-10) gene number of mono aligned reads. All bams
#   . Per (top-10) gene number of mono aligned UMIs. All bams
#   . Number of mono aligned reads to JAK1-wt / JAK1-mutant / other . Per bam
#   . Number of mono aligned UMIs to JAK1-wt / JAK1-mutant / other . Per bam
#   . Histogram of number of reads per cell. (All bams)
#   . Histogram of number of UMIs per cell. (All bams)
#   . Histogram of number of reads per cell/UMI. (All bams)
#   . Histogram of number of reads per cell/UMI. (All bams, JAK1-wildtype UMIs)
#   . Histogram of number of reads per cell/UMI. (All bams, JAK1-mutant UMIs)
#   . Within versus between cells average UMI distances (based on all UMIs, also those with n=1 read)
# - Output 
#   . Stats as mentioned above
#   . Per cell count of the number of wild-type and mutant JAK1 UMIs.
#   . Within versus between cells average UMI distances (based on final transcripts, ie kept UMIs)

# Logging:
# - Stats as described above

### runmodes ###
# - EnsemblMode=GRCh37           : Use /data/genomes/Homo_sapiens/Ensembl/GRCh37/Annotation/Genes/genes.gtf
#                                  /data/genomes/Homo_sapiens/Ensembl/GRCh37/Sequence/WholeGenomeFasta/genome.fa
# - EnsemblMode=GRCh38 (default) : Use /data/genomes/Homo_sapiens/Ensembl/GRCh38/Annotation/Genes/Homo_sapiens.GRCh38.88.chr.gtf and 
#                                  /data/genomes/Homo_sapiens/Ensembl/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens.GRCh38.dna.chromosomes.fa
# - outputDirectory=<apath>        : Use directory <apath> for output
# - outputDirectory='./' (default) : Use current directory for output
################# Detailed description of pipeline to go from alignments to per cell JAK1 (wt and mutant) transcript counts #################


### loading modules
import os,time,sys,string,numpy,HTSeq,gzip,itertools
import matplotlib
matplotlib.use('Agg')
import pylab
pylab.rcParams['pdf.fonttype']=42


### runmodes ###
# - outputDirectory=<apath>        : Use directory <apath> for output
# - outputDirectory='./' (default) : Use current directory for output, if outputDirectory is unspecified

outputDirectory=None
for argument in sys.argv:
	if 'outputDirectory=' in argument:
		if not outputDirectory==None:
			raise Exception('outputDirectory set twice? '+string.join(sys.argv,' '))
		else:
			outputDirectory=argument.replace('outputDirectory=','')

if outputDirectory==None or outputDirectory=='False':
	# setting outputDirectory to default if not specified
	outputDirectory=os.getcwd()+'/'

if not outputDirectory[-1]=='/':raise Exception('outputDirectory needs to end in /')

# - EnsemblMode=GRCh37           : Use /data/genomes/Homo_sapiens/Ensembl/GRCh37/Annotation/Genes/genes.gtf
#                                  /data/genomes/Homo_sapiens/Ensembl/GRCh37/Sequence/WholeGenomeFasta/genome.fa
# - EnsemblMode=GRCh38 (default) : Use /data/genomes/Homo_sapiens/Ensembl/GRCh38/Annotation/Genes/Homo_sapiens.GRCh38.88.chr.gtf and 
#                                  /data/genomes/Homo_sapiens/Ensembl/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens.GRCh38.dna.chromosomes.fa
EnsemblMode=None
for argument in sys.argv:
	if 'EnsemblMode=' in argument:
		if not EnsemblMode==None:
			raise Exception(argv)
		else:
			EnsemblMode={'GRCh38':'GRCh38','GRCh37':'GRCh37'}[argument.replace('EnsemblMode=','')] # this will create an error if it is not GRCh37 or GRCh38

if EnsemblMode==None:
	# setting EnsemblMode to default if not specified
	EnsemblMode='GRCh38'


print 'EnsemblMode: %s; outputDirectory: %s'%(EnsemblMode,outputDirectory)

### Setting the JAK1 interval of interest, and wildtype and patient sequence at this interval.
WTmutJAKinterval=HTSeq.GenomicInterval('1',64845518,64845521,'-') # JAK1 position (64845519) mutated/wildtype plus/minus 1 nt
WTmutPositionNTs={'wt':{64845518:'G',64845519:'C',64845520:'T'},'mut':{64845518:'G',64845519:'A',64845520:'T'}}# Genomic (+ strand) nt's. not the actual wt/mutant sequences

### fetching dataset names and bams
# Samples:
# - For each sample, a set of bam files needs to be present without specific naming:
#   . Specified in the command line per sample in the following way:
#     - "sample,sampleName,bam1path,bam2path,bamNpath" 
#     - A space separates different sample entries
#     - Each sample entry is comma-separated and starts with the text "sample", followed by the name of the sample
#     - The 3rd and subsequent items in the sample entry are paths to the bam-files for this sample
#     - No special signs in the sampleName, only a-z A-Z 0-9 _
# - Assume the alignments of individual read are in the same file, and that the alignments are name sorted
def givePerSampleBAMs():
	perSampleBams={}
	for argument in sys.argv:
		if 'sample' in argument and argument.split(',')[0]=='sample':
			argumentSPLIT=argument.split(',')
			if not len(argumentSPLIT)>=3:raise Exception(argument+'not properly set')
			sampleIndicator,sampleName=argumentSPLIT[:2]
			BAMfiles=argumentSPLIT[2:]
			if not all([bamfile[-4:]=='.bam' for bamfile in BAMfiles]):raise Exception(argument+' not all bams?')
			if not len(BAMfiles)==len(set(BAMfiles)):raise Exception(argument+' repeated bams?')
			for BAMfile in BAMfiles:
				if not os.path.isfile(BAMfile):
					raise Exception(BAMfile+' does not exist?')
			if perSampleBams.has_key(sampleName):raise Exception(sampleName+' repeated sample name?')
			if not all([character in string.ascii_letters+string.digits+'_' for character in sampleName]):raise Exception(sampleName+' should be only a-z A-Z 0-9 _')
			outputFNs=['JAK1analysis_histogram_%s_JAK1locusReadsPerCell.pdf'%(sampleName),'JAK1analysis_histogram_%s_readPerCell.pdf'%(sampleName),
				'JAK1analysis_histogram_%s_readsPerCellUMI_JAKmut.pdf'%(sampleName),'JAK1analysis_histogram_%s_readsPerCellUMI_JAKwt.pdf'%(sampleName),
				'JAK1analysis_histogram_%s_readsPerCellUMI.pdf'%(sampleName),'JAK1analysis_histogram_%s_UMIsPerCell.pdf'%(sampleName),
				'JAK1analysis_JAK1transcriptCounts_perCell_%s.txt'%(sampleName),'JAK1analysis_stats_%s.txt'%(sampleName)]
			for outputFN in outputFNs:
				outputStatsFile=outputDirectory+'stats_%s.txt'%(sampleName)
				if os.path.isfile(outputDirectory+outputFN):
					raise Exception(outputDirectory+outputFN+' already exists')
			perSampleBams[sampleName]=BAMfiles
	return perSampleBams


perSampleBams=givePerSampleBAMs()
sampleNames=perSampleBams.keys()




### some output to stdout for early trouble shooting
for sampleName in sampleNames:
	print '''############################ background for %(sampleName)s ############################
%(time)s
Output directory: %(outputDirectory)s
Sample name: %(sampleName)s 
BAM(s): 
 - %(BAMline)s

'''%{'sampleName':sampleName,'time':time.ctime(),'outputDirectory':outputDirectory,'BAMline':string.join([BAMfile for BAMfile in perSampleBams[sampleName]],'\n - ')}


### loading basic info on the human genome annotation
# - Loading standard bc1/bc2-combi to bc+4-letter abreviation codes (as done in inDrops pipeline from https://github.com/indrops/indrops)

# based on https://github.com/indrops/indrops/blob/master/indrops.py, from function rev_comp
#___tbl = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}
___tbl_revComplementBases = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}
def rev_comp(seq):
	# based on https://github.com/indrops/indrops/blob/master/indrops.py, from function rev_comp
	#return ''.join(___tbl[s] for s in seq[::-1])
	return string.join([___tbl_revComplementBases[s] for s in seq[::-1]],'')

def giveNucleotideBasedToIndropBasedCellBCs():
	# based on https://github.com/indrops/indrops/blob/master/indrops.py, from function stable_barcode_names
	# BC2s from https://github.com/indrops/indrops/blob/master/ref/barcode_lists/gel_barcode2_list.txt
	bc2sFileContent = 'AAACAAAC\nAAACACGG\nAAACACTA\nAAACCGCC\nAAACGATC\nAAACGTGA\nAAACTACA\nAAACTGTG\nAAAGAAAG\nAAAGAGGC\nAAAGCCCG\nAAAGTCAT\nAAATAGCA\nAAATTCCG\nAACAAATG\nAACAGAAC\nAACAGCGG\nAACGATTT\nAACGCCAA\nAACGGTAG\nAACGTTAC\nAACTCAGT\nAACTGCCT\nAAGAACAG\nAAGAAGGT\nAAGAGTAT\nAAGCCTTC\nAAGCTCCT\nAAGGATGA\nAAGGCGCT\nAAGGGACC\nAAGTATTG\nAAGTCCAA\nAAGTCGGG\nAAGTGAGA\nAAGTTGTC\nAATAAGGA\nAATACATC\nAATATGAC\nAATCCGGC\nAATCGAAG\nAATCGTTC\nAATGGCGT\nAATGTATG\nACAAAGAT\nACAAGTAG\nACAATCTT\nACACCAAG\nACAGATAA\nACAGGCCA\nACATCTCG\nACATGGAC\nACCAACCC\nACCAAGGG\nACCACAGA\nACCAGTTT\nACCCATGC\nACCCGATT\nACCCTCAA\nACCGTCGA\nACCTGAAG\nACCTTCCC\nACGAATTC\nACGACGAC\nACGCTTAA\nACGGAGCA\nACGGCAGT\nACGGGTTA\nACGGTTGG\nACGTAAAC\nACTAATTG\nACTACCCG\nACTAGAGC\nACTCATAC\nACTCGGAA\nACTGCTGG\nACTGGTCA\nACTTCGCT\nAGAAACCA\nAGAAAGTG\nAGAAGCTT\nAGAATCAA\nAGACCTCA\nAGACGAGG\nAGAGAGAC\nAGAGGTGC\nAGCAACGC\nAGCACGTA\nAGCATGCC\nAGCCATCT\nAGCGTGGT\nAGCTCCAC\nAGCTTCGA\nAGGACACA\nAGGAGTCG\nAGGCAATA\nAGGCCGAA\nAGGCGTTT\nAGGGACTG\nAGGGTAAA\nAGGTAAGC\nAGGTATAT\nAGGTTCCC\nAGTAATGG\nAGTAGTTA\nAGTCACAA\nAGTCCGTG\nAGTGCTTC\nAGTTGAAC\nAGTTGCGG\nAGTTTGTA\nATAACAGG\nATAAGCTA\nATACACCC\nATACTCTC\nATAGATGT\nATATGCAA\nATATGGGT\nATCAATCG\nATCAGGGA\nATCCCACC\nATCCGCAT\nATCCTAGT\nATCGCGCT\nATCGTAAC\nATCTTGGC\nATGACAAC\nATGACTTG\nATGCATAT\nATGCGGAG\nATGGGCTC\nATGGTCTG\nATGTGCCG\nATTACCTT\nATTATTCG\nATTCTGAG\nATTGAAGT\nATTGGCCC\nATTTCCAT\nATTTGTTG\nCAAACATT\nCAACGCAG\nCAAGGAAT\nCAAGGGTT\nCAAGGTAC\nCAATCTAG\nCAATTCTC\nCACAACCT\nCACAAGTA\nCACTAACC\nCACTTGAT\nCAGACTCG\nCAGATGGG\nCAGGTTGC\nCAGTTTAA\nCATGACGA\nCATGCTGC\nCATTCATT\nCATTCGGG\nCATTTCTA\nCCACCTCT\nCCACGTTG\nCCAGACAG\nCCAGCGAA\nCCATATGA\nCCATCCAC\nCCATCGTC\nCCATGCAT\nCCCGTAAG\nCCCGTTCT\nCCCTCTTG\nCCCTGTTT\nCCCTTGCA\nCCGACTTT\nCCGAGATC\nCCGATACG\nCCGGAAAT\nCCGTAGCT\nCCGTCTTA\nCCTACGCT\nCCTATTTA\nCCTCATGA\nCCTTTACA\nCCTTTGTC\nCGAAACTC\nCGAACCGA\nCGAAGAAG\nCGACATTT\nCGAGGCTA\nCGATCCAA\nCGATGGCA\nCGGACTAA\nCGGCTGTA\nCGGTGAGT\nCGTACCGA\nCGTCGAAT\nCGTGCAAC\nCGTGGGAT\nCGTGTACA\nCGTGTGTT\nCGTTGCCT\nCGTTTCGT\nCTAACGCC\nCTACGGGA\nCTAGACTA\nCTAGCACG\nCTAGTAGG\nCTCAAACA\nCTCACATC\nCTCCCAAA\nCTCCTCCA\nCTCGGTGA\nCTCTATAG\nCTCTGCGT\nCTGAAGGG\nCTGAGCGT\nCTGCGATG\nCTGCTAGA\nCTGGAACA\nCTGGGTAT\nCTGTCGCA\nCTGTGACC\nCTGTTAAA\nCTGTTGTG\nCTGTTTCC\nCTTAGGCC\nCTTAGTGT\nCTTCTACG\nCTTTATCC\nCTTTCACT\nCTTTGGAC\nGAAAGACA\nGAAATACG\nGAAGATAT\nGAATCCCA\nGAATGCGC\nGACACAAA\nGACACCTG\nGACTAGCG\nGAGAAACC\nGAGCGGAA\nGAGGAGTG\nGAGGGTCA\nGAGTGTAC\nGATACGCA\nGATGCAGA\nGATGGTTA\nGATGTGGC\nGATTAAAG\nGATTACTT\nGATTGGGA\nGATTTCCC\nGCAAACTG\nGCACTCAG\nGCATCACT\nGCATCGAG\nGCCAAAGC\nGCCAACAT\nGCCTGGTA\nGCCTTGTG\nGCGCTGAT\nGCGGTAAC\nGCGTATTC\nGCGTGCAA\nGCTAAGTT\nGCTACCGT\nGCTATGGG\nGCTCGTAG\nGCTTCTCC\nGGAACGAA\nGGAAGTCC\nGGACTGGA\nGGACTTCT\nGGAGGTTT\nGGAGTAAG\nGGATTGTT\nGGCAAGGT\nGGCACTTC\nGGCCCAAT\nGGCGACAA\nGGCTATAA\nGGCTTTGC\nGGGAGATG\nGGGATTAC\nGGGCATCA\nGGGTCATT\nGGGTCTAG\nGGTAAATC\nGGTAGCCA\nGGTCCTAA\nGGTCTTTC\nGGTGTCGA\nGGTTACAC\nGGTTAGGG\nGGTTGAGA\nGTAAACAA\nGTAAGCCG\nGTAATCTG\nGTACGCTT\nGTACGGAC\nGTATACGT\nGTATTGAC\nGTCAAGAG\nGTCAGACC\nGTCAGGTT\nGTCCACTA\nGTCCGTCA\nGTCCTTGC\nGTCTAATC\nGTCTGGAA\nGTCTTCCT\nGTGAACTC\nGTGAGGCA\nGTGATAAA\nGTGCCCAT\nGTGCGAAG\nGTGGTGCT\nGTGTCACC\nGTGTCAGG\nGTTACTAG\nGTTCTGCT\nGTTGTCCG\nTAAACCGA\nTAACTTCT\nTAAGGGCC\nTAATCCAT\nTAATGTGG\nTACCCTGC\nTACCGCTC\nTACCTAAG\nTACCTCCC\nTACGCGAG\nTACGTTCG\nTACTGAAT\nTAGATCAA\nTAGCCACA\nTAGCGGAT\nTAGGCTTT\nTAGGTACG\nTAGTAGCC\nTAGTCTCT\nTATCCACG\nTATCTGTC\nTATGTGAA\nTATTAGCG\nTCAAATGG\nTCAAGGCG\nTCAGCCTC\nTCATACCA\nTCATAGCT\nTCATTTCA\nTCCAGAAG\nTCCCTGGA\nTCCGACAC\nTCCGCTGT\nTCCTATAT\nTCGACTGC\nTCGAGTTT\nTCGCAATC\nTCGGTCAT\nTCGTGGGT\nTCGTTCCC\nTCTAAACT\nTCTATTCC\nTCTGATTT\nTCTTTGAC\nTGAATAGG\nTGAATCCT\nTGACGTCG\nTGAGAGCG\nTGAGCACA\nTGCACCAG\nTGCCGGTA\nTGCGACTA\nTGCTGACG\nTGCTTCAT\nTGCTTGGG\nTGGAAAGC\nTGGACGGA\nTGGCTAGT\nTGGGAATT\nTGGTGTCT\nTGGTTAAC\nTGTTATCA'
	bc2s = [bc2 for bc2 in bc2sFileContent.split('\n')]
	rev_bc2s = [rev_comp(bc2) for bc2 in bc2s]
	v3_names = {}
	#barcode_iter = product(bc2s, rev_bc2s)
	barcode_iter = itertools.product(bc2s, rev_bc2s)
	#name_iter = product(string.ascii_uppercase, repeat=4)
	name_iter = itertools.product(string.ascii_uppercase, repeat=4)
	for barcode, name in zip(barcode_iter, name_iter):
		#v3_names['-'.join(barcode)] = 'bc' + ''.join(name)
		v3_names[string.join(barcode,'-')] = 'bc' + string.join(name,'')
	return v3_names


nucleotideBasedToIndropBasedCellBCs=giveNucleotideBasedToIndropBasedCellBCs()

# - Loading the gene info from the gtf in global objects
# - EnsemblMode=GRCh37           : Use /data/genomes/Homo_sapiens/Ensembl/GRCh37/Annotation/Genes/genes.gtf
#                                  /data/genomes/Homo_sapiens/Ensembl/GRCh37/Sequence/WholeGenomeFasta/genome.fa
# - EnsemblMode=GRCh38 (default) : Use /data/genomes/Homo_sapiens/Ensembl/GRCh38/Annotation/Genes/Homo_sapiens.GRCh38.88.chr.gtf and 
#                                  /data/genomes/Homo_sapiens/Ensembl/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens.GRCh38.dna.chromosomes.fa

if EnsemblMode=='GRCh37':
	GTFfileName='/data/genomes/Homo_sapiens/Ensembl/GRCh37/Annotation/Genes/genes.gtf'
	genomePath='/data/genomes/Homo_sapiens/Ensembl/GRCh37/Sequence/WholeGenomeFasta/genome.fa'

if EnsemblMode=='GRCh38':
	GTFfileName='/data/genomes/Homo_sapiens/Ensembl/GRCh38/Annotation/Genes/Homo_sapiens.GRCh38.88.chr.gtf'
	genomePath='/data/genomes/Homo_sapiens/Ensembl/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens.GRCh38.dna.chromosomes.fa'

print 'Running in EnsemblMode %s:\n - GTFfileName: %s\n - genomePath (not used in this analysis): %s'%(EnsemblMode,GTFfileName,genomePath)


def SetUpGlobalObjectsForHumanOnlyMapping():
	global gtffile,AllFeats,StrandSwitcher,ENSTtoENSG,ENSGtoGeneName
	#gtffile = HTSeq.GFF_Reader("/data/genomes/Homo_sapiens/Ensembl/GRCh37/Annotation/Genes/genes.gtf")
	gtffile = HTSeq.GFF_Reader(GTFfileName)
	AllFeats = HTSeq.GenomicArrayOfSets( "auto", stranded=True )
	n=0
	ENSTtoENSG={}
	ENSGtoGeneName={}
	for feature in gtffile:
		n+=1
		if n%25==0:
			print n,'\r',
			sys.stdout.flush()
		AllFeats[ feature.iv ] += feature
		if feature.attr.has_key('gene_id'):
			if ENSGtoGeneName.has_key(feature.attr['gene_id']) and ENSGtoGeneName[feature.attr['gene_id']]!=feature.attr.get('gene_name',None):
				print
				print ENSGtoGeneName[feature.attr['gene_id']],feature.attr.get('gene_name',None)
				print
			ENSGtoGeneName[feature.attr['gene_id']]=feature.attr.get('gene_name',None)
			if feature.attr.has_key('transcript_id'):
				if ENSTtoENSG.has_key(feature.attr['transcript_id']) and ENSTtoENSG[feature.attr['transcript_id']]!=feature.attr['gene_id']:
					print
					print feature.attr['transcript_id'],feature.attr['gene_id']
					print
				ENSTtoENSG[feature.attr['transcript_id']]=feature.attr['gene_id']
	StrandSwitcher={'-':'+','+':'-'}
	print



### Functions for determining cell-UMI-barcodes information
# important standard functions
def giveReverseComplement(DNAseq):
    #convert all to capitals, than rev/compl, non-ATCG stays non-ATCG
    DNAseq=DNAseq.upper()
    DNAseq=DNAseq.replace('A','x')
    DNAseq=DNAseq.replace('T','A')
    DNAseq=DNAseq.replace('x','T')
    DNAseq=DNAseq.replace('C','x')
    DNAseq=DNAseq.replace('G','C')    
    DNAseq=DNAseq.replace('x','G')
    return DNAseq[::-1]

# ------------------------- Analyze alignments per sample ------------------------- #
# - Analyze each mono-aligned read, get cell-barcode, UMI, and alignment info
#   . alignment info can be JAK1 wt/mutant, or non-JAK alignment
# - Analyze per cell
#   . Number of JAK1 wt UMIs, number of JAK1 mutant UMIs
#   . Consider read-counts, UMI similarity, etc in detailed pipeline
# - Give alignment/workflow stats(+figures)
#   . Number of reads all/mono/multi-aligned. Per bam 
#   . Number of mono aligned reads to JAK1-wt / JAK1-mutant / JAK1-undefined / other . Per bam
#   . Number of mono aligned UMIs to JAK1-wt / JAK1-mutant / JAK1-undefined / other . All bams
#   . Per (top-10) gene number of mono aligned reads. All bams
#   . Per (top-10) gene number of mono aligned UMIs. All bams
#   . Histogram of number of reads per cell. (All bams)
#   . Histogram of number of UMIs per cell. (All bams)
#   . Histogram of number of reads per cell/UMI. (All bams)
#   . Histogram of number of reads per cell/UMI. (All bams, JAK1-wildtype UMIs)
#   . Histogram of number of reads per cell/UMI. (All bams, JAK1-mutant UMIs)
#   . Within versus between cells average UMI distances (based on all UMIs, also those with n=1 read)
#   . 
# - Output 
#   . Stats as mentioned above
#   . Per cell count of the number of wild-type and mutant JAK1 UMIs.
#   . Within versus between cells average UMI distances (based on final transcripts, ie kept UMIs)


def writeDownANDgivePerCellJAK1wildtypeAndMutantTranscriptCounts(perCellPerUMIinfo=None,sampleName=None,minReadCNTperUMI=None,minWTvsMUTreadCNTdifference=None):
	# getting per cell wt/mutant JAK1 transcript counts
	perCellJAK1transcriptCounts={}
	# - If cell/UMI has multiple assignment from different reads, pick the one with most reads for it. 
	# - only pick if at least minReadCNTperUMI read count
	for cellBC in perCellPerUMIinfo.keys():
		perCellJAK1transcriptCounts[cellBC]={'wt':0,'mutant':0}
		for UMI in perCellPerUMIinfo[cellBC].keys():
			# for "JAK1wtMonoCNT","JAK1differentSiteMonoCNT","JAK1mutMonoCNT","JAK1undefinedMonoCNT", other gene, select the one with most reads
			perQualificationReadCounts=[]
			for item in ["JAK1wtMonoCNT","JAK1differentSiteMonoCNT","JAK1mutMonoCNT","JAK1undefinedMonoCNT"]:
				perQualificationReadCounts.append((perCellPerUMIinfo[cellBC][UMI][item],item))
			for gn in perCellPerUMIinfo[cellBC][UMI]['nonJAK1MonoCNTperGene']:
				perQualificationReadCounts.append((perCellPerUMIinfo[cellBC][UMI]['nonJAK1MonoCNTperGene'][gn],gn))
			perQualificationReadCounts.sort(reverse=True)
			WTplusMUTreadTotal=sum([cnt for cnt,qualif in perQualificationReadCounts if qualif in ["JAK1wtMonoCNT","JAK1mutMonoCNT","JAK1undefinedMonoCNT"]])
			if (perQualificationReadCounts[0][1] in ["JAK1wtMonoCNT","JAK1mutMonoCNT","JAK1undefinedMonoCNT"]) and (not perQualificationReadCounts[0][0]>=(minWTvsMUTreadCNTdifference*WTplusMUTreadTotal)):
				continue
			if perQualificationReadCounts[0][0]>perQualificationReadCounts[1][0] and perQualificationReadCounts[0][0]>=minReadCNTperUMI: # one qualification with heighest read count, that is at least minReadCNTperUMI
				umiQualification=perQualificationReadCounts[0][1]
				if umiQualification=="JAK1wtMonoCNT":
					perCellJAK1transcriptCounts[cellBC]['wt']+=1
				if umiQualification=="JAK1mutMonoCNT":
					perCellJAK1transcriptCounts[cellBC]['mutant']+=1
	# writing per cell wt/mutant JAK1 transcript counts
	OUTfile=open(outputDirectory+'JAK1analysis_JAK1transcriptCounts_perCell_%s.txt'%(sampleName),'w')
	OUTfile.write("cellBC\tJAK1wtTXcnt\tJAK1mutTXcnt\n")
	for cellBC in sorted(perCellJAK1transcriptCounts.keys()):
		JAK1wtTXcnt=perCellJAK1transcriptCounts[cellBC]['wt']
		JAK1mutTXcnt=perCellJAK1transcriptCounts[cellBC]['mutant']
		OUTfile.write("%(cellBC)s\t%(JAK1wtTXcnt)s\t%(JAK1mutTXcnt)s\n"%{'cellBC':cellBC,'JAK1wtTXcnt':JAK1wtTXcnt,'JAK1mutTXcnt':JAK1mutTXcnt,})
	OUTfile.close()
	return perCellJAK1transcriptCounts


def writeDownTXidentificationStats(perCellJAK1transcriptCounts=None,logFile=None,minReadCNTperUMI=None,minWTvsMUTreadCNTdifference=None):
	logTEXT='''#### Final transcript identification stats based on all bam files ####
# UMI assigned to either JAK1 wt/mut/undefined/otherLocus or other gene based on read count
# UMI assigned if at least %(minReadCNTperUMI)s reads for this assignment
# UMI assigned to JAK1-wt or JAK1-mutant if at least %(minWTvsMUTreadCNTdifference)s fraction of JAK1-wt+JAK1-mutant reads are assigned to this UMI
'''%{'minReadCNTperUMI':minReadCNTperUMI,'minWTvsMUTreadCNTdifference':minWTvsMUTreadCNTdifference}
	for i in range(10):
		nWT=len([cl for cl in perCellJAK1transcriptCounts.keys() if perCellJAK1transcriptCounts[cl]['wt']==i])
		nMUT=len([cl for cl in perCellJAK1transcriptCounts.keys() if perCellJAK1transcriptCounts[cl]['mutant']==i])
		logTEXT+='\n - %s cells with %s JAK1-wt transcripts'%(nWT,i)
		logTEXT+='\n - %s cells with %s JAK1-mutant transcripts'%(nMUT,i)
	logTEXT+='\n - %s cells with 10+ JAK1-wt transcripts'%(len([cl for cl in perCellJAK1transcriptCounts.keys() if perCellJAK1transcriptCounts[cl]['wt']>=10]))
	logTEXT+='\n - %s cells with 10+ JAK1-mutant transcripts'%(len([cl for cl in perCellJAK1transcriptCounts.keys() if perCellJAK1transcriptCounts[cl]['mutant']>=10]))
	nWTonly=len([cl for cl in perCellJAK1transcriptCounts.keys() if perCellJAK1transcriptCounts[cl]['wt']>0 and perCellJAK1transcriptCounts[cl]['mutant']==0])
	nMUTonly=len([cl for cl in perCellJAK1transcriptCounts.keys() if perCellJAK1transcriptCounts[cl]['wt']==0 and perCellJAK1transcriptCounts[cl]['mutant']>0])
	nBOTH=len([cl for cl in perCellJAK1transcriptCounts.keys() if perCellJAK1transcriptCounts[cl]['wt']>0 and perCellJAK1transcriptCounts[cl]['mutant']>0])
	nNONE=len([cl for cl in perCellJAK1transcriptCounts.keys() if perCellJAK1transcriptCounts[cl]['wt']==0 and perCellJAK1transcriptCounts[cl]['mutant']==0])
	logTEXT+='\n\n %s cells with only JAK1-wt transcripts'%(nWTonly)
	logTEXT+='\n %s cells with only JAK1-mutant transcripts'%(nMUTonly)
	logTEXT+='\n %s cells both JAK1-mutant and JAK1-wt transcripts'%(nBOTH)
	logTEXT+='\n %s cells with no JAK1-mutant/wt transcripts\n\n'%(nNONE)
	print logTEXT
	logFile.write(logTEXT)


def makeAllBAMbasedHistograms_andMore(logFile=None,perCellPerUMIinfo=None,sampleName=None):
	#   . Histogram of number of reads per cell. (All bams)
	#   . Histogram of number of UMIs per cell. (All bams)
	#   . Histogram of number of reads per cell/UMI. (All bams)
	#   . Histogram of number of reads per cell/UMI. (All bams, JAK1-wildtype UMIs)
	#   . Histogram of number of reads per cell/UMI. (All bams, JAK1-mutant UMIs)
	# - Print number of cells with more than 1/10/100/1000/10000 reads
	readPerCell=[sum([perCellPerUMIinfo[cellBC][UMI]['monoCNT'] for UMI in perCellPerUMIinfo[cellBC].keys()]) for cellBC in perCellPerUMIinfo.keys()]
	JAK1locusReadsPerCell=[sum([perCellPerUMIinfo[cellBC][UMI]['JAK1wtMonoCNT']+perCellPerUMIinfo[cellBC][UMI]['JAK1mutMonoCNT']+perCellPerUMIinfo[cellBC][UMI]['JAK1undefinedMonoCNT'] for UMI in perCellPerUMIinfo[cellBC].keys()]) for cellBC in perCellPerUMIinfo.keys()]
	UMIsPerCell=[len(perCellPerUMIinfo[cellBC].keys()) for cellBC in perCellPerUMIinfo.keys()]
	readsPerCellUMI=[perCellPerUMIinfo[cellBC][UMI]['monoCNT'] for cellBC in perCellPerUMIinfo.keys() for UMI in perCellPerUMIinfo[cellBC].keys()]
	readsPerCellUMI_JAKwt=[perCellPerUMIinfo[cellBC][UMI]['JAK1wtMonoCNT'] for cellBC in perCellPerUMIinfo.keys() for UMI in perCellPerUMIinfo[cellBC].keys()]
	readsPerCellUMI_JAKmut=[perCellPerUMIinfo[cellBC][UMI]['JAK1mutMonoCNT'] for cellBC in perCellPerUMIinfo.keys() for UMI in perCellPerUMIinfo[cellBC].keys()]
	for vals,analysisName in [(JAK1locusReadsPerCell,'JAK1locusReadsPerCell'),(readPerCell,'readPerCell'),(UMIsPerCell,'UMIsPerCell'),(readsPerCellUMI,'readsPerCellUMI'),(readsPerCellUMI_JAKwt,'readsPerCellUMI_JAKwt'),(readsPerCellUMI_JAKmut,'readsPerCellUMI_JAKmut')]:
		pylab.close()
		pylab.hist([numpy.log10(val) for val in vals if val>0],bins=250)
		pylab.savefig(outputDirectory+'JAK1analysis_histogram_%s_%s.pdf'%(sampleName,analysisName))
	print 'scp jcalis@heracles.rockefeller.edu:%s/JAK1analysis_histogram_*pdf .;open JAK1analysis_histogram_*pdf'%(os.getcwd())
	logTEXT='#### Cell identification stats based on all bam files ####'
	for NreadLimit in [1,10,100,1000,10000]:
		Ncells=len([cellBC for cellBC in perCellPerUMIinfo.keys() if 
			sum([perCellPerUMIinfo[cellBC][UMI]['monoCNT'] for UMI in perCellPerUMIinfo[cellBC].keys()])>=NreadLimit])
		Nreads=sum([sum([perCellPerUMIinfo[cellBC][UMI]['monoCNT'] for UMI in perCellPerUMIinfo[cellBC].keys()]) for cellBC in perCellPerUMIinfo.keys() if 
			sum([perCellPerUMIinfo[cellBC][UMI]['monoCNT'] for UMI in perCellPerUMIinfo[cellBC].keys()])>=NreadLimit])
		logTEXT+='\n - %s cells with at least %s any reads (%s reads)'%(Ncells,NreadLimit,Nreads)
	for NreadLimit in [1,10,100,1000,10000]:
		Ncells=len([cellBC for cellBC in perCellPerUMIinfo.keys() if 
			sum([perCellPerUMIinfo[cellBC][UMI]['JAK1wtMonoCNT']+perCellPerUMIinfo[cellBC][UMI]['JAK1mutMonoCNT']+perCellPerUMIinfo[cellBC][UMI]['JAK1undefinedMonoCNT'] 
			for UMI in perCellPerUMIinfo[cellBC].keys()])>=NreadLimit])
		Nreads=sum([sum([perCellPerUMIinfo[cellBC][UMI]['JAK1wtMonoCNT']+perCellPerUMIinfo[cellBC][UMI]['JAK1mutMonoCNT']+perCellPerUMIinfo[cellBC][UMI]['JAK1undefinedMonoCNT'] for UMI in perCellPerUMIinfo[cellBC].keys()]) for cellBC in perCellPerUMIinfo.keys() if 
			sum([perCellPerUMIinfo[cellBC][UMI]['JAK1wtMonoCNT']+perCellPerUMIinfo[cellBC][UMI]['JAK1mutMonoCNT']+perCellPerUMIinfo[cellBC][UMI]['JAK1undefinedMonoCNT'] 
			for UMI in perCellPerUMIinfo[cellBC].keys()])>=NreadLimit])
		logTEXT+='\n - %s cells with at least %s JAK1locus reads (%s reads)'%(Ncells,NreadLimit,Nreads)
	print logTEXT
	logFile.write(logTEXT)



def writeDownAllBAMsStats(perCellPerUMIinfo=None,logFile=None,minReadCNTperUMI=None,minWTvsMUTreadCNTdifference=None):
	#   . Number of mono aligned UMIs to JAK1-wt / JAK1-mutant / JAK1-undefined / other . All bams
	#   . Per (top-10) gene number of mono aligned UMIs. All bams
	#   . Within versus between cells average UMI distances (based on all UMIs, also those with n=1 read)
	# - If cell/UMI has multiple assignment from different reads, pick the one with most reads for it. 
	UMIcntSTATs={}
	for item in ["monoCNT","multiCNT","nAlignedReads","JAK1wtMonoCNT","JAK1differentSiteMonoCNT","JAK1mutMonoCNT","JAK1undefinedMonoCNT"]:
		UMIcntSTATs[item]=0 
	UMIcntSTATs['nonJAK1MonoCNTperGene']={}
	for cellBC in perCellPerUMIinfo.keys():
		for UMI in perCellPerUMIinfo[cellBC].keys():
			if perCellPerUMIinfo[cellBC][UMI]['monoCNT']>=minReadCNTperUMI:UMIcntSTATs['monoCNT']+=1
			if perCellPerUMIinfo[cellBC][UMI]['multiCNT']>=minReadCNTperUMI:UMIcntSTATs['multiCNT']+=1
			if perCellPerUMIinfo[cellBC][UMI]['nAlignedReads']>=minReadCNTperUMI:UMIcntSTATs['nAlignedReads']+=1
			# for "JAK1wtMonoCNT","JAK1differentSiteMonoCNT","JAK1mutMonoCNT","JAK1undefinedMonoCNT", other gene, select the one with most reads
			perQualificationReadCounts=[]
			for item in ["JAK1wtMonoCNT","JAK1differentSiteMonoCNT","JAK1mutMonoCNT","JAK1undefinedMonoCNT"]:
				perQualificationReadCounts.append((perCellPerUMIinfo[cellBC][UMI][item],item))
			for gn in perCellPerUMIinfo[cellBC][UMI]['nonJAK1MonoCNTperGene']:
				perQualificationReadCounts.append((perCellPerUMIinfo[cellBC][UMI]['nonJAK1MonoCNTperGene'][gn],gn))
			perQualificationReadCounts.sort(reverse=True)
			WTplusMUTreadTotal=sum([cnt for cnt,qualif in perQualificationReadCounts if qualif in ["JAK1wtMonoCNT","JAK1mutMonoCNT","JAK1undefinedMonoCNT"]])
			if (perQualificationReadCounts[0][1] in ["JAK1wtMonoCNT","JAK1mutMonoCNT","JAK1undefinedMonoCNT"]) and (not perQualificationReadCounts[0][0]>=(minWTvsMUTreadCNTdifference*WTplusMUTreadTotal)):
				continue
			if perQualificationReadCounts[0][0]>perQualificationReadCounts[1][0] and perQualificationReadCounts[0][0]>=minReadCNTperUMI: # one qualification with heighest read count, that is at least minReadCNTperUMI
				umiQualification=perQualificationReadCounts[0][1]
				if umiQualification in ["JAK1wtMonoCNT","JAK1differentSiteMonoCNT","JAK1mutMonoCNT","JAK1undefinedMonoCNT"]:
					UMIcntSTATs[umiQualification]+=1
				else:
					if not UMIcntSTATs['nonJAK1MonoCNTperGene'].has_key(umiQualification):
						UMIcntSTATs['nonJAK1MonoCNTperGene'][umiQualification]=0
					UMIcntSTATs['nonJAK1MonoCNTperGene'][umiQualification]+=1
	nAlignedReads,monoCNT,multiCNT=UMIcntSTATs['nAlignedReads'],UMIcntSTATs['monoCNT'],UMIcntSTATs['multiCNT']
	JAK1differentSiteMonoCNT,JAK1undefinedMonoCNT=UMIcntSTATs['JAK1differentSiteMonoCNT'],UMIcntSTATs['JAK1undefinedMonoCNT']
	JAK1mutMonoCNT,JAK1wtMonoCNT=UMIcntSTATs['JAK1mutMonoCNT'],UMIcntSTATs['JAK1wtMonoCNT']
	nJAKs=JAK1differentSiteMonoCNT+JAK1undefinedMonoCNT+JAK1mutMonoCNT+JAK1wtMonoCNT
	nJAKwtMutlocus=JAK1mutMonoCNT+JAK1wtMonoCNT+JAK1undefinedMonoCNT
	nonJAK1cnt=sum(UMIcntSTATs['nonJAK1MonoCNTperGene'].values())
	cntPerGN=[(cnt,gn) for gn,cnt in UMIcntSTATs['nonJAK1MonoCNTperGene'].items()]
	cntPerGN.sort(reverse=True)
	top10line=string.join([' - %s %s'%(cnt,gn) for cnt,gn in cntPerGN[:10]],'\n')
	withinUMIdists,btwnUMIdists=[],[]
	randomOrderCellBCs=perCellPerUMIinfo.keys()
	numpy.random.shuffle(randomOrderCellBCs)
	randomOrderCellBCsI=0
	for cellBC in perCellPerUMIinfo.keys()[:10]:
		UMIs=perCellPerUMIinfo[cellBC].keys()
		for i in range(len(UMIs)-1):
			for j in range(i+1,len(UMIs)):
				dist=sum([UMIs[i][p]!=UMIs[j][p] for p in range(len(UMIs[i]))])
				withinUMIdists.append(dist)
		otherCellUMIs=perCellPerUMIinfo[randomOrderCellBCs[randomOrderCellBCsI]].keys()
		randomOrderCellBCsI+=1
		for UMI1 in UMIs:
			for UMI2 in otherCellUMIs:
				dist=sum([UMI1[p]!=UMI2[p] for p in range(len(UMI1))])
				btwnUMIdists.append(dist)
	UMIdistLine='''UMI-to-UMI distance
within cells average/median: %.2f / %s
between cells average/median: %.2f / %s
'''%(numpy.average(withinUMIdists),numpy.median(withinUMIdists),numpy.average(btwnUMIdists),numpy.median(btwnUMIdists))
	logTEXT='''
#### UMI stats based on all bam files ####
# UMI assigned to either JAK1 wt/mut/undefined/otherLocus or other gene based on read count
# UMI assigned if at least %(minReadCNTperUMI)s reads for this assignment
# UMI assigned to JAK1-wt or JAK1-mutant if at least %(minWTvsMUTreadCNTdifference)s fraction of JAK1-wt+JAK1-mutant reads are assigned to this UMI
%(nAlignedReads)s aligned UMIs
%(monoCNT)s mono-mapped aligned UMIs
%(multiCNT)s multi-mapped aligned UMIs

%(nJAKs)s JAK1 aligned UMIs
- %(nJAKwtMutlocus)s aligned to JAK1-wt/mutant locus
  . %(JAK1wtMonoCNT)s wild-type
  . %(JAK1mutMonoCNT)s mutant
  . %(JAK1undefinedMonoCNT)s undefined
- %(JAK1differentSiteMonoCNT)s aligned outside the JAK1-wt/mutant locus

%(nonJAK1cnt)s non-JAK1 aligned UMIs (showing top-10)
%(top10line)s

%(UMIdistLine)s
'''%{'nAlignedReads':nAlignedReads,'monoCNT':monoCNT,'multiCNT':multiCNT,'nJAKs':nJAKs,
'JAK1differentSiteMonoCNT':JAK1differentSiteMonoCNT,'JAK1undefinedMonoCNT':JAK1undefinedMonoCNT,'JAK1mutMonoCNT':JAK1mutMonoCNT,
'JAK1wtMonoCNT':JAK1wtMonoCNT,'nonJAK1cnt':nonJAK1cnt,'top10line':top10line,'nJAKwtMutlocus':nJAKwtMutlocus,
'UMIdistLine':UMIdistLine,'minReadCNTperUMI':minReadCNTperUMI,'minWTvsMUTreadCNTdifference':minWTvsMUTreadCNTdifference}
	print logTEXT
	logFile.write(logTEXT)


def writeDownPerBAMstats(perBAMstats=None,logFile=None):
	logTEXT=''
	for bamfile in sorted(perBAMstats.keys()):
		nAlignedReads,monoCNT,multiCNT=perBAMstats[bamfile]['nAlignedReads'],perBAMstats[bamfile]['monoCNT'],perBAMstats[bamfile]['multiCNT']
		JAK1differentSiteMonoCNT,JAK1undefinedMonoCNT=perBAMstats[bamfile]['JAK1differentSiteMonoCNT'],perBAMstats[bamfile]['JAK1undefinedMonoCNT']
		JAK1mutMonoCNT,JAK1wtMonoCNT=perBAMstats[bamfile]['JAK1mutMonoCNT'],perBAMstats[bamfile]['JAK1wtMonoCNT']
		nJAKs=JAK1differentSiteMonoCNT+JAK1undefinedMonoCNT+JAK1mutMonoCNT+JAK1wtMonoCNT
		nJAKwtMutlocus=JAK1mutMonoCNT+JAK1wtMonoCNT+JAK1undefinedMonoCNT
		nonJAK1cnt=sum(perBAMstats[bamfile]['nonJAK1MonoCNTperGene'].values())
		cntPerGN=[(cnt,gn) for gn,cnt in perBAMstats[bamfile]['nonJAK1MonoCNTperGene'].items()]
		cntPerGN.sort(reverse=True)
		top10line=string.join([' - %s %s'%(cnt,gn) for cnt,gn in cntPerGN[:10]],'\n')
		logTEXT+='''
#### %(bamfile)s ####
%(nAlignedReads)s aligned reads
%(monoCNT)s mono-mapped aligned reads
%(multiCNT)s multi-mapped aligned reads

%(nJAKs)s JAK1 aligned reads
- %(nJAKwtMutlocus)s aligned to JAK1-wt/mutant locus
  . %(JAK1wtMonoCNT)s wild-type
  . %(JAK1mutMonoCNT)s mutant
  . %(JAK1undefinedMonoCNT)s undefined
- %(JAK1differentSiteMonoCNT)s aligned outside the JAK1-wt/mutant locus

%(nonJAK1cnt)s non-JAK1 aligned reads (showing top-10)
%(top10line)s

'''%{'bamfile':bamfile,'nAlignedReads':nAlignedReads,'monoCNT':monoCNT,'multiCNT':multiCNT,'nJAKs':nJAKs,
'JAK1differentSiteMonoCNT':JAK1differentSiteMonoCNT,'JAK1undefinedMonoCNT':JAK1undefinedMonoCNT,'JAK1mutMonoCNT':JAK1mutMonoCNT,
'JAK1wtMonoCNT':JAK1wtMonoCNT,'nonJAK1cnt':nonJAK1cnt,'top10line':top10line,'nJAKwtMutlocus':nJAKwtMutlocus}
	print logTEXT
	logFile.write(logTEXT)


def analyzeAllRawAlignmentFiles(BAMfiles=None):
	# Determine for each alignment:
	# - cellBC, UMI
	# - aligned gene (separately qualify JAK1-wt / JAK1-mutant / JAK1-undefined), or intergenic
	# Give per cellBC, per UMI, per qualification read counts
	#   . qualification1: mono or multi mapping. 
	#   . qualification2 if mono: JAK1-wt / JAK1-mutant / JAK1-undefined / other gene / intergenic
	perBAMstats,allBAMsStats={},{}
	perCellPerUMIinfo={}
	global perPositionCoverage
	perPositionCoverage=HTSeq.GenomicArray("auto",stranded=True,typecode='i') # a single entry array, only 1 count at each position. 
	for BAMfile in BAMfiles:
		perBAMstats[BAMfile]={"monoCNT":0,"multiCNT":0,"nAlignedReads":0,"JAK1wtMonoCNT":0,"JAK1mutMonoCNT":0,"JAK1undefinedMonoCNT":0,"JAK1differentSiteMonoCNT":0,"nonJAK1MonoCNTperGene":{}}
		BAMrunner=HTSeq.bundle_multiple_alignments(HTSeq.BAM_Reader(BAMfile))
		cnt=0
		while 1:
			alignments=next(BAMrunner,None)
			if not alignments:
				BAMrunner.close()
				break
			cnt+=1
			if cnt%100==0:
				print BAMfile,cnt,'\r',
				sys.stdout.flush()
				#if cnt>10000:break
			cellBC,UMI=giveCellBCandUMIfromAlignment(alignment=alignments[0])
			if not perCellPerUMIinfo.has_key(cellBC):perCellPerUMIinfo[cellBC]={}
			if not perCellPerUMIinfo[cellBC].has_key(UMI):
				perCellPerUMIinfo[cellBC][UMI]={"monoCNT":0,"multiCNT":0,"nAlignedReads":0,"JAK1wtMonoCNT":0,"JAK1differentSiteMonoCNT":0,
				"JAK1mutMonoCNT":0,"JAK1undefinedMonoCNT":0,"nonJAK1MonoCNTperGene":{}}
			perBAMstats[BAMfile]["nAlignedReads"]+=1
			perCellPerUMIinfo[cellBC][UMI]["nAlignedReads"]+=1
			if len(alignments)>1: # multi-mapper
				perBAMstats[BAMfile]["multiCNT"]+=1
				perCellPerUMIinfo[cellBC][UMI]["multiCNT"]+=1
			else:
				perBAMstats[BAMfile]["monoCNT"]+=1
				perCellPerUMIinfo[cellBC][UMI]["monoCNT"]+=1
				alignment=alignments[0]
				Qualification=giveGeneQualificationFromAlignment(alignment=alignment)
				for cigarPiece in alignment.cigar:
					if cigarPiece.type=='M':
						perPositionCoverage[cigarPiece.ref_iv]+=1
				#if Qualification in ["JAK1wtMonoCNT","JAK1mutMonoCNT","JAK1undefinedMonoCNT"]:
				#	return alignment
				if Qualification in ["JAK1wtMonoCNT","JAK1mutMonoCNT","JAK1undefinedMonoCNT","JAK1differentSiteMonoCNT"]:
					perBAMstats[BAMfile][Qualification]+=1
					perCellPerUMIinfo[cellBC][UMI][Qualification]+=1
				else:
					# in this case Qualification becomes either a gene or the designation as "intergenic"
					if not perBAMstats[BAMfile]["nonJAK1MonoCNTperGene"].has_key(Qualification):
						perBAMstats[BAMfile]["nonJAK1MonoCNTperGene"][Qualification]=0
					if not perCellPerUMIinfo[cellBC][UMI]["nonJAK1MonoCNTperGene"].has_key(Qualification):
						perCellPerUMIinfo[cellBC][UMI]["nonJAK1MonoCNTperGene"][Qualification]=0
					perBAMstats[BAMfile]["nonJAK1MonoCNTperGene"][Qualification]+=1
					perCellPerUMIinfo[cellBC][UMI]["nonJAK1MonoCNTperGene"][Qualification]+=1
	return perCellPerUMIinfo,perBAMstats,allBAMsStats


def giveGeneQualificationFromAlignment(alignment=None):
	# analyze the matched alignment (cigar) parts, which genes do they match to?
	cigarPieces=[cigarPiece for cigarPiece in alignment.cigar if cigarPiece.type=='M']
	geneNAMEs=set([feat.attr.get('gene_name') for cigarPiece in cigarPieces for IV,feats in AllFeats[cigarPiece.ref_iv].steps() for feat in feats])
	if 'JAK1' in geneNAMEs:
		# check if the read covered the mutant/wt interval
		mutJAK1cigarPieces=[cigarPiece for cigarPiece in cigarPieces if cigarPiece.ref_iv.contains(WTmutJAKinterval)]
		if len(mutJAK1cigarPieces)>1:raise Exception('fail')
		if len(mutJAK1cigarPieces)==0:
			Qualification="JAK1differentSiteMonoCNT" 
		else:
			mutJAK1cigarPiece=mutJAK1cigarPieces[0]
			readSEQ=alignment.read_as_aligned.seq[mutJAK1cigarPiece.query_from:mutJAK1cigarPiece.query_to] # this (read_as_aligned) makes it into genomic (+ strand) nt's
			genomicPosits=range(mutJAK1cigarPiece.ref_iv.start,mutJAK1cigarPiece.ref_iv.end)
			if all([WTmutPositionNTs['wt'][posit]==readSEQ[genomicPosits.index(posit)] for posit in WTmutPositionNTs['wt'].keys()]):
				Qualification="JAK1wtMonoCNT"
			elif all([WTmutPositionNTs['mut'][posit]==readSEQ[genomicPosits.index(posit)] for posit in WTmutPositionNTs['mut'].keys()]):
				Qualification="JAK1mutMonoCNT"
			else:
				Qualification="JAK1undefinedMonoCNT"
	elif not geneNAMEs:
		Qualification="intergenic"
	elif len(geneNAMEs)>1:
		Qualification='multipleGenenames'
	else:
		Qualification=list(geneNAMEs)[0]
	return Qualification


def giveCellBCandUMIfromAlignment(alignment=None):
	# based on the name, get cell barcode and UMI in nucleotide format
	# from XB:Z:*, get cellBC in indrop-format
	# from XU:Z:*, get UMI in indrop-format
	# check if nt and indrop format UMI are the same
	# check if nucleotide format cell barcode translates into indrop-format cell barcode as expected
	# return cellBC/UMI in indrop-format
	cellBC_nt,UMI_nt=alignment.read.name.split(':')[:2]
	cellBC_indrop=dict(alignment.optional_fields)['XB']
	UMI_indrop=dict(alignment.optional_fields)['XU']
	if not UMI_nt==UMI_indrop:raise Exception(alignment.get_sam_line()+' UMI fail')
	if not nucleotideBasedToIndropBasedCellBCs[cellBC_nt]==cellBC_indrop:raise Exception(alignment.get_sam_line()+' cellBC fail')
	return cellBC_indrop,UMI_indrop


def makeFigANDwriteDownTXidentificationVSreadDepth(perCellPerUMIinfo=None,sampleName=None,minReadCNTperUMI=None,logFile=None,minWTvsMUTreadCNTdifference=None):
	#9. Reads are downsamples, and number of identified transcripts are assessed. This is plotted against each other to determine library saturation. 
	#   - Shown as plot JAK1analysis_XYplot_readCountVStranscriptCount_<sampleName>.pdf
	# - Use downsampling fractions 0.01, 0.02, 0.04, 0.08, 0.16, 0.25, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0
	logFile.write("#### read depth versus transcript identification ####\n")
	Xs,Ys=[],[]
	for dsFraction in 0.01, 0.02, 0.04, 0.08, 0.16, 0.25, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0:
		perCellPerUMIinfo_downSampled={}
		readCNT=0
		#cnt=0
		for cellBC in perCellPerUMIinfo.keys():
			#cnt+=1
			#print dsFraction,cnt,len(perCellPerUMIinfo),'\r',
			#sys.stdout.flush()
			perCellPerUMIinfo_downSampled[cellBC]={}
			for UMI in perCellPerUMIinfo[cellBC].keys():
				perCellPerUMIinfo_downSampled[cellBC][UMI]={'nonJAK1MonoCNTperGene':{}}
				for item in ["JAK1wtMonoCNT","JAK1differentSiteMonoCNT","JAK1mutMonoCNT","JAK1undefinedMonoCNT"]:
					currentCNT=perCellPerUMIinfo[cellBC][UMI][item]
					newCNT=sum(numpy.random.random(currentCNT)<dsFraction)
					readCNT+=newCNT
					perCellPerUMIinfo_downSampled[cellBC][UMI][item]=newCNT
				for gn in perCellPerUMIinfo_downSampled[cellBC][UMI]['nonJAK1MonoCNTperGene']:
					currentCNT=perCellPerUMIinfo[cellBC][UMI]['nonJAK1MonoCNTperGene'][gn]
					newCNT=sum(numpy.random.random(currentCNT)<dsFraction)
					readCNT+=newCNT
					perCellPerUMIinfo_downSampled[cellBC][UMI]['nonJAK1MonoCNTperGene'][gn]=newCNT
		# getting per cell wt/mutant JAK1 transcript counts
		perCellJAK1transcriptCounts={}
		# - If cell/UMI has multiple assignment from different reads, pick the one with most reads for it. 
		# - only pick if at least minReadCNTperUMI read count
		for cellBC in perCellPerUMIinfo.keys():
			perCellJAK1transcriptCounts[cellBC]={'wt':0,'mutant':0}
			for UMI in perCellPerUMIinfo[cellBC].keys():
				# for "JAK1wtMonoCNT","JAK1differentSiteMonoCNT","JAK1mutMonoCNT","JAK1undefinedMonoCNT", other gene, select the one with most reads
				perQualificationReadCounts=[]
				for item in ["JAK1wtMonoCNT","JAK1differentSiteMonoCNT","JAK1mutMonoCNT","JAK1undefinedMonoCNT"]:
					perQualificationReadCounts.append((perCellPerUMIinfo_downSampled[cellBC][UMI][item],item))
				for gn in perCellPerUMIinfo_downSampled[cellBC][UMI]['nonJAK1MonoCNTperGene']:
					perQualificationReadCounts.append((perCellPerUMIinfo_downSampled[cellBC][UMI]['nonJAK1MonoCNTperGene'][gn],gn))
				perQualificationReadCounts.sort(reverse=True)
				WTplusMUTreadTotal=sum([cnt for cnt,qualif in perQualificationReadCounts if qualif in ["JAK1wtMonoCNT","JAK1mutMonoCNT","JAK1undefinedMonoCNT"]])
				if (perQualificationReadCounts[0][1] in ["JAK1wtMonoCNT","JAK1mutMonoCNT","JAK1undefinedMonoCNT"]) and (not perQualificationReadCounts[0][0]>=(minWTvsMUTreadCNTdifference*WTplusMUTreadTotal)):
					continue
				if perQualificationReadCounts[0][0]>perQualificationReadCounts[1][0] and perQualificationReadCounts[0][0]>=minReadCNTperUMI: # one qualification with heighest read count, that is at least minReadCNTperUMI
					umiQualification=perQualificationReadCounts[0][1]
					if umiQualification=="JAK1wtMonoCNT":
						perCellJAK1transcriptCounts[cellBC]['wt']+=1
					if umiQualification=="JAK1mutMonoCNT":
						perCellJAK1transcriptCounts[cellBC]['mutant']+=1
		txCNT=sum([perCellJAK1transcriptCounts[cellBC]['wt']+perCellJAK1transcriptCounts[cellBC]['mutant'] for cellBC in perCellJAK1transcriptCounts.keys()])
		Xs.append(readCNT)
		Ys.append(txCNT)
		logFile.write("%s fraction of reads (n=%s): %s transcripts\n"%(dsFraction,readCNT,txCNT))
		print "%s fraction of reads (n=%s): %s transcripts\n"%(dsFraction,readCNT,txCNT)
	pylab.close()
	pylab.plot(Xs,Ys,'r.',lw=None)
	pylab.xlabel('read depth')
	pylab.ylabel('transcript count')
	pylab.savefig(outputDirectory+'JAK1analysis_XYplot_readCountVStranscriptCount_%s.pdf'%(sampleName))
	print 'scp jcalis@heracles.rockefeller.edu:%sJAK1analysis_XYplot_readCountVStranscriptCount_%s.pdf .'%(outputDirectory,sampleName)



def plotJAK1mutantVersusJAK1wildtypeTranscriptCounts(sampleName=None,outputDirectory=None):
	WTcnts,PATcnts=[],[]
	afile=open(outputDirectory+'JAK1analysis_JAK1transcriptCounts_perCell_%s.txt'%(sampleName),'r')
	SKIPheaderLine=afile.readline()
	perXYcnts={}
	while 1:
		aline=afile.readline()
		if not aline:
			afile.close()
			break
		cellID,wtCnt,patCnt=aline.replace('\n','').split('\t')
		wtCnt,patCnt=string.atoi(wtCnt),string.atoi(patCnt)
		wtCnt,patCnt=max(0.5,wtCnt),max(0.5,patCnt)
		WTcnts.append(wtCnt)
		PATcnts.append(patCnt)
		if not perXYcnts.has_key((wtCnt,patCnt)):
			perXYcnts[(wtCnt,patCnt)]=0
		perXYcnts[(wtCnt,patCnt)]+=1
	print len(PATcnts)
	XYkeys=perXYcnts.keys()
	#XYkeys.remove((0.5,0.5))
	Xs,Ys=[XYkey[0] for XYkey in XYkeys],[XYkey[1] for XYkey in XYkeys]
	cellCNTs=[10*perXYcnts[XYkey] for XYkey in XYkeys]
	pylab.close()
	pylab.scatter(Xs,Ys,c='k',s=cellCNTs,alpha=0.3,linewidths=0)
	pylab.xlabel('JAK1 wild-type count')
	pylab.ylabel('JAK1 patient count')
	pylab.xlim((0.3,100))
	pylab.ylim((0.3,100))
	pylab.loglog()
	pylab.savefig(outputDirectory+'JAK1analysis_XYscatterPlot_perCellJAK1wtVSpatientTranscriptCount_%s.pdf'%(sampleName))


if __name__ == '__main__':
	# loading gtf data 
	print time.ctime(),'Loading gtf file'
	SetUpGlobalObjectsForHumanOnlyMapping()
	print time.ctime(),'Done loading gtf file'
	### Run each sample
	LogFiles={}
	for sampleName in sampleNames:
		outputStatsFile=outputDirectory+'JAK1analysis_stats_%s.txt'%(sampleName)
		LogFiles[sampleName]=open(outputStatsFile,'w')
		LogFiles[sampleName].write('''############################ Stats for run on %(sampleName)s ############################
Run start: %(time)s
Output directory: %(outputDirectory)s
Sample name: %(sampleName)s 
BAM(s): 
 - %(BAMline)s

	'''%{'sampleName':sampleName,'time':time.ctime(),'outputDirectory':outputDirectory,'BAMline':string.join([BAMfile for BAMfile in perSampleBams[sampleName]],'\n - ')})
		LogFiles[sampleName].write('Running in EnsemblMode %s:\n - GTFfileName: %s\n'%(EnsemblMode,GTFfileName))
		perCellPerUMIinfo,perBAMstats,allBAMsStats=analyzeAllRawAlignmentFiles(BAMfiles=perSampleBams[sampleName])
		writeDownPerBAMstats(perBAMstats=perBAMstats,logFile=LogFiles[sampleName])
		minReadCNTperUMI=2
		minWTvsMUTreadCNTdifference=0.9 # interpret this as at least 0.9 (ie 90%) of the wt+patient reads must be either from the patient or from the mutant
		writeDownAllBAMsStats(perCellPerUMIinfo=perCellPerUMIinfo,logFile=LogFiles[sampleName],minReadCNTperUMI=minReadCNTperUMI,minWTvsMUTreadCNTdifference=minWTvsMUTreadCNTdifference)
		makeAllBAMbasedHistograms_andMore(logFile=LogFiles[sampleName],perCellPerUMIinfo=perCellPerUMIinfo,sampleName=sampleName)
		perCellJAK1transcriptCounts=writeDownANDgivePerCellJAK1wildtypeAndMutantTranscriptCounts(perCellPerUMIinfo=perCellPerUMIinfo,sampleName=sampleName,minReadCNTperUMI=minReadCNTperUMI,minWTvsMUTreadCNTdifference=minWTvsMUTreadCNTdifference)
		writeDownTXidentificationStats(perCellJAK1transcriptCounts=perCellJAK1transcriptCounts,logFile=LogFiles[sampleName],minReadCNTperUMI=minReadCNTperUMI,minWTvsMUTreadCNTdifference=minWTvsMUTreadCNTdifference)
		makeFigANDwriteDownTXidentificationVSreadDepth(perCellPerUMIinfo=perCellPerUMIinfo,sampleName=sampleName,minReadCNTperUMI=minReadCNTperUMI,logFile=LogFiles[sampleName],minWTvsMUTreadCNTdifference=minWTvsMUTreadCNTdifference)
		LogFiles[sampleName].close()
		plotJAK1mutantVersusJAK1wildtypeTranscriptCounts(sampleName=sampleName,outputDirectory=outputDirectory)
		print time.ctime(),'finished running sample ',sampleName
	# done
	print time.ctime(),'finished running all samples'









