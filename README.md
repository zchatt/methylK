# methylK 
_methylK_ is a directory of scripts used within [(Chatterton et al, Front. Mol. Neurosci., 2021)](https://www.frontiersin.org/articles/10.3389/fnmol.2021.672614/full)

Background: DNA methylation is an epigenetic modification critical for cell-specification and cell-function that is covlently bound to cytosine.Bisulfite conversion of DNA followed by Next Generation Sequencing (NGS) produces single-base level DNA methylation information on NGS reads. Cell-specific DNA methylation patterns can be used to identify the cell-of-origin of DNA molecules, such as the origin of cell-free DNA (cfDNA) molecules. To acheive this, we binarize DNA methylation of primary cells-types within genomic context, creating cell-specific DNA methylation reference for the assignment of NGS reads mixed source (i.e. cfDNA) to their cell-of-origin.

## Instructions

### File Formats 

targets: a tab-delimited .txt file with headers "sample_seqname" "tissue" "type". One sample per row. The "sample_seqname" refers to the sample names eg. "sample_seqname".R{1/2}.fastq.gz. The "tissue" defines the tissue group eg. "PBMC". The "type" column defines the samples use in the analysis and are either "interest" or "contrast" (tissues to be deconvoluted) or "identify" eg. cfDNA.

genome: sorted FASTA reference file. This can be done using seqkit (seqkit sort -i in.fa -o out.fa)

fastq: files need to be named "sample_seqname".R{1/2}.fastq.gz. 

 note - If you are starting from ".bed" or ".bedGraph" files need to be in the format chr,start,end,methylation_percentage

### Software Required
[bedtools](https://bedtools.readthedocs.io/en/latest/),
[bismark](https://www.bioinformatics.babraham.ac.uk/projects/bismark/)  (0.18.2),
[bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) (2.2.5),
[gatk](https://gatk.broadinstitute.org/hc/en-us) (3.6.0),
[kallisto](https://github.com/pachterlab/kallisto_paper_analysis) (0.46.0),
[parallel](https://www.gnu.org/software/parallel/) (20160222),
[picard](https://broadinstitute.github.io/picard/) (2.7.1),
[python](https://www.python.org/) (2.7.9),
[samtools](http://www.htslib.org/) (1.9),
[trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) (0.36) &
[R](https://www.r-project.org/) (3.6.3). In addition the following [R](https://www.r-project.org/) packages should be installed;

	install.packages(c("plyr","dplyr","tidyr","data.table","ggplot2","gridExtra","scales","caret"))

note - All of the software should all be accesible from $PATH. We have run the analysis on a Linux OS (CentOS release 6.9, 2.6.32-696.16.1.el6.x86_64).

	GATK_JAR=gatk
	PICARD_JAR=picard
	njobs=4 # Also set to make server happy when using parallel

## Download
	git clone https://github.com/zchatt/methylK.git
	chmod +x methylK/*

## Quickrun
	# Inputs #
	methylK_dir=$(readlink -f methylK)
	genome=$methylK_dir/tNGBS_n33_lambda1.3.fa
	# note - the methylK/test directory contains .fastq files and targets.txt that can be used to test scripts
	sdir=$methylK_dir/test # sample directory containing all PE .fastq files 
	targets=$methylK_dir/test/targets.txt # location of targets file
	odir=$methylK_dir/test/output # ouput directory is where all results will be written
	mkdir $odir
	genome_bismark=$odir/bismark_genome/Bisulfite_Genome
	bed_location=$odir/bismark_results

	# hard-coded values
	cov_threshold=5 # coverage threshold in order to include cytosine in analysis
	bin_threshold=50 # methylation % threshold for binarizing DNA methylation i.e. < 50% = T, > 50% = C
	read_length=25 # length that reads will be truncated. This should be the shortest read length

	# 0. Set output directory
	cd $odir

	# 1. FASTQ to meth
	$methylK_dir/fastq_to_bed.sh $methylK_dir $genome $sdir $odir $targets
	
	# 2. meth to FASTA
	$methylK_dir/bismark_to_fasta.sh $odir $targets $cov_threshold $genome $genome_bismark $bed_location $bin_threshold $GATK_JAR $PICARD_JAR
	
	# 3. Assignment of tissue:interest & tissue:contrast
	$methylK_dir/quant_mk_cell.sh $methylK_dir $targets $sdir $odir $read_length
	
	# 4. Signal-to-noise calculation
	Rscript --vanilla $methylK_dir/snr_calc.R $odir $methylK_dir $targets
	
	# 5. Assignment of tissue:identify
	$methylK_dir/quant_mk_cfdna.sh $methylK_dir $targets $sdir $odir $read_length
	
	# 6. Quantify tissue:interest from tissue:identify over signal-to-noise
	Rscript --vanilla $methylK_dir/snr_quant.R $odir $methylK_dir $targets

note 1. If the test data has run correctly then the fraction of reads assigned to our tissue:interest (eg."NeuNpos") within our "unknown" samples = 0.00564325
	
	grep -v "NA" $odir/NeuNpos_mkdf.txt | awk '{sum+=$8;} END{print sum;}'

note 2. We provide the .kidx file (cell_methylotype.kidx) that can be used to quantify Neuron and Glia-cfDNA following bisulfite amplicon sequencing using the assays described in [(Chatterton et al,Front. Mol. Neurosci., 2021)](https://www.frontiersin.org/articles/10.3389/fnmol.2021.672614/full). To do so uncomment/ comment line 107/109 of quant_mk_cfdna.sh

## Step-by-step

## 1. FASTQ to meth
Overview - Alignment and DNA methylation extraction is performed on tissue:interest & tissue:contrast. If DNA methylation is already in either .tsv, .bedGraph or .bismark.cov format then begin at [meth to FASTA](#2-meth-to-fasta). For convenience we provide scripts used within [(Chatterton et al,Front. Mol. Neurosci., 2021)](https://www.frontiersin.org/articles/10.3389/fnmol.2021.672614/full) for the preparation of bisulfite genome [genome_prepare.sh](https://github.com/zchatt/methylK/blob/master/genome_prepare.sh) alignment and methylation calling using [Bismark](https://www.bioinformatics.babraham.ac.uk/projects/bismark/) [meth_trim_align_call.sh](https://github.com/zchatt/methylK/blob/master/meth_trim_align_call.sh). 

note 1. NGS read quality should be ascertained prior to running using programs such as fastqc. 
note 2. The .fastq files need to be named using the following convention "sample_seqname".R{1/2}.fastq.gz 
note 3. The scripts perform alignment and calling of all samples tissue:interest & tissue:contrast, unless 6th argument is set to identify which runs tissue:identify.
note 4. The test data should run in ~20 mins using 4 CPU and 16Gb RAM.

	# Inputs #
	methylK_dir=methylK/
	genome=$methylK_dir/tNGBS_n33_lambda1.3.fa
	sdir=$methylK_dir/test # sample directory containing all PE .fastq files 
	targets=$methylK_dir/test/targets.txt # location of targets file
	odir=$methylK_dir/test/output # ouput directory is where all results will be written

	# Run #
	cd $odir # run from output directory
	$methylK_dir/fastq_to_bed.sh $methylK_dir $genome $sdir $odir $targets

	# alignment of tissue:identify can be performed by setting the sixth argument;
	# $methylK_dir/fastq_to_bed.sh $methylK_dir $genome $sdir $odir $targets identify 

## 2. meth to FASTA
Overview - DNA methylation fractional measurements are binarized eg. bin_threshold=50 results in any cytosine with DNA methylation < 50% being translated into a "T" within FASTA format. The script will produce a FASTA file (methylotype.fasta) for each tissue:interest & tissue:contrast that are combined and k-mers are indexed using [Kallisto](https://github.com/pachterlab/kallisto_paper_analysis) software.

note 1. Binarisation and insertion of DNA methylation contextualises it within the genome sequence that enable tissue/ cell-specific DNA methylation k-mers. 
note 2. It is critical to only use cytosines with coverage accross all tissue:interest & tissue:contrast. Cytosines cannot be masked as [Kallisto](https://github.com/pachterlab/kallisto_paper_analysis) looks for ACGUT comptibility and if none is found, such as in the case of "N", the nucleotide is replaced with with a random nucleotide. To avoid this we break the FASTA sequences at cytosine without coverage accross all tissue:interest & tissue:contrast.

	# Inputs #
	genome_bismark=$odir/bismark_genome/Bisulfite_Genome
	bed_location=$odir/bismark_results
	bin_threshold=50 # methylation % threshold for binarizing DNA methylation i.e. < 50% = T, > 50% = C
	cov_threshold=5 # threshold for coverage to include cytosine in analysis

	# Run #
	$methylK_dir/bismark_to_fasta.sh $odir $targets $cov_threshold $genome $genome_bismark $bed_location $bin_threshold $GATK_JAR $PICARD_JAR

## 3. Assignment of tissue:interest & tissue:contrast
Overview - Here we index k-mers within the master_methylotype.fasta file using the [Kallisto](https://github.com/pachterlab/kallisto_paper_analysis) "index" function. Trimmed, paired and truncated reads (.fastq) from the tissue type "interest" and "contrast" are then assigned (pseudoaligned) using the [Kallisto](https://github.com/pachterlab/kallisto_paper_analysis) "quant" function. The shortest read length sequenced within the population needs to be supplied. This serves two functions; firstly, we use this value to truncate all reads to the same length prior to running "quant". Secondly, Kallisto's maximum k-mer length is 31. As some groups may have sequenced using 2 x 26bp we use this value to define the k-mer size eg. 26bp=25-mer & >31bp=31-mer. 

	# Inputs #
	read_length=25

	# Run #
	$methylK_dir/quant_mk_cell.sh $methylK_dir $targets $sdir $odir $read_length

## 4. Signal-to-noise calculation
Overview -  The number of cytosine's and guanines are counted within each read assigned to a CT/GA reference sequences of tissue type "interest" and "contrast".Signal-to-noise thresholds are calculated for each co-methylation event (eg. 1:26[read length]) for all reads assigned to tissue type "interest" and "contrast". The optimal signal-to-noise threshold for quantifying reads derived from tissue type "interest" are determined by the point of inflection of the F1-statistic vs signal-to-noise. Two output files will be produced for each tissue type "contrast" 1) a signal-to-noise ratio file that contains information on the assay, co-methylation, true & false-positive rates used for F1 and the signal-to-noise ratio eg. NeuNpos_snr.txt, and 2) signal-to-noise ratio file that is subset with the signal-to-noise threshold, eg. NeuNpos_thresh-2950-snr.txt, that is used in downstream quantification of cfDNA.

	# Run #
	Rscript --vanilla $methylK_dir/snr_calc.R $odir $methylK_dir $targets

## 5. Assignment of tissue:identify
Overview - The paried-end .fastq files from the tisue type "identify" are truncated to $read_length and are then assigned to tissue type "interest" and "contrast" using the k-mer index (above). The number of methylated cytosines within each assigned read are then counted and the proportion of assigned reads with co-methylation events above the signal-to-noise threshold are defined as tissue:interest.

	# Run #
	$methylK_dir/quant_mk_cfdna.sh $methylK_dir $targets $sdir $odir $read_length

## 6. Quantify tissue:interest from tissue:identify over signal-to-noise
Overview - The number of methylated cytosines within each assigned read are counted and the proportion of assigned reads with co-methylation events above the signal-to-noise threshold are defined within the "frac_upsdfq_thresh" for all tissue:interest results files "\_mkdf.txt". 

	# Run #
	Rscript --vanilla $methylK_dir/snr_quant.R $odir $methylK_dir $targets



