# methylK 
_methylK_ is a directory of scripts used within _(Chatterton et al, Methods for detecting Brain-Cell derived Cell-Free DNA, 2020)_

Background: DNA methylation is an epigenetic modification that is intricately involved in cell-specification and cell-function. Bisulfite conversion of DNA followed by Next Generation Sequencing (NGS) allows the analysis of DNA methylation at the single molecule level. The cell-specificity of DNA methylation patterns can be used to identify the cell-of-origin of DNA molecules, such as deconvolution of cell-free DNA (cfDNA). To acheive this, we binarize DNA methylation patterns of primary cells-types within genomic context, creating cell-specific DNA methylation reference genomes for the assignment of DNA fragments of unknown source (i.e. cfDNA) to their cell-of-origin.

## Instructions

### File Formats 

 targets: a tab-delimited table (.txt) with the headers "sample_seqname" "tissue" "type". Each sample to be analysed is represented by one row. The sample_seqname column refers to the sample names of the .fastq files i.e. "sample_seqname".R{1/2}.fastq.gz. The "tissue" column is used to define which tissue group this sample belongs and can be any value eg. "PBMC". The "type" column defines how that sample will be used within the analysis and can be one of 3 values "interest", "contrast" or "identify". For type "interest" and "contrast" a methylotype.fasta file will be made for each "tissue" x "type". 

 genome: reference genome file (.fa). Needs to sorted. This can be done using seqkit (seqkit sort -i in.fa -o out.fa)
 fastq: files need to be named "sample_seqname".R{1/2}.fastq.gz. 

 note - If you are starting from ".bed" or ".bedGraph" files need to be in the format chr,start,end,methylation_percentage

### Software 

[bedtools](https://bedtools.readthedocs.io/en/latest/),
[bismark](https://www.bioinformatics.babraham.ac.uk/projects/bismark/)  (0.18.2),
[bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) (2.2.5),
[gatk](https://gatk.broadinstitute.org/hc/en-us) (3.6.0),
[kallisto](https://github.com/pachterlab/kallisto_paper_analysis) (0.46.0),
[parallel](https://www.gnu.org/software/parallel/) (20160222),
[picard](https://broadinstitute.github.io/picard/) (2.7.1),
[python](https://www.python.org/) (2.7.9),
[R](https://www.r-project.org/) (3.6.3),
[samtools](http://www.htslib.org/) (1.9),
[trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) (0.36)

note - We have run the analysis on a Linux operating system. All of the software should all be accesible from $PATH

	export PATH=ngsutils/bin:$PATH # ngsutils git
	GATK_JAR=gatk
	PICARD_JAR=picard
	njobs=4 # Also set to make server happy when using parallel

## Quickrun

	# Inputs #
	methylK_dir=/methylK
	genome=$methylK_dir/tNGBS_n33_lambda1.3.fa

	# test directory contains .fastq files and targets.txt
	sdir=$methylK_dir/test # sample directory containing all PE .fastq files 
	targets=$methylK_dir/test/targets.txt # location of targets file
	odir=$methylK_dir/test/output # ouput directory is where all results will be written
	genome_bismark=$odir/Bisulfite_Genome
	bed_location=$odir/bismark_results

	# hard-coded thresholds/ values
	cov_threshold=5 # coverage threshold in order to include cytosine in analysis
	bin_threshold=50 # methylation % threshold for binarizing DNA methylation i.e. < 50% = T, > 50% = C
	read_length=25 # length that reads will be truncated. This should be the shortest read length

	# Run #
	cd $odir
	$methylK_dir/fastq_to_bed.sh $methylK_dir $genome $sdir $odir $targets
	$methylK_dir/bismark_to_fasta.sh $odir $targets $cov_threshold $genome $genome_bismark $bed_location $bin_threshold $GATK_JAR $PICARD_JAR
	$methylK_dir/quant_mk_cell.sh $methylK_dir $targets $sdir $odir $read_length
	Rscript --vanilla $methylK_dir/snr_calc.R $odir $methylK_dir $targets
	$methylK_dir/quant_mk_cfdna.sh $methylK_dir $targets $sdir $odir $read_length
	Rscript --vanilla $methylK_dir/snr_quant.R $odir $methylK_dir $targets

	note - we provide the .kidx file (cell_methylotype.kidx) that can be used to quantify Neuron and Glia-cfDNA if the assays used within (Chatterton et al, Methods for detecting Brain-Cell derived Cell-Free DNA, 2020) have been applied to cfDNA. To do so uncomment line 112 in quant_mk_cfdna.sh.


## Step-by-step

## 1. FASTQ to meth
Overview - Alignment and DNA methylation extraction needs to be performed on all sample types "interest" or "contrast". If this step has already been completed and DNA methylation has been extracted and stored in either .tsv, .bedGraph or .bismark.cov format then we can move to "meth to FASTA". For convenience we provide scripts used within _(Chatterton et al, Methods for detecting Brain-Cell derived Cell-Free DNA, 2020)_ for the preparation of bisulfite genome (genome_prepare.sh) alignment and methylation calling using [Bismark](https://www.bioinformatics.babraham.ac.uk/projects/bismark/) (meth_trim_align_call.sh). Note: 1. NGS read quality should be ascertained prior to running using programs such as fastqc. 2. The .fastq files need to be named using the following convention "sample_seqname".R{1/2}.fastq.gz 3. The scripts perform alignment and calling of all samples of "types" (column 3 of targets.txt) "interest" and "contrast". 4. The NGS libraries were created with Illumina Nextera and are trimmed using these adapter sequences. Using the test data this should run in ~20mins using 4 CPU and 16Gb RAM.

	# Inputs #
	methylK_dir=methylK/
	genome=$methylK_dir/tNGBS_n33_lambda1.3.fa

	# test directory contains .fastq files and targets.txt
	sdir=$methylK_dir/test # sample directory containing all PE .fastq files 
	targets=$methylK_dir/test/targets.txt # location of targets file
	odir=$methylK_dir/test/output # ouput directory is where all results will be written

	# Run #
	cd $odir # run from output directory
	$methylK_dir/fastq_to_bed.sh $methylK_dir $genome $sdir $odir $targets


## 2. meth to FASTA
Overview - Here we binerize the DNA methylation values depending on our "bin_threshold". For instance bin_threshold=50 results in any cytosine with DNA methylation >=50% being translated into a "C" with FASTA format, converesly cytosines with DNA methylation <50% are tranlsated into "T".The script will produce a FASTA file (methylotype.fasta) for each "tissue" (column 2 of targets.txt) of types (column 3 of targets.txt) "interest" and "contrast". In order to create a k-mer index of sample types "interest" and "contrast" we need to combine the methylotype.fasta files from each "tissue" into master_methylotype.fasta file. Base changes that are dependent on the DNA methylation of each tissue are introduced into the .fasta file. The base changes enable the DNA methylation to be contextualised within the genome sequence so that DNA methylation based k-mers can be characterised between "tissues". 

Note. We use [Kallisto](https://github.com/pachterlab/kallisto_paper_analysis) software to index k-mers from the master_methylotype.fasta file. It is therefore critical to only use cytosines with coverage accross all "tissues". Cytosines cannot be masked as [Kallisto](https://github.com/pachterlab/kallisto_paper_analysis) looks for ACGUT comptibility and if none is found, such as in the case of "N", the nucleotide is replaced with with a random nucleotide. Therefore we break-up the .fasta sequences at cytosine that do not have coverage accross all sample types "interest" and "contrast".

	# Inputs #
	genome_bismark=$odir/Bisulfite_Genome
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
Overview - The number of methylated cytosines within each assigned read are counted and the proportion of assigned reads with co-methylation events above the signal-to-noise threshold are defined "frac_upsdfq_thresh" for all tissue:interest (mkdf.txt).

	# Run #
	Rscript --vanilla $methylK_dir/snr_quant.R $odir $methylK_dir $targets






