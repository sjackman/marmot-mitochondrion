# Yellow-bellied marmot (*Marmota flaviventris*)

# Draft genome assembly
draft=marmot-mt

# Number of threads
t=16

# Parallel gzip
gzip=pigz -p$t

# Paths to programs.
abyss_bin=/gsc/btl/linuxbrew/Cellar/abyss/2.1.1-k128/bin
pilon_jar=/gsc/btl/linuxbrew/Cellar/pilon/1.22/pilon-1.22.jar

# Report run time and memory usage
export SHELL=zsh -opipefail
export REPORTTIME=1
export TIMEFMT=time user=%U system=%S elapsed=%E cpu=%P memory=%M job=%J
time=command time -v -o $@.time

.DELETE_ON_ERROR:
.SECONDARY:

all: \
	marmot1-20M.bx.unicycler.fa \
	marmot-mt.mitos.fix.gbf.pdf \
	marmot-mt.bwa.NC_018367.1.cds.sort.bam.bai

# NCBI

# Download CDS amino acid sequence from NCBI using Entrez Direct.
NC_%.cds.faa:
	efetch -db nuccore -id NC_$* -format fasta_cds_aa | seqtk seq >$@

# Download CDS DNA sequence from NCBI using Entrez Direct.
NC_%.cds.fa:
	efetch -db nuccore -id NC_$* -format fasta_cds_na | seqtk seq >$@

# Download a nucleotide sequence from NCBI using Entrez Direct.
NC_%.fa:
	efetch -db nuccore -id NC_$* -format fasta | seqtk seq >$@

# Download a nucleotide sequence from NCBI using ncbi-acc-download.
NC_%_0.fa:
	ncbi-acc-download -Ffasta -o NC_$* NC_$*

# Data

# Take the first ten million read pairs.
marmot1-20M.fq.gz: \
		data/lane4_ACAGAGGT--TATAGTTG--CGGTCCCA--GTCCTAAC_S1_L004_R1_001.fastq.gz \
		data/lane4_ACAGAGGT--TATAGTTG--CGGTCCCA--GTCCTAAC_S1_L004_R2_001.fastq.gz
	seqtk mergepe $^ | head -n80000000 | $(gzip) >$@

# samtools

# Index a FASTA file.
%.fa.fai: %.fa
	samtools faidx $<

# Sort a SAM file and produce a sorted BAM file.
%.sort.bam: %.sam.gz
	samtools sort -@$t -T$$(mktemp -u -t $@.XXXXXX) -o $@ $<

# Index a BAM file.
%.bam.bai: %.bam
	samtools index -@$t $<

# BWA

# Index the target genome.
%.fa.bwt: %.fa
	bwa index $<

# Align a FASTA file to the draft assembly and sort by position.
$(draft).bwa.%.sort.bam: %.fa $(draft).fa.bwt
	bwa mem $(draft).fa $< | samtools sort -@$t -T$$(mktemp -u -t $@.XXXXXX) -o $@

# Align linked reads to the draft genome and sort by position.
%.$(lr).bx.sort.bam: %.fa.bwt $(lr).bx.fq.gz
	bwa mem -t$t -pC $*.fa $(lr).bx.fq.gz | samtools sort -@$t -T$$(mktemp -u -t $@.XXXXXX) -o $@

# EMA

# Download the barcode white list.
4M-with-alts-february-2016.txt:
	curl -o $@ https://raw.githubusercontent.com/10XGenomics/supernova/master/tenkit/lib/python/tenkit/barcodes/4M-with-alts-february-2016.txt

# Count barcodes.
%.ema-ncnt: %.fq.gz 4M-with-alts-february-2016.txt
	ema count -w 4M-with-alts-february-2016.txt -o $* $<

# Extract the barcode to BX:Z tag using ema preproc.
%.bx.fq.gz: %.fq.gz %.ema-ncnt
	gunzip -c $< | ema preproc -t$t -b -n1 -w 4M-with-alts-february-2016.txt -o $*.ema $*.ema-ncnt
	$(gzip) <$*.ema/ema-bin-000 >$@
	rm -rf $*.ema

# Align linked reads to the draft genome using EMA and sort by position.
%.$(lr).bx.ema.sort.bam: $(lr).bx.fq.gz %.fa.bwt
	$(time) ema align -t$t -r $*.fa -1 $< \
	| samtools view -@$t -u -F4 \
	| samtools sort -@$t -T$$(mktemp -u -t $@.XXXXXX) -o $@

# ntHash

# Count k-mers using ntCard.
%.ntcard_k32.hist: %.fq.gz
	ntcard -t$t -c1000 -k 32,64,96,128 -p $*.ntcard $<

# Convert a .hist to a .histo file for GenomeScope.
%.histo: %.hist
	sed -n 's/^f//p' $< | tr '\t' ' ' >$@

# ABySS
k=128

# Assemble short reads using ABySS.
%.k$k.abyss.fa: %.fq.gz
	mkdir $*.k$k.abyss
	$(abyss_bin)/abyss-pe -C $*.k$k.abyss mpirun=mpirun np=$t k=$k v=-v e=5 E=1 c=5 name=marmot lib=pe1 pe1=../$<
	ln -s $*.k$k.abyss/marmot-scaffolds.fa $@

# Unicycler

# Select the first read of the read pair.
%.1.fq.gz: %.fq.gz
	seqtk seq -1 $< | $(gzip) >$@

# Select the second read of the read pair.
%.2.fq.gz: %.fq.gz
	seqtk seq -2 $< | $(gzip) >$@

# Assemble short reads using Unicycler.
%.unicycler.fa: %.1.fq.gz %.2.fq.gz
	unicycler -t$t -o $*.unicycler -1 $*.1.fq.gz -2 $*.2.fq.gz
	seqtk seq $*.unicycler/assembly.fasta >$@

# Polish the assembly using Illumina reads.
# Correct only single nucleotide errors and small indels.
# Do not correct local misassemblies.
%.unicycler-polish.stamp: %.fa \
		marmot1-20M.bx.1.fq.gz \
		marmot1-20M.bx.2.fq.gz
	rm -rf $*.unicycler-polish
	mkdir $*.unicycler-polish
	cd $*.unicycler-polish \
	&& $(time) unicycler_polish --threads $t --no_fix_local --pilon $(pilon_jar) \
		-a ../$< \
		-1 ../marmot1-20M.bx.1.fq.gz \
		-2 ../marmot1-20M.bx.2.fq.gz
	touch $@

# Copy the final FASTA file from unicycler_polish.
%.unicycler-polish.fa: %.unicycler-polish.stamp
	seqtk seq $*.unicycler-polish/???_final_polish.fasta | sed 's/ LN:i:[0-9]*//' >$@

# MAFFT

# Multiple sequence alignment using MAFFT.
%.mafft.fa: %.fa NC_018367.1.fa NC_025277.1.fa NC_026705.1.fa NC_026706.1.fa NC_027278.1.fa NC_027283.1.fa NC_031209.1.fa NC_031210.1.fa NC_032370.1.fa NC_032371.1.fa NC_032372.1.fa NC_032373.1.fa NC_032374.1.fa NC_032375.1.fa NC_032376.1.fa NC_035576.1.fa
	cat $^ >$*.mafft.in.fa
	mafft --thread $t --auto $*.mafft.in.fa >$*.mafft.out.fa
	seqtk seq $*.mafft.out.fa >$*.mafft.fa
	rm -f $*.mafft.in.fa $*.mafft.out.fa

# FastTree

# Construct a phylogentic tree using FastTree.
%.fasttree.tree: %.fa
	FastTree -nt -quote $< >$@

# IQ-TREE

# Construct a phylogentic tree using IQ-TREE.
%.iqtree.tree: %.fa
	sed 's/ /_/g' $< >$*.iqtree.fa
	iqtree -nt $t -seed 0 -bb 1000 -s $*.iqtree.fa -pre $*.iqtree
	mv $*.iqtree.treefile $*.iqtree.tree
	mv $*.iqtree.contree $*.iqtree.con.tree
	rm $*.iqtree.fa

# tbl2asn

# Add structured comments to the FASTA file.
%.fsa: %.fa
	sed '/^>/s/ .*/ [organism=Marmota flaviventris] [location=mitochondrion] [completeness=complete] [topology=circular] [mgcode=2]/' $< >$@

# Convert annotations from Genbank TBL to GBF.
%.gbf: %.tbl marmot-mt.fsa marmot-mt.cmt marmot-mt.sbt
	ln -sf marmot-mt.fsa $*.fsa
	tbl2asn -a s -i $*.fsa -w marmot-mt.cmt -t marmot-mt.sbt -Z $*.discrep -Vbv
	mv errorsummary.val $*.errorsummary.val

# OrganellarGenomeDRAW

# Render the genome annotation to PNG.
%.gbf.png: %.gbf
	drawgenemap --format png --infile $< --outfile $< --gc --density 126
	mogrify -units PixelsPerInch -density 300 $@

# Render the genome annotation to Postscript.
%.gbf.ps: %.gbf
	drawgenemap --format ps --infile $< --outfile $< --gc --density 126

# Convert Postscript to PDF.
%.pdf: %.ps
	ps2pdf $< $@