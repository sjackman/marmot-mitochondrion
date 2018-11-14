# Yellow-bellied marmot (*Marmota flaviventris*)

# Draft genome assembly
draft=marmot-mt

# Linked reads
lr=marmot1-20M

# Nanopore reads
nanopore=SJ83_1

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
	SRR7878800-20M.bx.k59.unicycler.mt.rot.fa \
	marmot-mt.mitos.fix.gbf.pdf \
	marmot-mt.bwa.NC_018367.1.cds.sort.bam.bai

nanopore: SJ83_1.marmot-xx.paf.gz.pdf

# NCBI

# Download the Woodchuck hepatitis virus DNA sequence from NCBI.
M60765.1.fa:
	curl 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?retmode=text&id=M60765.1&db=nucleotide&rettype=fasta' | seqtk seq >$@

# Download a nucleotide sequence from NCBI.
NC_%.fa:
	curl 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?retmode=text&id=NC_$*&db=nucleotide&rettype=fasta' | seqtk seq >$@

# Download CDS amino acid sequence from NCBI.
NC_%.cds.faa:
	curl 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?retmode=text&id=NC_$*&db=nucleotide&rettype=fasta_cds_aa' | seqtk seq >$@

# Download CDS DNA sequence from NCBI.
NC_%.cds.fa:
	curl 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?retmode=text&id=NC_$*&db=nucleotide&rettype=fasta_cds_na' | seqtk seq >$@

# Download a GBF (GenBank flat file) from NCBI.
NC_%.gbf:
	curl 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?retmode=text&id=NC_$*&db=nucleotide&rettype=gb' >$@

# Download a TBL (feature table) from NCBI.
NC_%.tbl:
	curl 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?retmode=text&id=NC_$*&db=nucleotide&rettype=ft' >$@

# Download a GFF file from NCBI.
NC_%.gff:
	curl 'https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=gff3&id=NC_$*' >$@

# sratoolkit

# Download the first ten million read pairs of the 10x Chromium sequencing data from NCBI SRA.
SRR7878800-20M.fq.gz:
	fastq-dump -Z --split-spot SRR7878800 | sed 's/ .*//;s/^+.*/+/' | head -n 80000000 | $(gzip) >$@

# Download the 10x Chromium sequencing data from NCBI SRA.
SRR7878800.fq.gz:
	fastq-dump -Z --split-spot SRR7878800 | sed 's/ .*//;s/^+.*/+/' | $(gzip) >$@

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

# Compute depth of coverage.
%.bam.depth: %.bam
	samtools depth -a $< >$@

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
k=59

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

# Assemble short reads using Unicycler with a specified value of k.
%.k$k.unicycler.fa: %.1.fq.gz %.2.fq.gz
	unicycler -t$t --kmers=$k --no_correct -o $*.k$k.unicycler -1 $*.1.fq.gz -2 $*.2.fq.gz
	seqtk seq $*.k$k.unicycler/assembly.fasta >$@
	ln -s $*.k$k.unicycler/assembly.gfa $*.k$k.unicycler.gfa

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

# Extract the mitochondrial genome from the assembly graph.
%.mt.gfa: %.gfa NC_018367.1.fa
	Bandage reduce $< $@ --scope aroundblast --query NC_018367.1.fa

# Convert GFA to FASTA.
%.mt.fa: %.mt.gfa
	awk '/^S/ { print ">marmot-mt Marmota flaviventris mitochondrion, complete genome\n" $$3 }' $< >$@

# Rotate the mitochondrial genome to trnF-GAA.
%.mt.rot.fa: %.mt.fa
	sed 's/GTTAATGTAGCTTAATCT/ &/' $< | awk '/^>/ { print; next } {print $$2 $$1}' >$@

# Rotate the unknown circular molecule to the first ORF following trnV-TAC.
%.xx.rot.fa: %.xx.fa
	sed 's/ATGGTTTACTGTTTTGGA/ &/' $< | awk '/^>/ { print; next } {print $$2 $$1}' >$@

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

# Prokka

# Convert the FASTA file to the Prokka FASTA format
%.cds.prokka.faa: %.cds.faa
	gsed -E \
		-e 's/^>([^ ]*) .*gene=([^]]*).*protein=([^]]*).*$$/>\1 ~~~\2~~~\3/' \
		-e 's/^-//' \
		-e 's/ATP/atp/' -e 's/COX/cox/' -e 's/CYTB/cob/' -e 's/ND/nad/' \
		$< >$@

# Annotate genes using Prokka.
prokka/%.gff: %.fa NC_018367.1.cds.prokka.faa
	prokka --cpus $t --kingdom Mitochondria --gcode 2 --addgenes \
		--genus Marmota --species 'flaviventris mitochondrion' \
		--centre BCGSC \
		--rfam --proteins NC_018367.1.cds.prokka.faa \
		--force --outdir prokka --prefix $* \
		$<

# Annotate bacterial genes using Prokka with genetic code 4.
%.gcode4.prokka/marmot.gff: %.fa
	prokka --cpus $t --kingdom Bacteria --gcode 4 --addgenes \
		--genus Unknown --species 'circular DNA molecule' \
		--centre BCGSC \
		--force --outdir $*.gcode4.prokka --prefix marmot \
		$<

# Fix up the Prokka GFF file.
%.prokka.gff: prokka/%.gff
	gsed -e '/^##FASTA/,$$d' -e 's/gnl[^\t]*_1/marmot-mt/' $< >$@

# ARAGORN

# Annotate tRNA using ARAGORN and output a plain text report.
%.aragorn.txt: %.fa
	aragorn -gcstd -c -o $@ $<

# tbl2asn

# Add structured comments to the FASTA file.
%.fsa: %.fa
	sed '/^>/s/ .*/ [organism=Marmota flaviventris] [location=mitochondrion] [completeness=complete] [topology=circular] [mgcode=2]/' $< >$@

# Convert annotations from Genbank TBL to GBF.
marmot-mt.mitos.fix.gbf: %.gbf: %.tbl marmot-mt.fsa marmot-mt.cmt marmot-mt.sbt
	ln -sf marmot-mt.fsa $*.fsa
	tbl2asn -a s -i $*.fsa -w marmot-mt.cmt -t marmot-mt.sbt -Z $*.discrep -Vbv
	mv errorsummary.val $*.errorsummary.val

# Convert annotations from Genbank TBL to GBF.
%.gbf: %.tbl %.fsa %.cmt marmot-mt.sbt
	tbl2asn -a s -i $*.fsa -w $*.cmt -t marmot-mt.sbt -Z $*.discrep -Vbv
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

# minimap2

%.$(ref).paf.gz: $(ref).fa %.fq.gz
	minimap2 -t$t -c -xmap-ont $^ | $(gzip) >$@

# R

# Plot a PAF file.
%.paf.gz.pdf: %.paf.gz
	Rscript -e 'rmarkdown::render("plot-paf.rmd", "html_document", "$*.plot-paf.html", params = list(input_paf="$<"))'
