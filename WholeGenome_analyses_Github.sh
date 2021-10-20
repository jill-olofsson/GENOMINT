### This schript explain how all the sequenced whole genomes have been analysed

#### Clean raw reads with AdaptorRemoval
SAMP=/PATH/to/sample_file
data=/PATH/to/raw_data
clean=/PATH/to/clean_data

cat $SAMP/sample_file.txt | while read line ; do AdapterRemoval --file1 $data/"$line"_L001_R1_001.fastq.gz --file2 $data/"$line"_L001_R2_001.fastq.gz --settings $clean/"$line".settings --output1 $clean/"$line"_R1_truncated.fastq.gz --output2 $clean/"$line"_R2_truncated.fastq.gz --singleton $clean/"$line"_singleton_truncated.fastq.gz --outputcollapsed $clean/"$line"_collapsed.fastq.gz --outputcollapsedtruncated $clean/"$line"_collapsed_truncated.fastq.gz --discarded $clean/"$line"_discarded.fastq.gz --adapter1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG --adapter2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTGNNNNNNGTGTAGATCTCGGTGGTCGCCGTATCATT --trimns --trimqualities --minquality 20 --minlength 25 --maxns 20 --collapse --gzip ; done

####  Plastome analyses #######
	# Map to the chloroplast genome
SAMP=/PATH/to/sample_file
clean=/PATH/to/clean_data
map_cp=/PATH/to/mapped_data
REF_cp=/PATH/to/cp_reference_genome
			# map
cat $SAMP/sample_file.txt | while read line ; do bowtie2 -x ${REF_cp}/Mlongifolia_cp_noIRa.fasta --no-unal -1 ${clean}/"$line"/"$line"_R1_truncated.fastq.gz -2 ${clean}/"$line"/"$line"_R2_truncated.fastq.gz -U ${clean}/"$line"/"$line"_collapsed.fastq.gz,${clean}/"$line"/"$line"_collapsed_truncated.fastq.gz | samtools view -bS - | samtools sort - > ${map_cp}/"$line"_cpMlong.bam ; done
			# clean
				# remove low quality mapping (q20)
cat ${SAMP}/sample_file.txt | while read line ; do samtools view -b -q 20 ${map_cp}/"$line"_cpMlong.bam > ${map_cp}/"$line"_cpMlong_q20.bam ; done
				# remove PCR duplicates
cat ${SAMP}/sample_file.txt | while read line ; do java -jar picard.jar MarkDuplicates I=${map_cp}/"$line"_cpMlong_q20.bam O=${map_cp}/"$line"_cpMlong_q20_rmDup.bam M=${map_cp}/"$line"_cpMlong_Dup.txt ; done

	# Extract multi-samples fasta file
vcf=/PATH/to/vcf_files
fasta=/PATH/to/final_fasta_files
		# call all positions
cat ${SAMP}/sample_file.txt | while read line ; do bcftools mpileup -q 20 -I -B -Ou -f ${REF}/Mlongifolia_cp_noIRa.fasta ${map_cp}/"$line"_cpMlong_q20_rmDup.bam | bcftools call -c -Ov - > ${vcf}/"$line"_samtools_cp_mapMlong.vcf ; done
		# convert to fasta files 
			#last lane is to change the name from ref to sample in each file
cat ${SAMP}/sample_file.txt | while read line ; do bgzip ${vcf}/"$line"_samtools_cp_mapMlong.vcf ; done
cat ${SAMP}/sample_file.txt | while read line ; do bcftools index ${vcf}/"$line"_samtools_cp_mapMlong.vcf.gz ; done
cat ${SAMP}/sample_file.txt | while read line ; do bcftools query -f [%IUPACGT] ${vcf}/"$line"_samtools_cp_mapMlong.vcf.gz > ${vcf}/"$line"_samtools_cp_mapMlong.fasta ; done
cat ${SAMP}/sample_file.txt | while read line ; do  echo "" >> ${vcf}/"$line"_samtools_cp_mapMlong.fasta ; done
cat ${SAMP}/sample_file.txt| while read line ; do sed "1 i\>$line" ${vcf}/"$line"_samtools_cp_mapMlong.fasta > temp && mv temp ${vcf}/"$line"_samtools_cp_mapMlong.fasta ; done
		# change ambiguities to N
cat ${SAMP}/sample_file.txt| while read line ; do bcftools query -f [%IUPACGT] ${vcf}/"$line"_samtools_cp_mapMlong.vcf.gz | sed 's/R/N/g' | sed 's/Y/N/g' | sed 's/S/N/g' | sed 's/W/N/g' | sed 's/K/N/g' | sed 's/M/N/g' | sed 's/B/N/g' | sed 's/D/N/g' | sed 's/H/N/g' | sed 's/V/N/g' > ${vcf}/"$line"_samtools_cp_mapMlong_noAmb.fasta ; done
cat ${SAMP}/sample_file.txt| while read line ; do  echo "" >> ${vcf}/"$line"_samtools_cp_mapMlong_noAmb.fasta ; done
cat ${SAMP}/sample_file.txt | while read line ; do sed "1 i >$line" ${vcf}/"$line"_samtools_cp_mapMlong_noAmb.fasta > temp && mv temp ${vcf}/"$line"_samtools_cp_mapMlong_noAmb.fasta ; done
#		# cat all fasta files togeter for a multi-sample fasta to use in alignments programs
cat ${SAMP}/sample_file.txt| while read line ; do cat ${vcf}/"$line"_samtools_cp_mapMlong.fasta >> ${fasta}/cp_mapMlong.fasta ; done
cat ${SAMP}/sample_file.txt| while read line ; do cat ${vcf}/"$line""$line"_samtools_cp_mapMlong_noAmb.fasta >> ${fasta}/cp_mapMlong_noAmb.fasta ; done

	# Align all plastomes in mafft and build a phylogenetic tree
aln=/PATH/to/alignments
script=/PATH/to/scripts
phyl=/PATH/to/phylogenetics
		#Align
mafft --retree 1 --maxiterate 0 ${fasta}/cp_mapMlong_noAmb_allrefs.fasta > mafft_cp_mapMlong_noAmb_allrefs.fasta
			# convert to aligned fasta to phylip
perl ${script}/Fasta2Phylip.pl ${aln}/mafft_cp_mapMlong_noAmb_allrefs.fasta ${aln}/mafft_cp_mapMlong_noAmb_allrefs.phy
		# Build a raxml bootstrapped tree
raxmlHPC -s mafft_cp_mapMlong_noAmb_allrefs.phy -n mafft_cp_mapMlong_noAmb_allrefs -m GTRGAMMA -e 0.001 -p 12345
raxmlHPC -s mafft_cp_mapMlong_noAmb_allrefs.phy -n mafft_cp_mapMlong_noAmb_allrefs_boot -m GTRGAMMA -e 0.001 -p 12345 -b 12345 -# 100
raxmlHPC -m GTRGAMMA -e 0.001 -p 12345 -f b -t RAxML_bestTree.mafft_cp_mapMlong_noAmb_allrefs -z RAxML_bootstrap.mafft_cp_mapMlong_noAmb_allrefs_boot -n mafft_cp_mapMlong_noAmb_allrefs_bootresults

###########	Whole Genome Analyses	#######################
	# Map to the M. longifolia genome
SAMP=/PATH/to/sample_file
clean=/PATH/to/clean_data
map_wg=/PATH/to/mapped_data
REF_wg=/PATH/to/cp_reference_genome
	# map
cat $SAMP/sample_fie.txt | while read line ; do bowtie2 -x ${REF_wg}/GCA_001642375.1_Mlong1.0_genomic.fna --no-unal -1 ${clean}/"$line"/"$line"_R1_truncated.fastq.gz -2 ${clean}/"$line"/"$line"_R2_truncated.fastq.gz -U ${clean}/"$line"/"$line"_collapsed.fastq.gz,${clean}/"$line"/"$line"_collapsed_truncated.fastq.gz | samtools view -bS - | samtools sort - > ${map_wg}/"$line"_cpMlong.bam ; done
	# clean
		# remove low quality mapping (q20)
cat ${SAMP}/sample_file.txt | while read line ; do samtools view -b -q 20 ${map_wg}/"$line"_Mlong10.bam > ${map_wg}/"$line"_Mlong10_q20.bam ; done
		# remove PCR duplicates
cat ${SAMP}/sample_file.txt | while read line ; do java -jar picard.jar MarkDuplicates -I ${map_wg}/"$line"_Mlong10_q20.bam -O ${map_wg}/"$line"_Mlong10_q20_rmDup.bam -M ${map_wg}/"$line"_Mlong10_Dup.txt ; done
		# remove reads that mapped to the cp
			# find the ID of reads mapping to the cp
cat ${SAMP}/sample_file.txt | while read line ; do samtools view -F 4 ${map_cp}/"$line"_cpMlong.bam | cut -f1 | sort -u > ${map_cp}/"$line"_cpMlong_mapped_reads.txt ; done
			# remove reads using picard
cat ${SAMP}/sample_file.txt | while read line ; do java -jar picard.jar FilterSamReads I=${map_wg}/"$line"_Mlong10_q20_rmDup.bam O=${map_wg}/"$line"_Mlong10_q20_rmDup_nocp.bam READ_LIST_FILE=${map_cp}/"$line"_cpMlong_mapped_reads.txt FILTER=excludeReadList ; done

	# Call SNPs in ANGSD
geno=/PATH/to/geno_files
		# make the bam file list
ls ${map_wg}/*_Mlong10_q20_rmDup_nocp.bam > ${SAMP}/bamfiles.txt

		# call snps on different groups of samples this is just an example
		# call SNPs with a post cutoff set to 0.95

angsd  -b ${SAMP}/bamfiles.txt -ref ${REF}/GCA_001642375.1_Mlong1.0_genomic.fna -out ${geno}/SNP_angsd_post95 -uniqueOnly 1 -remove_bads 1 -trim 0 -C 50 -baq 1 -minMapQ 20 -minQ 36 -minInd "set to fit" -docounts 1 -setMaxDepthInd "set to fit" -gl 1 -domajorminor 1 -domaf 1 -doglf 2 -dopost 1 -SNP_pval 1e-6 -dogeno 5 --ignore-RG 0 -geno_minDepth 3 -geno_maxDepth 50 -postCutoff 0.95

		###############################

###########	Whole Genome Analyses	#######################
	# PCangsd (this is an example)
pcangsd=/PATH/to/output

python pcangsd.py -beagle ${geno}/SNP_angsd_post95.beagle.gz -admix -o ${pcangsd}/SNP_angsd_post95

		# plot the results in R
R
samp<-read.table("samples_metadata.txt", header=T, dec=".", sep="\t")
data <- as.matrix(read.table("SNP_angsd_post95.cov"))
eig <- eigen(data)

pdf("PCangsd.pdf", height=4, width=12)
par(mfrow=c(1,3))
plot(eig$vectors[,1],eig$vectors[,2], col=as.vector(samp[,6]), pch=samp[,8], 
     xlab=paste("PC1 (",round((eig$values[1]/sum(eig$values))*100),"%)", sep=""),
     ylab=paste("PC2 (",round((eig$values[2]/sum(eig$values))*100),"%)", sep=""))
abline(v=0, lty=2, col="grey30")
abline(h=0, lty=2, col="grey30")

plot(eig$vectors[,1],eig$vectors[,3], col=as.vector(samp[,6]), pch=samp[,8], 
     xlab=paste("PC1 (",round((eig$values[1]/sum(eig$values))*100),"%)", sep=""),
     ylab=paste("PC3 (",round((eig$values[3]/sum(eig$values))*100),"%)", sep=""))
abline(v=0, lty=2, col="grey30")
abline(h=0, lty=2, col="grey30")

plot(eig$vectors[,2],eig$vectors[,3], col=as.vector(samp[,6]), pch=samp[,8], 
     xlab=paste("PC2 (",round((eig$values[2]/sum(eig$values))*100),"%)", sep=""),
     ylab=paste("PC3 (",round((eig$values[3]/sum(eig$values))*100),"%)", sep=""))
abline(v=0, lty=2, col="grey30")
abline(h=0, lty=2, col="grey30")
dev.off()

	# Build a phylogenetic tree
		# thin down the bcf file and turn into a phy file
script=/PATH/to/scripts
geno=/PATH/to/geno_files
phy=/PATH/to/phylogenetics

bcftools view -O v ${geno}/SNP_angsd_post95.bcf -o ${geno}/SNP_angsd_post95.vcf
vcftools --vcf ${geno}/SNP_angsd_post95.vcf --thin 1000 --recode --out ${geno}/SNP_angsd_post95_thin1kb
bcftools view -h ${geno}/SNP_angsd_post95_thin1kb.recode.vcf | tail -1 | cut -f 1,2,4,10- > ${phy}/SNP_angsd_post95_thin1kb_header
bcftools query -f '%CHROM\t%POS\t%REF[\t%TGT]\n' ${geno}/SNP_angsd_post95_thin1kb.recode.vcf | cat SNP_angsd_post95_thin1kb_header - > ${phy}/SNP_angsd_post95_thin1kb.tab
perl ${script}/vcf_tab_to_fasta_alignment_JOL.pl -i ${phy}/SNP_angsd_post95_thin1kb.tab > ${phy}/SNP_angsd_post95_thin1kb.fasta
perl ${scropy}/Fasta2Phylip.pl ${phy}/SNP_angsd_post95_thin1kb.fasta ${phy}/SNP_angsd_post95_thin1kb.phy
	# build a raxml bootstrapped tree
		# do a rapid bootstrap
raxmlHPC -s ${phy}/SNP_angsd_post95_thin1kb.phy -n SNP_angsd_post95_thin1kb -m GTRGAMMA -p 12345 -f a -# 100 -x 12345

	# Extract the RNA and align them
		# use PhyloFlash on the two paired and the collapsed sepearat and then use everything in spades
clean=/PATH/to/clean_data
SAMP=/PATH/to/sample_file
			# make the db
perl phyloFlash_makedb_custom.pl --univec_file UniVec/UniVec.fasta --silva_file SILVA_all_ribosomal_PhyloFlash.fasta
			# run PhyloFlash
cat ${SAMP}/sample_file.txt | while read line ; do perl phyloFlash.pl -lib "$line" -read1 ${clean}/"$line"_R1_truncated.fastq.gz -read2 ${clean}/"$line"_R2_truncated.fastq.gz -readlength 150 ; done
cat ${SAMP}/sample_file.txt | while read line ; do perl phyloFlash.pl -lib "$line"_col_trunc -read1 ${clean}/"$line"_collapsed_truncated.fastq.gz -readlength 150 ; done
cat ${SAMP}/sample_file.txt| while read line ; do perl phyloFlash.pl -lib "$line"_col -read1 ${clean}/"$line"_collapsed.fastq.gz -readlength 150 ; done
			# run spades seperatly
cat ${SAMP}/sample_file.txt | while read line ; do mkdir "$line"/spades ; done
cat ${SAMP}/sample_file.txt | while read line ; do spades.py -o "$line"/spades -t 8 -m 20 -k 53,63,73 -1 "$line"/"$line"."$line"_R1_truncated.fastq.gz.SSU.1.fq -2 "$line"/"$line"."$line"_R1_truncated.fastq.gz.SSU.2.fq -s "$line"/"$line"_col_trunc."$line"_collapsed_truncated.fastq.gz.SSU.1.fq -s "$line"/"$line"_col."$line"_collapsed.fastq.gz.SSU.1.fq ; done
			# turn output into sequential
cat ${SAMP}/sample_file.txt| while read line ; do awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' "$line"/spades/scaffolds.fasta > "$line"/spades/scaffolds_seq.fasta ; done
			# remove short scaffolds <300 bp (read length)
cat ${SAMP}/sample_file.txt | while read line ; do awk 'BEGIN {RS = ">" ; ORS = ""} length($2) >= 300 {print ">"$0}' "$line"/spades/scaffolds_seq.fasta > "$line"/spades/scaffolds_seq_300.fasta ; done
			# add sample names
cat ${SAMP}/sample_file.txt | while read line ; do sed "s/NODE/'$line'_NODE/g" "$line"/spades/scaffolds_seq_300.fasta | sed "s/'//g"> temp && mv temp "$line"/spades/scaffolds_seq_300.fasta ; done
			# align within each individual
cat ${SAMP}/sample_file.txt | while read line ; do mafft --retree 1 --maxiterate 0 --adjustdirection "$line"/spades/scaffolds_seq_300.fasta > "$line"/spades/mafft_"$line"_scaffolds_seq_300.fasta ; done
			# build tree
raxmlHPC-SSE3 -s rRNA.phy -n rRNA -m GTRGAMMA -e 0.001 -p 12345
raxmlHPC-SSE3 -s rRNA_pilot.phy -n rRNA_boot -m GTRGAMMA -e 0.001 -p 12345 -b 12345 -# 100
raxmlHPC-SSE3 -m GTRGAMMA -e 0.001 -p 12345 -f b -t RAxML_bestTree.rRNA -z RAxML_bootstrap.rRNA_boot -n rRNA_bootresults

	# NGS admix (10 runs per K)
geno=/PATH/to/geno_files
out=/PATH/to/output

for var in $(seq 1 10);
  do
$NGSadmix -likes ${geno}/SNP_angsd_post95.beagle.gz -K 1 -outfiles ${out}/admix_SNP_angsd_post95_01-"$var"
done
for var in $(seq 1 10);
  do
$NGSadmix -likes ${geno}/SNP_angsd_post95.beagle.gz -K 2 -outfiles ${out}/admix_SNP_angsd_post95_02-"$var"
done
for var in $(seq 1 10);
  do
$NGSadmix -likes ${geno}/SNP_angsd_post95.beagle.gz -K 3 -outfiles ${out}/admix_SNP_angsd_post95_03-"$var" 
done
for var in $(seq 1 10);
  do
$NGSadmix -likes ${geno}/SNP_angsd_post95.beagle.gz -K 4 -outfiles ${out}/admix_SNP_angsd_5post95_04-"$var" 
done
for var in $(seq 1 10);
  do
$NGSadmix -likes ${geno}/SNP_angsd_post95.beagle.gz -K 5 -outfiles ${out}/admix_SNP_angsd_post95_05-"$var" 
done
for var in $(seq 1 10);
  do
$NGSadmix -likes ${geno}/SNP_angsd_post95.beagle.gz -K 6 -outfiles ${out}/admix_SNP_angsd_post95_06-"$var" 
done
for var in $(seq 1 10);
  do
$NGSadmix -likes ${geno}/SNP_angsd_post95.beagle.gz -K 7 -outfiles ${out}/admix_SNP_angsd_post95_07-"$var" 
done
for var in $(seq 1 10);
  do
$NGSadmix -likes ${geno}/SNP_angsd_post95.beagle.gz -K 8 -outfiles ${out}/admix_SNP_angsd_post95_08-"$var" 
done
for var in $(seq 1 10);
  do
$NGSadmix -likes ${geno}/SNP_angsd_post95.beagle.gz -K 9 -outfiles ${out}/admix_SNP_angsd_post95_09-"$var" 
done
for var in $(seq 1 10);
  do
$NGSadmix -likes ${geno}/SNP_angsd_post95.beagle.gz -K 10 -outfiles ${out}/admix_SNP_angsd_post95_10-"$var" 
done

	# plot results (example)
R
samp<-read.table("samples_metadata.txt", header=T, dec=".", sep="\t")
sampA<-samp[order(samp[,3]),]
#
admix_k2<-t(as.matrix(read.table("admix_SNP_angsd_post95_02-01.qopt")))
admix2_k2<-admix_k2[,order(samp[,3])]
pdf("K2.pdf", height=3, width=10)
barplot(admix2_k2,col=c("cyan","orange"),space=0,border=NA,xlab="Individuals",ylab="admixture")
dev.off()
#
admix_k3<-t(as.matrix(read.table("admix_SNP_angsd_post95_03-01.qopt")))
admix2_k3<-admix_k3[,order(samp[,3])]
pdf("K3.pdf", height=3, width=10)
barplot(admix2_k3,col=c("orange","cyan","green"),space=0,border=NA,xlab="Individuals",ylab="admixture")
dev.off()
#
admix_k4<-t(as.matrix(read.table("admix_SNP_angsd_post95_04-01.qopt")))
admix2_k4<-admix_k4[,order(samp[,3])]
pdf("K4.pdf", height=3, width=10)
barplot(admix2_k4,col=c("green","orange","cyan","seagreen"),space=0,border=NA,xlab="Individuals",ylab="admixture")
dev.off()
#
admix_k5<-t(as.matrix(read.table("admix_SNP_angsd_post95_05-01.qopt")))
admix2_k5<-admix_k5[,order(samp[,3])]
pdf("K5.pdf", height=3, width=10)
barplot(admix2_k5,col=c("cyan","coral","seagreen","green","orange"),space=0,border=NA,xlab="Individuals",ylab="admixture")
dev.off()

	# Analyse pairwise (Fst example)
		# convert bcf to vcf
geno=/PATH/to/genotypes
SAMP=/PATH/to/sample_file
fst=/PATH/to/output

bcftools view ${geno}/SNP_angsd_50missing_post95.bcf -O v -o ${geno}/SNP_angsd__post95.vcf
bgzip ${geno}/SNP_angsd_post95.vcf
		# calculate Fst between groups of samples using vcftools
vcftools --gzvcf ${geno}/SNP_angsd_post95.vcf.gz --weir-fst-pop ${SAMP}/group1.txt --weir-fst-pop ${SAMP}/group2.txt --maf 0.05 --out ${fst}/group1_group2_maf0.05 



