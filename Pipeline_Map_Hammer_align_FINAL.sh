#!/bin/bash -l
REF=/media/A/SHARK/REF/sphyrna/GCA_024679065.1_ASM2467906v1_genomic.fna #Reefrence genome of the great
BAM=/media/A/SHARK/MAPHAMMER/BAM # path to  BAM files 
DEPTH=/media/A/SHARK/MAPHAMMER/DEPTH # path to store depth data
SAF=/media/A/SHARK/MAPHAMMER/FINAL/SAF # path for allele frequencies
IND_SAF=/media/A/SHARK/MAPHAMMER/FINAL/IND_SAF # path for allele frequencies for individuals
SFS=/media/A/SHARK/MAPHAMMER/FINAL/SFS # path for SFS
IND_SFS=/media/A/SHARK/MAPHAMMER/FINAL/IND_SFS # path for SFS for indiiduals (to calculate heterozygosity
POPS=/media/A/SHARK/POPMAPS # path to population maps, each file in the folder contains list of individuals per pop
FQ=/media/A/SHARK/FASTQ # path to fastq files
BL=/media/A/SHARK/MAPHAMMER/BAMLISTS # path to bamlists for each pop for ANGSD
FILTERS=" -uniqueOnly 1 -remove_bads 1  -skipTriallelic 1 -minMapQ 30 -minQ 20 -doHWE 1 -hetbias_pval 1e-4 -hwe_pval 1e-4" # filters as per Laturgie paper
SITES=/media/A/SHARK/MAPHAMMER/FINAL/SITES_SFS # here we will gerenate a list of high confidence sites for every pop from which we can generate the SFS
# Index the genome

# bowtie2-build $REF $REF
# samtools faidx $REF
# bwa-mem index $REF


#cd /media/A/SHARK
#align and get individual depth (sample roughly 1 base per locus
#for i in $POPS/*
# do while read a
#  do /home/marevo/Softwares/bowtie2-2.5.0-linux-x86_64/bowtie2 --no-unal --very-sensitive-local -p $1 -x $REF -U $FQ/${a}.fq.gz 2> $BAM/${a}.stat | samtools view -S -b -@ $1 -q 20 | samtools sort  -@  $1 > $BAM/${a}.bam
#    samtools depth $BAM/${a}.bam | awk 'NR % 70 ==0 ' > $DEPTH/${a}.depth
#    echo $BAM/${a}.bam >> $BL/$(basename $i)
# done < $i
#done

# now use about 44 saamples from two populations to gte some basic sumary QC to select sites to obtain the SFS. We can use two very broadly diverging populations, chagos and indonesia, with Fst of 
# about 0.4. This is useful as we know most alleles will have very different frequencies, and we can use a filter for maximum frequency of heterozygous individuals of 0.5 to remove potential
# We output counts (to get a distribution of read depth) as well as quality scores.  Also, we remove sites with average read depth > 2*average read depth* number of individuals, here 2*11.2*44=985
# we also retain only sites sequenced in at least 25% of individuals (later, we will usemore stringent filtering at the population level, .e 50% of individuals)

angsd -P $1 -bam /media/A/SHARK/MAPHAMMER/cha22_indo22 $FILTERS -maxHetFreq 0.5 -minInd 11 -setMaxDepth 985  -doMajorMinor 1 -doCounts 1 -dumpCounts 1 -doQsDist 1 -doDepth 1 -doMaf 1 -dosnpstat 1 -GL 1 -doPost 2 -doGeno 11 -doGlf 2 -out $SITES/Sites
zcat $SITES/Sites.mafs.gz | cut -f 1,2 | tail -n +2 > $SITES/Sites
angsd sites index $SITES/Sites
#now do angsd, retain only sites sequenced for at least 50% of individuals within pops. Calculate unfolded and folded SFSs, using both the samtool and GATK GL models
for i in $BL/*
 do a=$(cat $i | wc -l) # get number of individuals in each pop
    b=$(echo $a*0.25 | bc | awk '{print int($1+0.5)}') # equal to 50% of individuals (to be used as the minimum individuals seuqenced at depth >1 as filter to retain site)
	echo "retain sites sequenced in $b individuals"
    angsd -P $1 -bam $i -doSaf 1 -doMajorMinor 1 -doCounts 1 -doMaf 1 -dosnpstat 1   $FILTERS  -maxHetFreq 0.5   -minInd $b -sites $SITES/Sites -out $SAF/$(basename $i)_GL1 -anc $REF -GL 1 -ref $REF
    winsfs $SAF/$(basename $i)_GL1.saf.idx -t $1 > $SFS/$(basename $i)_GL1_winSFS.sfs
done

# now generate 2D SFSs for chagos vs indo, ningaloo and south GBR 
winsfs $SAF/chagos_GL1.saf.idx $SAF/ningaloo_GL1.saf.idx -t $1 > $SFS/chagos-ningaloo_GL1_winSFS.sfs
winsfs $SAF/chagos_GL1.saf.idx $SAF/indo_GL1.saf.idx -t $1 > $SFS/chagos-indo_GL1_winSFS.sfs
winsfs $SAF/chagos_GL1.saf.idx $SAF/south_gbr_GL1.saf.idx -t $1 > $SFS/chagos-south_gbr_GL1_winSFS.sfs
winsfs $SAF/chagos_GL1.saf.idx $SAF/rowley_GL1.saf.idx -t $1 > $SFS/chagos-rowley_GL1_winSFS.sfs
winsfs $SAF/chagos_GL1.saf.idx $SAF/north_gbr_GL1.saf.idx -t $1 > $SFS/chagos-north_gbr_GL1_winSFS.sfs
winsfs $SAF/chagos_GL1.saf.idx $SAF/scott_GL1.saf.idx -t $1 > $SFS/chagos-scott_GL1_winSFS.sfs
winsfs $SAF/chagos_GL1.saf.idx $SAF/ningaloo_GL1.saf.idx $SAF/indo_GL1.saf.idx -t $1 > $SFS/chagos-ningaloo-indo_GL1_winSFS.sfs

# now generate indivdual SAF and SFSs to  calculate individual heterozygosity

ls $BAM/*bam > /media/A/SHARK/MAPHAMMER/bamlist


for i in $BL/*
while read a
 do angsd  -doSaf 1 -sites $SITES/Sites -C 50 -baq 1 -minQ 20 -minMapQ 30 -GL 1  -ref $REF  -anc $REF -i $a  -P $1  -out  $IND_SAF/$(basename "$a")_GL1
    winsfs $IND_SAF/$(basename "$a")_GL1.saf.idx -fold 1 > $IND_SFS/$(basename "$a")_GL1.sfs
done < $i




