#!/bin/bash
#plasmid checker - written Nov-Jan 2020 - Brad Hart
#
#Instuctions:
# Add assembled fasta files to same folder that you run the script from
# Edit script to give location of downloaded plsdb databases
# As part of NCBI E-utilities, you are required to have an API key in your bash_profile or bashrc. To obtain this key, log into an NCBI account and click on your name. Scroll down to generate the key.
# configure you bashrc or bash_profile to include export NCBI_API_KEY=xxxxxxxxxxxxxxxxxxxxxxxxxx (x's being the api key)
# Edit script to give location of fastq.gz files for the fasta files that are in the folder you run the script from.
# Naming rules apply that require A) the fasta file to be the exact same name of the fastq.gz files before the underscore, also, fastq.gz files must have an 'R' before the number of the pair
# For example, If your fasta file is named: K999k.fasta, then your fastq files must be named K999k_R1.fastq.gz K999k_R2.fastq.gz
# blast-2.7, plsdb blast and mash database, bedtools, samtools, bowtie2, mash, ncbi-eutils, abricate (with latest plasmid finder database), seqtk, bcftools, prokka, gffread,
### fasta.p files are assemblies created with denovo plasmid assembly pipeline which is done separately at the moment. Will look towards integrating it.
## I also remove contigs of size smaller than 1000 bases across all fasta files.
# TODO - find a way to verify hits
# Attempt to integrate mobsuite in some way
#
dbs=/home/harbj019/dbs/plsdb
fastqs=/home/harbj019/fastqs/Enterococci
#
#working=$(cd -P -- "$(dirname -- "$0")" && pwd -P)
samps=$(echo *.fasta | sed 's/.fasta//g')
#
#start of attempt to integrate denovo plasmid assembly script# https://academic.oup.com/nar/article/48/18/e106/5901968
#for i in samps; do
#  mkdir -p ${i}temp
#  cp ${fastqs}/${i}*.fastq.gz ${i}temp/
#  gunzip ${i}temp/${i}_R{1,2}.fastq.gz
#
####

for i in $samps; do
  mkdir -p $i
  seqtk seq -L 1000 ${i}.fasta > ${i}/${i}.fasta
  seqtk seq -L 1000 ${i}.fasta.p > ${i}/${i}.fasta.p
  mkdir -p ${i}/searchresults
done
#
# This step can have the fasta files combined into one file for subsequent searches?
#
#blastn
#
for i in $samps; do
  echo -e "\nRunning blast query for ${i}.fasta"
  blastn -query ${i}/${i}.fasta -task blastn -db ${dbs}/plsdb.fna -out ${i}/${i}.blast -num_threads 20 -num_alignments 20 -evalue 0.01 -perc_identity 90 -qcov_hsp_perc 90 -outfmt 6
  awk '{b[$1]="0"; e[$1]="";if (a[$1,$2]=="0") a[$1,$2]=$12; else {score=a[$1,$2]+$12; a[$1,$2]=score}}END{for (i in b) for (j in a) {split(j,c,SUBSEP); if (c[1]==i && a[j]>b[i]) {b[i]=a[j];e[i]=c[2]}}; for (i in b) print i"\t"e[i]"\t"b[i]}' ${i}/${i}.blast
cat ${i}/${i}.blast | cut -f 2 | head -n 7 | tee -a ${i}/${i}.plasmids >> ${i}/searchresults/${i}.blastsorted
mv ${i}/${i}.blast ${i}/searchresults/
echo -e "\nBlast query for ${i}.fasta complete."
echo -e "\nRunning blast query for ${i}.fasta.p"
blastn -query ${i}/${i}.fasta.p -task blastn -db ${dbs}/plsdb.fna -out ${i}/${i}.blast.p -num_threads 20 -num_alignments 20 -evalue 0.01 -perc_identity 90 -qcov_hsp_perc 90 -outfmt 6
awk '{b[$1]="0"; e[$1]="";if (a[$1,$2]=="0") a[$1,$2]=$12; else {score=a[$1,$2]+$12; a[$1,$2]=score}}END{for (i in b) for (j in a) {split(j,c,SUBSEP); if (c[1]==i && a[j]>b[i]) {b[i]=a[j];e[i]=c[2]}}; for (i in b) print i"\t"e[i]"\t"b[i]}' ${i}/${i}.blast.p
cat ${i}/${i}.blast.p | cut -f 2 | head -n 7 | tee -a ${i}/${i}.plasmids >> ${i}/searchresults/${i}.blastsorted.p
mv ${i}/${i}.blast.p ${i}/searchresults/
echo -e "\nBlast query for ${i}.fasta.p complete."
done
echo -e "\nBlast completed for all samples"
#
#mash
#
for i in $samps; do
  echo -e "\nRunning mash screening of ${i} against the plsdb database"
  mash screen ${dbs}/plsdb.msh ${i}/${i}.fasta -v 0.1 -i 0.90 -w | tee -a ${i}/searchresults/${i}.mashresults > ${i}/${i}.mashout
  cat ${i}/${i}.mashout | sort -k 1 -n -r | awk '{print $5}' | head -n 5 >> ${i}/${i}.plasmido
  cat ${i}/${i}.mashout | sort -k 2 -n -r | awk '{print $5}' | head -n 5 >> ${i}/${i}.plasmido
  mash screen ${dbs}/plsdb.msh ${i}/${i}.fasta.p -v 0.1 -i 0.90 -w | tee -a ${i}/searchresults/${i}.mashresults.p > ${i}/${i}.mashout.p
  cat ${i}/${i}.mashout.p | sort -k 1 -n -r | awk '{print $5}' | head -n 5 >> ${i}/${i}.plasmido
  cat ${i}/${i}.mashout.p | sort -k 2 -n -r | awk '{print $5}' | head -n 5 >> ${i}/${i}.plasmido
  cat ${i}/${i}.plasmido | sort -u >> ${i}/${i}.plasmids
  rm ${i}/${i}.plasmido
  echo -e "\nmash screen for ${i} completed"
done
echo -e "\nmash completed for all samples"
#
#abricate
#
for i in $samps; do
#abricate --db plsdb --threads 24 ${i}/${i}.fasta | cut -f 13 | sed 1d >> ${i}/${i}.plasmids
echo -e "\nQuerying plasmidfinder database utilising Abricate for ${i}\n"
  abricate --db plasmidfinder --threads 24 ${i}/${i}.fasta | cut -f 13 | sed 1d | tee -a ${i}/${i}.plasmids >> ${i}/searchresults/${i}.abricate
  abricate --db plasmidfinder --threads 24 ${i}/${i}.fasta.p | cut -f 13 | sed 1d | tee -a ${i}/${i}.plasmids >> ${i}/searchresults/${i}.abricate.p
echo -e "\nAbricate for ${i} complete"
done
#
echo -e "\nAbricate completed for all samples"
#
# symbolic links to fastq.gz files
#
for i in $samps; do
  ln -s ${fastqs}/${i}_R* $(pwd)/${i}/
done
#
# ensure list is unique and remove NZ prefix
#
for i in $samps; do
  cat ${i}/${i}.plasmids | sort -u | sed 's/NZ_//g' > ${i}/${i}.plasmid
done
#
# make bams
#
for i in $samps; do
  mkdir -p ${i}/refs
while read line; do
  FA=${i}/refs/${line}.fasta
  #GB=${i}/refs/${line}.gbk
  R1=${i}/${i}_R1.fastq.gz
  R2=${i}/${i}_R2.fastq.gz
  name=${i}/refs/${line}${i}
  efetch -db nucleotide -format=fasta -id=$line > $FA
  #efetch -db nucleotide -format=gb -id=$line > $GB
  bowtie2-build $FA $FA
  bowtie2 -p 24 -x $FA -1 $R1 -2 $R2 --no-unal | samtools sort -m 16G > "$name".bowtie.bam
  if [[ $(bedtools genomecov -ibam ${name}.bowtie.bam -bga | awk '$4==0') ]]; then
    zero=$(bedtools genomecov -ibam "$name".bowtie.bam -bga | awk '$4==0 {bpCountZero+=($3-$2)} {print bpCountZero}' | tail -1)
    nonzero=$(bedtools genomecov -ibam "$name".bowtie.bam -bga | awk '$4>0 {bpCountNonZero+=($3-$2)} {print bpCountNonZero}' | tail -1)
    percent=$(bc <<< "scale=6; ($nonzero / ($zero + $nonzero))*100")
  else
    percent=100.00
  fi
    echo -e "\n$i covers $percent percent of plasmid $line" >> ${i}/${i}.results
done < ${i}/${i}.plasmid
done
#
#Calculate matches > 80% and create matched plasmid information file.
#
for i in $samps; do
  cat ${i}/${i}.results | awk '{print $3"\t"$7}' | sed 's/\..*\t/\t/' | sed -r '/^\s*$/d' > ${i}/plink.txt
  while read line; do
  percent=$(echo $line | awk '{print $1}')
  if [[ "$percent" -ge "80" ]]; then
  echo ${line} | awk '{print $2}' >> ${i}/refs.txt
  fi
  done < ${i}/plink.txt
  rm ${i}/plink.txt
  echo -e "Accession\tSeq_length\tTitle" > ${i}/${i}_Matched_Plasmids_Details.txt
  while IFS= read -r -u "$fd_num" line; do
  esearch -db nuccore -query "$line" | efetch -format=docsum | xtract -pattern DocumentSummary -element AccessionVersion Slen Title >> ${i}/${i}_Matched_Plasmids_Details.txt
done {fd_num}<${i}/refs.txt
done
#
#Begin consensus process
#
#Tidy up files and move .bam files for use with consensus file creation
#
for i in $samps; do
  while read line; do
mkdir -p ${i}/consensus/${line}
mv ${i}/refs/${line}.fasta ${i}/consensus/${line}/
mv ${i}/refs/${line}${i}.bowtie.bam ${i}/consensus/${line}/
done < ${i}/refs.txt
rm -rf ${i}/refs/
mkdir ${i}/searchresults
mv ${i}/${i}.{blast,mashout,plasmid,plasmido,plasmids} ${i}/searchresults/
done
#
# Get consensus fastq file - vcfutils.pl is part of bcftools
for i in $samps; do
  while read line; do
    bcftools mpileup -Oz -f ${i}/consensus/${line}/${line}.fasta ${i}/consensus/${line}/${line}${i}.bowtie.bam | bcftools call -c | vcfutils.pl vcf2fq > ${i}/consensus/${line}/${i}${line}.fastq
#
# Convert .fastq to .fasta and set bases of quality lower than 20 to N
    seqtk seq -aQ64 -q20 -n N ${i}/consensus/${line}/${i}${line}.fastq > ${i}/consensus/${line}/${i}${line}.fasta
#
#retrieve genbank file for reference sequence
    efetch --db=nucleotide --format=gb --id=${line} > ${i}/consensus/${line}/${line}.gbk
#
# use prokka to retrieve gene ids and annotate gff file
    prokka --cpus 24 --locustag ${line} --proteins ${i}/consensus/${line}/${line}.gbk --prefix ${i}${line} --coverage 40 --outdir ${i}/consensus/${line}/annotation ${i}/consensus/${line}/${i}${line}.fasta
#
# use gffread to merge gff with consensus along with a few formatting changes:
    gffread ${i}/consensus/${line}/annotation/${i}${line}.gff -G | sed 's/gene_name/geneID/' | gffread - -g ${i}/consensus/${line}/${i}${line}.fasta --table @geneid -W -x ${i}/consensus/${line}/consensus.txt
    cat ${i}/consensus/${line}/consensus.txt | sed 's/>/\n>/' | sed 's/segs:.*\t/Geneid: /' > ${i}/${i}.on.${line}_consensus.txt
    rm -rf ${i}/consensus/${line}/consensus.txt
  done < ${i}/refs.txt
done
#
# to extract gene information from the prokka output:
#!/bin/bash
for i in $samps; do
  while IFS= read -r -u "$fd_num" line; do
    echo -e "**PLASMID**" >> ${i}/${i}_Plasmid_and_Gene_Summary.tsv
    echo -e "Accession\tSeq_length\tTitle" >> ${i}/${i}_Plasmid_and_Gene_Summary.tsv
    grep ${line} ${i}/${i}_Matched_Plasmids_Details.txt >> ${i}/${i}_Plasmid_and_Gene_Summary.tsv
    grep ${line} ${i}/${i}.results >> ${i}/${i}_Plasmid_and_Gene_Summary.tsv
    echo -e "**ANNOTATED GENES**" >> ${i}/${i}_Plasmid_and_Gene_Summary.tsv
      if test -f "${i}/consensus/${line}/annotation/${i}${line}.tsv"; then
        cat ${i}/consensus/${line}/annotation/${i}${line}.tsv | grep -v ypothetical | awk -F "\t" '{if($4 == "") print $7; else print $4"\t"$7}' >> ${i}/${i}_Plasmid_and_Gene_Summary.tsv
      else
        echo -e "****SEQUENCE DID NOT ANNOTATE***\nEither prokka found the gbk incomplete or too small or there were no matching genes" >> ${i}/${i}_Plasmid_and_Gene_Summary.tsv
      fi
    echo -e "\n" >> ${i}/${i}_Plasmid_and_Gene_Summary.tsv
  done {fd_num}<${i}/refs.txt
rm -rf ${i}/*.fastq.gz
rm -rf ${i}/${i}.fasta
rm -rf ${i}/refs.txt
done
echo -e "\nPlasmidBot complete!"
