#! /bin/bash

#Input checks

#check if the the argument is less than 1
# $# is the number of arguments(input) passed to the script; so if $# is less than 1(=0), echo..
# $1 is first argument passed to the script.
if [ "$#" -lt 1 ]; then
    echo -e "\nusage: ./analyze_GFF_features.sh < chromosomeID >\n"
    exit 1
fi

#Validity check
chromosomeID="$1"
array=("X", "Y", "MT", 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23)

#Not valid #
if [[ ! ${array[@]} =~ $chromosomeID ]]; then

    echo -e "\n"$chromosomeID" is not a valid chromosomeID (possible values : 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 MT X Y )\n"
    exit

fi
#Retrieving the data
    #Check if file exist: "test -f" tests file exist and regular or not
if ! test -f ./Homo_sapiens.GRCh38.110.chromosome."$1".gff3; then  #[ -e .//Homo_sapiens.GRCh38.110.chromosome."$1".gff3 ] also works
       
    #Download file
    curl -s https://ftp.ensembl.org/pub/current_gff3/homo_sapiens/Homo_sapiens.GRCh38.110.chromosome."$1".gff3.gz | gunzip -c > Homo_sapiens.GRCh38.110.chromosome."$1".gff3
fi
    #-s: silent, -c: the output should be written to the standard output (stdout) rather than replacing the input file
    #get a single file without giving output name: -O(Capital;not zero)

####################################

# Feature count
echo "Feature count chromosome "$1":"
echo "--------------------------"

#Display and count all the different types of features 
gff3=Homo_sapiens.GRCh38.110.chromosome."$1".gff3
features=$(awk '$1!~/^#/ {print $3}' $gff3)
echo "$features" | sort -f | uniq -c | awk '{printf "%-5s %s\n", $1, $2}'

echo " "

#####################################

# Top 10 lists

#exon

#Count exons
grep exon "$gff3" | cut -f 7-9 | cut -f 1 -d ';' > exons_"$1".txt # grep exon: Filter for exon and display line that contain exon
np=$(grep '+' exons_"$1".txt | wc -l)
nm=$(grep '-' exons_"$1".txt | wc -l)
nt=$(expr $np + $nm)

echo "Top 10: chromosome "$1":"
echo -e "--------------------------\n"

#Determine top 10 transcripts
echo -e ">>>>transcriptIDs with the highest number of exon"
#make a lit of top 10 transcripts
cut -f 3 exons_"$1".txt | uniq -c | sort -k 1 -n -r | head -10 | cut -f 2 -d ':' > top10_exon_id_"$1".txt

#Display top 10 transcripts
while read id_exon;
do
    
    n_exon=$(grep -c "$id_exon" exons_"$1".txt) #The number of exon
    gene=$(grep "ID=transcript:$id_exon" "$gff3" | cut -f 9 | cut -f 2 -d ';' | cut -f 2 -d ":") #get a gene of transcript
    discription=$(grep "ID=gene:$gene" "$gff3" | cut -f 9 | cut -f 4 -d ';' | cut -f 2 -d "=") #get a discription from gene line
    echo "Transcript "$id_exon">>> #exon: "$n_exon"  gene:"$gene"  "$discription""
done <top10_exon_id_"$1".txt

rm exons_"$1".txt
rm top10_exon_id_"$1".txt
echo " "


#CDS

grep CDS "$gff3" | cut -f 9 | cut -f 2 -d ';' > CDS_"$1".txt

echo -e ">>>>transcriptIDs with the highest number of CDS"   
cut -f 1 CDS_"$1".txt | uniq -c | sort -k 1 -n -r | head -10 | cut -f 2 -d ':' > top10_CDS_id_"$1".txt

while read id_CDS;
do
    n_CDS=$(grep -c "$id_CDS" CDS_"$1".txt)
    gene_CDS=$(grep "ID=transcript:$id_CDS" "$gff3" | cut -f 9 | cut -f 2 -d ';' | cut -f 2 -d ":")
    discription_CDS=$(grep "ID=gene:$gene_CDS" "$gff3" | cut -f 9 | cut -f 4 -d ';' | cut -f 2 -d "=")
    echo "Transcript "$id_CDS">>> #CDS: "$n_CDS"  gene:"$gene_CDS"  "$discription_CDS""
    
done <top10_CDS_id_"$1".txt

rm CDS_"$1".txt
rm top10_CDS_id_"$1".txt
echo " "


#five-prime-UTR
grep five_prime_UTR "$gff3" | cut -f 9 > fpU_"$1".txt

echo -e ">>>>transcriptIDs with the highest number of five-prime-UTR"   
cut -f 1 fpU_"$1".txt | uniq -c | sort -k 1 -n -r | head -10 | cut -f 2 -d ':' > top10_fpU_id_"$1".txt

while read id_fpU;
do
    n_fpU=$(grep -c "$id_fpU" fpU_"$1".txt)
    gene_fpU=$(grep "ID=transcript:$id_fpU" "$gff3" | cut -f 9 | cut -f 2 -d ';' | cut -f 2 -d ":")
    discription_fpU=$(grep "ID=gene:$gene_fpU" "$gff3" | cut -f 9 | cut -f 4 -d ';' | cut -f 2 -d "=")
    echo "Transcript "$id_fpU">>> #fpU: "$n_fpU"  gene:"$gene_fpU"  "$discription_fpU""
    
done <top10_fpU_id_"$1".txt

rm fpU_"$1".txt
rm top10_fpU_id_"$1".txt
echo " "


#three-prime-UTR
grep three_prime_UTR "$gff3" | cut -f 9 > tpU_"$1".txt

echo -e ">>>>transcriptIDs with the highest number of three-prime-UTR"   
cut -f 1 tpU_"$1".txt | uniq -c | sort -k 1 -n -r | head -10 | cut -f 2 -d ':' > top10_tpU_id_"$1".txt

while read id_tpU;
do
    n_tpU=$(grep -c "$id_tpU" tpU_"$1".txt)
    gene_tpU=$(grep "ID=transcript:$id_tpU" "$gff3" | cut -f 9 | cut -f 2 -d ';' | cut -f 2 -d ":")
    discription_tpU=$(grep "ID=gene:$gene_tpU" "$gff3" | cut -f 9 | cut -f 4 -d ';' | cut -f 2 -d "=")
    echo "Transcript "$id_tpU">>> #tpU: "$n_tpU"  gene:"$gene_tpU"  "$discription_tpU""
    
done <top10_tpU_id_"$1".txt
rm tpU_"$1".txt
rm top10_tpU_id_"$1".txt
echo " "
    

