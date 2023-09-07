#!/bin/bash

#############################
#####   Input details   #####
#############################
# Tab separated. See examples:

### Phylogenetic profile ###
### all these OGs belong to the same cluster ###
# Sorted according to first column, where the OG names should appear in order of appearance in the reference.
# Exclude paralogs from the profile.
# One way to prepare this file using the fdog output and NCBI data:
# cat cluster_prots.list | while read i; do cat /PATH/fdog_output/$i/*.phyloprofile | grep "|1" | sed 's/|/\t/g' | cut -f1,4,5 | sed 's/@.*@/./'; done > cluster_prots.profile
# Result:
# OG1	Tax_1	Prot_1
# OG1	Tax_2	Prot_2
# OG1	Tax_3	Prot_3
# OG2	Tax_1	Prot_4
# OG2	Tax_2	Prot_5
# OG2	Tax_3	Prot_6


### GFF dir ###
# This is a directory that collects TAXON_NAME.gff
# TAXON_NAME should be identical as it appears in the phylogenetic profile
# For instance, adjust the following sample script:
# grep "CDS" Tax_1_genomic.gff | grep -v "pseudo=true" | sed -e 's/ID=.*protein_id=//' -e 's/;.*//' | cut -f1,4,5,7,9 > mod_gff/Tax_1.gff
# Replicon 	Start 	End 	Strand 	Protein ID
# JQLE01000001 	141 	959 	+ 	GCA_000746395_00001
# JQLE01000001 	1253 	1987 	+ 	GCA_000746395_00002
# JQLE01000001 	2209 	2589 	+ 	GCA_000746395_00003


### Window size ###
# Number of allowed insertions in the cluster.

### Taxon order ###
# Ordered list in which taxa should be analyzed.
# The reference would ideally be first, followed by taxa on evolutionary distance to it.
# Further implementations will summarize information lineage-wise.
# Example:
# Tax_1
# Tax_2
# Tax_3

### Output name ###
# A file with said name will appear in the two output directories.
#
# Synteny/results.txt
# This shows if orthologs per-taxon are neighboring each other.
# GCF_000737145.1	CAB76983.1	CAB76984.1	CAB76985.1
# GCF_000278665.1		          CAB76984.1  CAB76985.1
# GCF_000069245.1		          CAB76984.1
# GCF_000472005.1	CAB76983.1	CAB76984.1	CAB76985.1
#
# Frequency/results.txt
# This will show, how often we see each of the compositions.
# Note that a cluster composition is accounted only as long as a taxon holds two or more orthologs.
# This means, a single protein is not able to form a cluster.
# Example:
#    2  CAB76983.1	CAB76984.1	CAB76985.1
#    1  		        CAB76984.1  CAB76985.1

################################################################################################
###   Sample command:                                                                        ###
###   bash vicinatrix_v3.sh <profile> <gff_dir> <window_size> <taxon_order> <output_name>    ###
###   bash vicinatrix_v3.sh cluster_prots.profile mod_gff/ 3 taxon_order.txt results         ###
################################################################################################

##### Global variables #####
profile=$(echo "$1")
gff=$(echo "$2")
window=$(echo "$3")
order=$(echo "$4")
out=$(echo "$5")

########################
#####   Warnings   #####
########################

if [[ $(cut -f1,2 $profile | sort | uniq -c | grep -v " 1 " | wc -l) -ge 1 ]] ;
    then echo "Co-orthologs found in profile. Exiting process" ;
    cut -f1,2 $profile | sort | uniq -c | grep -v " 1 " | sed 's/.* //' | while read n k;
        do grep "$n" $profile | grep "$k" ;
        done ;
    exit 1 ; # stops the whole script
    fi

############################
#####   Calculations   #####
############################

##### Making directories #####
if [[ $(ls | grep -w "tmp" | wc -l) -eq 0 ]] ; then mkdir tmp; fi
if [[ $(ls | grep -w "Synteny" | wc -l) -eq 0 ]] ; then mkdir Synteny; fi
if [[ $(ls | grep -w "Frequency" | wc -l) -eq 0 ]] ; then mkdir Frequency; fi

### Effective taxa = taxa with at least one ortholog in a cluster ###
join -t $'\t' <(cat -n $order | sed 's/^ \+//' | awk -F "\t" '{print $2 "\t" $1}' | sort -k1,1) <(cut -f2 $profile | sort -u)   | awk -F "\t" '{print $2 "\t" $1}' | sort -g | cut -f2 > tmp/effective_taxa.txt

### Effective window size = number of genes + allowed insertions - 1 (starting gene) ###
rewind=$(cut -f1 $profile | sort -u | wc -l | awk -v window=$window '{print $0 + window - 1 }')

# 1. Produces a template of the cluster for each taxon
# 2. Gets a starting gene and the rest of the genes go into file rest.txt
# 3. Finds replicon of starting gene in the gff
# 4. Finds adjacent positions of starting gene ( number of genes + window size ) and see how many of the rest are there.
# 5. Selects the next gene as the starting gene, and so on
# 6. If the neigborhood when selecting staring genes is the same, they will be merged.
# 7. The results of one cluster give one or more lines of each taxon and the respective clusters that were found together. Orthologs are assigned to the respective column.

cat tmp/effective_taxa.txt | while read taxon;
    do
    join -t $'\t' -a1 <(cut -f1 $profile | uniq | cat -n | sed 's/ //g' | awk -F "\t" '{print $2 "\t" $1}' | sort -k1,1) <(grep "$taxon" $profile | cut -f1,3 | sort -k1,1)  | awk -F "\t" '{print $2 "\t" $0}' | sort -k1,1 -g | cut -f2,3,4 | sed 's/$/\t/' | cut -f1,2,3 > tmp/full.txt ;
    awk -F "\t" '$3 != "" {print $1 "\t" $3}' tmp/full.txt > tmp/partial.txt;
    cat tmp/partial.txt | while read og start;
        do
        grep -vw "$start" tmp/partial.txt  | cut -f2 > tmp/rest.txt;
        replicon=$(grep "$start" $gff/$taxon.gff | cut -f1);
        cat <(echo "$start") <(grep -a$rewind "$start" $gff/$taxon.gff | grep -f tmp/rest.txt | grep "$replicon" |cut -f5) > tmp/hits.txt;
        join -t $'\t' -a1 <(cut -f1,2 tmp/full.txt | sort -k1,1) <(cat tmp/hits.txt | while read hh; do grep -w "$hh" tmp/full.txt | cut -f1,3; done | sort -k1,1) | sed 's/$/\t/' | cut -f1,2,3 | sed 's/.*\t$//' | cut -f1 | sed 's/$/XXX\t/' | tr -d '\n' | sed -e 's/\t$/\n/' -e 's/XXX//g' -e "s/^/$taxon\t/" ;
        done | uniq ;
    done > Synteny/$out.txt

### Empirical frequency ###
# Of the times that at least 2 candidate cluster members are present, how many of them still neighbor any other?

orph=$(cut -f2 $profile | sort | uniq -c | grep " 1 " | sed -e 's/.* //' -e 's/$/|/' | tr -d '\n' | sed -e 's/|$/\n/')
if [[ $(echo $orph | wc -c) -ne 1 ]] ;  then
     grep -vE "$orph" Synteny/$out.txt | cut -f1 --complement | sort | uniq -c | sed -e 's/ \+//' -e 's/ /\t/' > Frequency/$out.txt
     else
     cut -f1 --complement Synteny/$out.txt | sort | uniq -c | sed -e 's/ \+//' -e 's/ /\t/' > Frequency/$out.txt;
fi
echo "$out"
