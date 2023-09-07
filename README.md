# Vicinatrix
Summarizes conservation of synteny by extracting information from GFF files from a given set of genomes

###   Sample command:   ###

    bash vicinatrix_v3.sh <profile> <gff_dir> <window_size> <taxon_order> <output_name>
    bash vicinatrix_v3.sh cluster_prots.profile mod_gff/ 3 taxon_order.txt results

###   Input details   ###
    Tab separated. See examples:

### Phylogenetic profile ###
    A collection of orthologous groups (OGs) that belong to the same gene cluster, such as a bacterial operon ###
    Sorted according to first column, where the OG names should appear in order of appearance in your reference genome.
    Exclude paralogs from the profile.
    One way to prepare this file using the fdog output and NCBI data:
    
    OG1	Tax_1	Prot_1
    OG1	Tax_2	Prot_2
    OG1	Tax_3	Prot_3
    OG2	Tax_1	Prot_4
    OG2	Tax_2	Prot_5
    OG2	Tax_3	Prot_6

### GFF dir ###
    This is a directory that collects TAXON_NAME.gff
    TAXON_NAME should be identical as it appears in the phylogenetic profile
    For instance, adjust the following sample script from a genome assembly from GenBank:
    grep "CDS" Tax_1_genomic.gff | grep -v "pseudo=true" | sed -e 's/ID=.*protein_id=//' -e 's/;.*//' | cut -f1,4,5,7,9 > mod_gff/Tax_1.gff
    
    Replicon 	Start 	End 	Strand 	Protein ID
    JQLE01000001 	141 	959 	+ 	GCA_000746395_00001
    JQLE01000001 	1253 	1987 	+ 	GCA_000746395_00002
    JQLE01000001 	2209 	2589 	+ 	GCA_000746395_00003

### Window size ###
    Number of allowed insertions in the cluster.

### Taxon order ###
    Ordered list in which taxa should be analyzed.
    The reference would ideally be first, followed by taxa on evolutionary distance to it.

    Example:
    Tax_1
    Tax_2
    Tax_3

### Output name ###
    A file with said name will appear in the two output directories.

    Synteny/results.txt
    This shows if orthologs per-taxon are neighboring each other.
    Protein IDs are substituted for the OG name.
    Inversions and other changes of gene order are not considered. Only if they remain in the same neighborhood given the window size
    
    GCF_000737145.1	CAB76983.1	CAB76984.1	CAB76985.1
    GCF_000278665.1		        CAB76984.1	CAB76985.1
    GCF_000069245.1		        CAB76984.1
    GCF_000472005.1	CAB76983.1	CAB76984.1	CAB76985.1
    
    Frequency/results.txt
    This will show, how often we see each of the compositions. Useful to see the more dominant gene clusters (first column)
    Note that a cluster composition is accounted only as long as a taxon holds two or more orthologs.
    This means, a single protein is not able to form a cluster.
    
    Example:
    2  CAB76983.1	CAB76984.1	CAB76985.1
    1		CAB76984.1	CAB76985.1
