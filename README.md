# Mainly contains two types resoureces: 
#####################################################################################<br>  
#####################################################################################<br>  
####### 1. Databases<br>  
####### 2. Tools<br>  
#####################################################################################<br>  

## 1. Bioinformatics databases
### Brief description of functional genomic data

It is essential to understand the biological context of different types of data. The readers are referred to the figure on https://screen.encodeproject.org/, which provides an intuitive understanding of genomic data and the related biology.<br>



C1: DNA<br>
SNP: the change of a single base in DNA.<br>
INDEL: small insertions or deletions.<br>
CNV: copy number variations. The number of a DNA segment. For example, a gene may have different number of copies in different individuals.<br>
Hi-C: measuring physical interactions among DNA segments.<br>

C2: expression<br>
Gene expression: amount of RNA.<br>
Protein expression: amount of proteins.<br>
Single cell gene/protein expression: expression for single cells. Can be used for identifying cells or cell types where genes are expressed/functional.<br>
FISH: can be used to localize or quantify the expression of genes.<br>

C3: epigenetics<br>
Promoters: the regions where transcription factors bind to in order to turn on a gene.<br>
Bis-seq: measuring DNA methylation level.<br>
ATAC-seq: data from which the open chromatin regions of genes could be identified. The existence of open chromatin regions indicate that the associated genes are active/expressed. The open regions contain DNA motifs for transcription factors.<br>
DNase-seq: another tech for identifying open chromatin regions.<br>
ChIP-seq: there are two types of ChIP-seq methods. One for measuring histone modification, the other for measuring TF binding regions.<br>

C4: networks<br>
Protein-protein interactions: describing physical interactions between proteins.<br>
TF-target gene interaction: describig TF and their target genes.<br>
Functional gene networks: describing functional relationship between genes. Each edge represents the co-functional probability (CFP) between two genes.<br>
QTLs: a DNA variant (often a SNP) that is correlated with a molecular phenotype, including gene expression (eQTL), gene splicing (sQTL), DNA methylation (meQTLs) transcript ratio (trQTL), etc. <br>

C5: gene functions<br>
Gene function: what do gene do? GO, KEGG or MSigDB. <br>
Disease-gene association: disease realted to genes. <br>


## Genomes and gene models
#### Genomes (DNA)
Description: DNA sequences of genomes. <br>
Ref: <br>
Web: http://asia.ensembl.org/info/data/ftp/index.html<br>

#### Transcriptomes (RNA)
Description: RNA sequences of transcriptomes. <br>
Ref: <br>
Web: http://asia.ensembl.org/info/data/ftp/index.html<br>



#### Gene annotation in GENCODE
Description: annotation of exon structures of genes and isoforms. Format: GTF. <br>
Ref: <br>
Web: https://www.gencodegenes.org/human/<br>


#### Proteomes in UniProt 
Description: proteomes. <br>
Ref: <br>
Web: https://www.uniprot.org/proteomes/UP000005640<br>

## Gene functions and pathways

#### Function annotation in Gene Ontology 
Description: gene function. <br>
Ref: <br>
Web: http://current.geneontology.org/products/pages/downloads.html<br>


#### Pathways in KEGG 
Description: pathways. <br>
Ref: <br>
Web:https://www.genome.jp/kegg/pathway.html <br>


#### Gene sets in MSigDB
Description: a comprehensive compendium of gene sets, including data from GO, KEGG, etc. <br>
Ref: <br>
Web: https://www.gsea-msigdb.org/gsea/downloads.jsp   (registration needed)<br>



## Disease-associated gene databases
#### GWAS Catalog
Description: a database of disease-associated SNPs discovered by GWAS<br>
Ref: <br>
Web: https://www.ebi.ac.uk/gwas/docs/file-downloads<br>

#### OMIM
Description: A comprehensive, authoritative compendium of human genes and genetic phenotypes that is freely available and updated daily. <br>
Ref: <br>
Web: https://www.omim.org/downloads<br>

####  DIDA
Description: a database of genes/variants for digenic diseases<br>
Ref: DIDA: A curated and annotated digenic diseases database, NAR, 2016<br>
Web: http://dida.ibsquare.be/explore/<br>

####  SFARI
Description: a database of genes/variants implicated in autism.<br>
Ref: <br>
Web: https://gene.sfari.org/<br>



## Regulatory genomic features

#### Ensembl Regulatory Build
Description: Containing 926,535 cis-regulatory elements (CREs) across 839 cell types in humans and 339,815 cCREs across 157 cell types in mice.<br>
Ref: The Ensembl Regulatory Build, Genome Biology, 2015 <br>
Web: http://asia.ensembl.org/info/data/ftp/index.html <br>


#### ENCODE3
Description: coordinates of regulatory features such as promoters.<br>
Ref:  <br>
Web: https://screen.encodeproject.org/ <br>

#### Tissue/cell type-specific GRNs
Description: a comprehensive resource of 394 cell type– and tissue-specific gene regulatory networks for human <br>
Ref: Marbach et al., tissue-specific regulatory circuits reveal variable modular perturbations across complex diseases, Nat Genet, 2016 <br>
Web: http://regulatorycircuits.org <br>


#### GeneEnhancers
Description: a database of enhancers and their target genes<br>
Ref: Fishilevich, S. et al. GeneHancer: genome-wide integration of enhancers and target genes in GeneCards. Database (Oxford) 2017, bax028<br>
Web: https://www.genecards.org/GeneHancer_version_4-4.<br>



#### RegulomeDB
Description:  <br>
Ref:  <br>
Web: https://regulomedb.org/regulome-search/<br>


#### BOCA
Description: open chromatin by ATAC-seq assay in 14 brain regions and 2 cell types from 5 postmortem human brains.<br>
Ref: Fullard et al., An atlas of chromatin accessibility in the adult human brain. Genome Research, 2018. <br>
Web: https://bendlj01.u.hpc.mssm.edu/multireg/ <br>

#### Tissue-resolved TFBS based on DNAse-seq
Description: Atlas of Transcription Factor Binding Sites from ENCODE DNase Hypersensitivity Data Across 27 Tissue Types.<br>
Ref: Cory C. Funk et al., Atlas of Transcription Factor Binding Sites from ENCODE DNase Hypersensitivity Data Across 27 Tissue Types, Cell Rep, 2020. <br>
Web: http://data.nemoarchive.org/other/grant/sament/sament/footprint_atlas <br>


#### HOCOMOCO database of Transcription factor binding site models
Description: A collection of TFBS models for humans and mouse.<br>
Ref:.<br>
Web: https://hocomoco11.autosome.ru/.<br>


## Gene expresion data
#### GEO
Description: The most comprehensive databases of gene expression (both microarray and RNA-seq, and across many species).<br>
Ref: .<br>
Web: https://www.ncbi.nlm.nih.gov/geo/.<br>



#### GTEx
Description: an ongoing effort to build a comprehensive public resource to study tissue-specific gene expression and regulation. Samples were collected from 54 non-diseased tissue sites across nearly 1000 individuals as of Aug. 2020, primarily for molecular assays including WGS, WES, and RNA-Seq.<br>
Ref: .<br>
Web: https://gtexportal.org/home/datasets.<br>



#### TCGA
Description: Gene expression and genotype data for 33 cancers with a total of 11,315 cases.<br>
Ref: .<br>
Web: https://portal.gdc.cancer.gov/.<br>



#### DEE2
Description: Uniformly processed RNA-seq gene expression data.<br>
Ref: Digital expression explorer 2: a repository of uniformly processed RNA sequencing data. Gigascience. 2019.<br>
Web: http://dee2.io/.<br>


#### scQuery
Description: Containing processed scRNA-seq data from over 500 different studies with over 300 unique cell types.<br>
Ref: A web server for comparative analysis of single-cell RNA-seq data, Nature Commun, 2018 .<br>
Web: https://scquery.cs.cmu.edu/.<br>

#### Tabula Muris
Description: a compendium of single cell transcriptome data from the model organism Mus musculus, containing nearly 100,000 cells from 20 organs and tissues.<br>
Ref:  .<br>
Web: https://tabula-muris.ds.czbiohub.org/.<br>


## Projects/Datasets with multi-omic data
#### Geuvadis project
Description: containing RNA sequencing and genotype data of 462 unrelated human lymphoblastoid cell line samples from the CEU, FIN, GBR, TSI and YRI populations from the 1000 Genomes.<br>
Data: gene expression, SNPs.<br>
Ref: Transcriptome and genome sequencing uncovers functional variation in humans, Nature, 2013 .<br>
Web: https://www.ebi.ac.uk/arrayexpress/experiments/E-GEUV-1/.<br>



##################################################
## 2. Tools
##################################################

## Disease-related
Pascal:  summarize SNP-association P values from GWAS at the level of genes using the Pascal tool


##










