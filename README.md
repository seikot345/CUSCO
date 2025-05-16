# CUSCO : <ins>Cu</ins>rator of <ins>S</ins>ingle-<ins>C</ins>opy <ins>O</ins>rthologs

# Installation & introduction
1. how to build a conda environment that can run cusco pipeline.  
Install [Polyphest][Polyphest] first. And install [Astral][Astral] from [Astral github website][AstralHP], Install also [pangene][pangene] using following command. These software should add a directory to the PATH. Then, create new conda environment with environment.yaml. This section describes the installation procedure assuming software version management using Manbaforge (Miniforge).
```
# Create new environment
conda create -n cusco python==3.10.16 pip
conda activate cusco

# Install polyphest
git clone https://github.com/NakhlehLab/Polyphest.git
cd Polyphest
pip install -r requirements.txt

# update environment
git clone https://github.com/seikot345/CUSCO.git
conda env update --file environment.yaml

# Install pangene
curl -L https://github.com/lh3/pangene/releases/download/v1.1/pangene-1.1-bin.tar.bz2|tar jxf -

# Unzip Astral
unzip Astral.5.7.8.zip

# Add PATH of both softwares
export PATH=$PATH:$(pwd)/pangene-1.1-bin/bin_x64-linux:$(pwd)/pangene-1.1-bin/scripts:$(pwd)/Astral:$(pwd)/Polyphest

# Create new environment
conda env create -n cusco -f environment.yaml
```
2. After the environment of analys is set up, prepare a reference protein sequences at current directory and a species directory that contains the whole genome sequences to be used cusco pipline, like following files:
```
./reference.faa
./species/sample1.fna, sample2.fna, ... 
```
# Preparing input for cusco 
This process can be done with "prepare" command if you prepare following files:
```
./cusco.py
./pangene-1.1-bin
./reference.faa
./species/sample1.fna, sample2.fna, ... 

# Prepare input file for cusco  
python cusco.py prepare -t INT -r reference.faa (or -p)
```
"prepare" command requires a reference sequence file name with -r. If you want the amino acid sequence of a single-copy gene CDS in the output, add -p (The default is the nucleic acid sequence). -t is common to all following commands and will get all available cores unless the number of cores is specified.  
After running "pupepare", the following files and directories should be created in the current directory: 
```
./geneCNV.Rtab
./gff/sample1.gff, sample2.gff, ...
./graph.gfa
./input/sample1.fna, sample2.fna, ... (or sample1.faa, sample2.faa, ...)
./paf/sample1.paf, sample2.paf, ...
```

# Making single-copy orthologs fasta and table file
This process can be done  with "curate" command if you prepare input as follows:
```
./cusco.py
./geneCNV.Rtab
./gff/sample1.gff, sample2.gff, ...
./input/sammple1.fna, sample2.fna, ...  

# Make single-copy orthologs fasta and table file
python cusco.py curate -t INT 
```
This process perform extracting single-copy orthologs from Rtab file and trimming twice by gff and fna or faa.  
After running "curate", the following files and directories should be created in the current directory: 
```
./single_copy_genes/gene1.fasta, gene2.fasta, ...
./single_copy_genes.tsv
```

# Making phylogenetic tree of each gene and samples using single-copy orthologs.
This process can be done with "phylo" command if you prepare following files:
```
./Astral
./cusco.py
./single_copy_genes/gene1.fasta, gene2.fasta, ...

# Make phylogenetic tree
python cusco.py phylo -t INT (or -m STR -bs INT)
```
The fasta id must be edited to include only the species name using the edid function described below. You can chose substitution models using -m (default: -m MFP). "phylo" runs 1000 times of ultrafast bootstrap by default. You can specify the number of UFbootstrap with -bs. Specifying 0 will not execute UFbootstrap.  
After running "phylo", the following files and directories should be created in the current directory: 
```
./gene_tree/gene1.bionj, gene1.ckp.gz, gene1.contree, gene1.iqtree, gene1.log, gene1.mldist, gene1.model.gz, gene1.splits.nex, gene1.treefile, gene1.ufboot, gene1.uniqueseq.phy, ...
./species_tree.tree

if bs 0,
./gene_tree/gene1.bionj, gene1.ckp.gz, gene1.iqtree, gene1.log, gene1.mldist, gene1.model.gz, gene1.treefile, gene1.uniqueseq.phy, ...
./species_tree.tree

# Confirm the number of single-copy genes used in constructing the phylogenetic tree
cat species_*_tree.log | grep "Number of gene trees:" | head -n 1
```

# Making a list of candidate marker single-copy orthologs
This process can be done with "marker" command if you prepare following files:
```
./species_tree.tree
./gene_tree/gene1.contree, gene2.contree, ... (or ./gene_tree/gene1.treefile, gene2.treefile, ...)

# Make marker single-copy orthologs list
python cusco.py marker

# The command if you have original species tree like below
./your_original_species_tree.tree
./gene_tree/gene1.contree, gene2.contree, ... (or ./gene_tree/gene1.treefile, gene2.treefile, ...)

python cusco.py marker -sp your_original_species_tree.tree
```
"marker" command makes list file that candidate ranking of useful marker site for additonal samples. Based on the Robinson–Foulds distance normalized to range between 0 and 1 (the closer to 0 the more similar the species trees are), marker genes are determined which gene phylogenetic trees are similar to the species phylogenetic tree.  
After running "marker", the following files and directories should be created in the current directory: 
```
./marker.list
```
# Primer design
This process can be done with "primer" command if you prepare following files:
```
./cusco.py
./gff/sample1.gff, sample2.gff, ...
./single_copy_genes.tsv
./species/sample1.fna, sample2.fna, ... 


# make multi fasta for detect candedate primaer resion
python cusco.py primer -g gene_name -i species_name1 species_name2 ...
```
"primer" command makes a multifasta file to create primers. It requires a reference marker gene name with -g and also requires a reference name of species with -i. 
After running "primer", the following files and directories should be created in the current directory:
```
./primer_design/gene_name_primerdesign.fasta
```
# Pipeline "prepare" to "marker" 
```
python cusco.py pipeline -t INT -r reference.faa (or -p, -bs INT, -sp)
```
# Modify the fasta ID 
curatescg provides the process to change the fasta ID of all fasta file in a directory specific to each argument.
```
# >MP00001 to >sample1_MP00001 of sample1.fasta in input directory
python cusco.py edid -insert
# >sample1_MP00001 to >sample1 of gene1.fasta in single_copy_genes directory
python cusco.py edid -eject
# >sample1 to >sample1_MP00001 of sample1.fasta in single_copy_genes directory using single_copy_genes.tsv
python cusco.py edid -restore
```
# Citation
Seiko, T., Nagasawa, K., Naito, K. (2025). CUSCO: a tool for curating single-copy orthologs and extracting marker genes for phylogenetic tree construction with extra samples. <i>Authorea</i>

[Polyphest]:https://github.com/NakhlehLab/Polyphest
[pangene]: https://github.com/lh3/pangene
[Astral]: https://github.com/smirarab/ASTRAL/raw/master/Astral.5.7.8.zip
[AstralHP]: https://github.com/smirarab/ASTRAL?tab=readme-ov-file
