# Molecular Phylogenetics Guide
This is a practical guide for identifying homologous genes and constructing phylogenetic trees using command line tools. 

This guide is intended to be from basics, and is aimed at students working with plant genomes - but the workflow is applicable to **any gene family**. This guide also assumes a Linux environment (`Ubuntu/WSL` or `MobaXTerm` on Windows), although most commands should also work on the MacOS terminal.
This pipeline is used regularly by David Hoey.

# Software installation
## Very basic bash command line
- Suppose input file, e.g., `input.txt`, is placed in `C:\Users\yourname\DATA`
- This folder can be accessed from Ubuntu as `/mnt/c/Users/yourname/DATA`
- In short,  `C:` on Windows = `/mnt/c` on Ubuntu. 
- `C:\Users\david\Documents\Applications` becomes
```
cd /mnt/c/users/david/documents/applications
```
basic commands that make life easier:
```
rm file.name
```
Use above to delete files in case they are incorrectly installed.
```
ls
```
Use above to see environment - i.e. the contents of the directory that you are in.
It means "list", but I like to remember it as "let's see" (!)

## Installations (conda/homebrew-based)
Install the following programs before beginning:

- MAFFT (multiple sequence alignment)
- trimAL (alignment trimming)
- IQ-TREE2 (phylogenetic inference)
- HMMER (homology searches)
- SeqKit (sequence extraction)
- EMBOSS (translation of CDS)

The easiest method is via Anaconda. Install Miniconda - either download directly from site and enter below code in directory containing file:
```
bash Miniconda3-latest-Linux-x86_64.sh
```
Or input below code directly
```
curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
sh Miniconda3-latest-Linux-x86_64.sh
```
Add channel configurations (bioconda contains `HMMER`, `mafft`, `trimAL`, and `iqtree`)
```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --add channels biocore
```
Install necessary functions with the following command:
```
conda install mafft trimal iqtree hmmer seqkit emboss pal2nal
```
Other important bioinformatics programs, such as blast, can also be installed with the `conda install` command.
```
conda install blast
```

# Phylogenetics pipeline

## Typical phylogenetics workflow
Most gene family analyses follow the same basic workflow:
```
Identify homologs
        │
        ▼
Collect sequences
        │
        ▼
Curate dataset
        │
        ▼
Multiple sequence alignment
        │
        ▼
Trim poorly aligned regions
        │
        ▼
Infer phylogenetic tree
        │
        ▼
Interpret evolutionary relationships
```

# 1. Finding homologs
There are several ways to collect sequences. For small datasets, sequences can be downloaded from databases such as:
- JGI Phytozome (https://phytozome-next.jgi.doe.gov/)
- PlantTFDB (https://planttfdb.gao-lab.org/family.php?fam)
- NCBI (https://www.ncbi.nlm.nih.gov/)
- SYMDB (https://www.polebio.lrsv.ups-tlse.fr/symdb/)
(Last accessed July 2026).
However, manually collecting sequences is time-intensive and rapidly becomes impractical.

For a robust method for determining homologs in a genome, HMMER can be used. Most useful is the `phmmer` command, which can be used to find homologs of a query in a given proteome database:
```
phmmer sequence_query.fasta proteome.fasta > query_phmmer_output.txt
```
However, for projects involving many genomes, I recommend using **PhyloDig**, which automates homology searches, optional Pfam filtering, and sequence extraction.

# 2. Sequence curation
Good phylogenies depend more on careful sequence selection than on sophisticated tree-building methods. 
Some general recommendations:
- It is best to get representatives from a broad selection of species, the less over-representation of one group, the better.
- If a particular genome is poorly annotated, this can cause problems with trimming (usually too much trimming if sequences are too divergent). Therefore, make sure all genomes being included are of a decent quality.
- Ideally have one representative splice isoform per gene.
- Check that predicted proteins contain the expected conserved domains (`PhyloDig` Pfam filtering can do this).
- Small trees can run fairly quickly on a fairly standard laptop or computer. Large phylogenies (>1000 sequences) may require running overnight or on a computer cluster.

## Combining fasta files
For producing alignments and trees you will need to know how to combine `.fasta` files.
The easiest way to do this is to use the `cat` command, to concatenate your files (i.e. to append them one after another).
To concatenate all `.fasta` files in a directory into a single file:
```
cat *.fasta > combined.fasta
```
The following script can be used for combining all fasta files in a directory, for a large phylogeny. I find this very helpful when playing around with my species set. But if you are happy with using the above command instead, go ahead and do that. I prefer the below as it gives you some control over directories.
```
#!/bin/bash

input_dir="path/to/your/input/directory"
output_file="combined.fasta"

# Clear the output file if it exists
> "$output_file"

for input_file in "$input_dir"/*.fasta; do
    cat "$input_file" >> "$output_file"
    echo "Appended $input_file to $output_file"
done

echo "All files combined into $output_file"
```
**Make sure that each input file has a return line at the bottom, or fasta files may be combined improperly.**

# 3. Multiple sequence alignment
To make your gene tree, you need to align your sequences using MAFFT or a similar multiple sequence alignment software.
For most datasets the following is sufficient:
```
mafft input.file > output.file
```
For larger or more divergent gene families, I usually obtain better results using the following (takes longer):
```
mafft --localpair --maxiterate 1000 input.fasta > alignment.fasta
```
# 4. Alignment trimming
You **can** manually trim, although I never got the hang of this. My preferred method is using `TrimAL`:
```
trimal -in input_align.fasta -out output_trim.fasta -fasta -gappyout
```
You can also try the automated settings or other options.
```
trimal -in input.file -out output.file -fasta -automated1
```
With the gap threshold `-gt` option, you can manually set how much you want to trim off. 
For example, `-gt 0.8` will remove a column if >20% are gaps.
```
trimal -in input.file -out output.file -gt 0.8
```
No trimming method is universally optimal. Always inspect the alignment if the results appear unexpected. There are other trimming softwares out there such as `phyx` but I mainly use trimAL 

# 5. Phylogenetic inference
Once you are happy with your trimming and sequence selection, you can now produce your gene tree using IQ-TREE2.

IQ-TREE is the biggest bottleneck so do curate your sequences carefully before this.

Ensure all sequences desired are included before this step as the other steps can be quite fast & trivial once correctly installed.

-T command controls how many CPU cores iqtree uses. If unsure, run -T AUTO. Sometimes IQ-TREE does not like running short alignments with too many cores, AUTO tests the optimum number of threads. For very large phylogenies, these can be carried out on the HPC.
## 5.a) Protein trees
Trees using protein sequences are the most straightforward to produce, and depending on the gene family, you can get all the information you need from it.
A typical command for a protein alignment:
```
iqtree2 -s output_trim.fasta -m MFP -bb 10000 -ninit 10000 -nm 10000 -T AUTO
```
## 5.b) CDS trees
- Sometimes, domains are small and automated trimming removes too much info from amino acid alignments. 
- It is sometimes necessary to produce an alignment using coding sequences, which preserve 3x more information (3 nucleotides = 1 AA). Furthermore, each nucleotide in a codon evolves at different rates, so different evolutionary models can be applied to each one in order to improve your tree inference.
- You need a `.fasta` file with all of your coding sequences which you would like to align in it, and a matching file with all of the protein sequences using the same sequence headers. (**Hint**: you can use `PhyloDig` to produce this!).
- If you have the coding sequence file you can use `transeq` (EMBOSS) to produce it. Make sure to include the `-trim` option (important for `pal2nal` compatibility).
```
  transeq -sequence cds.fasta -outseq prot.fasta -trim -frame 1
```
- Run mafft normally on the protein file (if you are doing this you probably want the long one):
```
mafft --localpair --maxiterate 1000 prot.fasta > prot_align.fasta
```
- Now you have an alignment file, you can generate a codon alignment using pal2nal as follows:
```
pal2nal.pl prot_align.fasta cds.fasta -output fasta > cds_align.fasta
```
- Now you can trim the codon alignment. I still like `-gappyout` for this but try some other options. In my experience, the longer the post-trim alignment, the better. 
```
trimal -in cds_align.fasta -out cds_trim.fasta -gappyout
```
- Now, you need to make a nexus.part file to help IQ-TREE do the evolutionary model testing.
- Run the following on your trimmed CDS alignment:
```
seqkit stats cds_trim.fasta
```
- Run it also on a protein file (trimmed with the same settings) to check the number of sites. The CDS file should be 3x longer!
- Let's say your trimmed CDS file is 615 in length according to `seqkit stats`. Your nexus file should looks something like this:
```
#nexus
begin sets;
  charset pos1 = DNA, 1-613\3;
  charset pos2 = DNA, 2-614\3;
  charset pos3 = DNA, 3-615\3;
end;
```
- Save it as part.nex
- Now you can run IQ-TREE on it.
```
iqtree2 -s cds_trim.fasta -p part.nex -m MFP+MERGE -bb 10000 -nm 10000 -alrt 10000 -T AUTO
```
This should improve branch lengths and bootstrap values substantially.

## 5.c) Optimising alignments and trimming for better trees (making the best of a bad bunch!)
Making a phylogeny is an iterative process.

Let's say you have done all of the above and your bootstraps are still a bit rubbish, or you don't think you have resolved the tree fully yet.

Try first adding or removing some genomes (perhaps some of the less well established genomes which you added are causing problems with aligning or trimming?).

Sometimes a single mid-quality genome has messed up my alignments completely. Usually the problem is that they are too short after trimming for decent inference.

But - some gene families are indeed just a bit tricky, and require some more systematic analysis of alignments and trimming algorithms.

So, let's do that. First, try a few different algorithms from `mafft`.
```
mkdir align_test
mafft --auto input_protein_combined.fasta > align_test/auto.fasta
mafft --localpair --maxiterate 1000 input_protein_combined.fasta > align_test/linsi.fasta
mafft --globalpair --maxiterate 1000 input_protein_combined.fasta > align_test/ginsi.fasta
mafft --ep 0 --genafpair --maxiterate 1000 input_protein_combined.fasta > align_test/einsi.fasta
```
Once that is finished, try different trimming algorithms on all of these.
```
cd align_test
for f in *.fasta; do
    base=${f%.fasta}
    trimal -in "$f" -out "${base}_gappyout.fasta" -gappyout
    trimal -in "$f" -out "${base}_automated1.fasta" -automated1
done
```
Now you can have a look at the length of each alignment: see which ones are over trimmed, poorly aligned, etc.
```
seqkit stats *_gappyout.fasta *_automated1.fasta
```
In general the alignment which is preserving most of the sequence will produce the better tree. Better still if you use this alignment for a CDS tree. 
Look at this example - only the `--globalpair` alignment will produce a good tree for this gene family!
```
file                    format  type     num_seqs  sum_len  min_len  avg_len  max_len
auto_gappyout.fasta     FASTA   Protein       359   10,052       28       28       28
einsi_gappyout.fasta    FASTA   Protein       361   11,552       32       32       32
ginsi_gappyout.fasta    FASTA   Protein       370  304,880      824      824      824
linsi_gappyout.fasta    FASTA   Protein       360   11,160       31       31       31
auto_automated1.fasta   FASTA   Protein       357    4,284       12       12       12
einsi_automated1.fasta  FASTA   Protein       360    7,920       22       22       22
ginsi_automated1.fasta  FASTA   Protein       370  204,980      554      554      554
linsi_automated1.fasta  FASTA   Protein       359    8,616       24       24       24
```

# 6. Visualising trees

Visualise treefiles with `iTOL` (interactive Tree Of Life). https://itol.embl.de/personal_page.cgi
You can produce some nice images with `iTOL` web server, its most recent iterations are very good even for large trees.
`FigTree` is also a nice alternative although I use `iTOL` more often.

## 6.a) Rooting trees
'iqtree' produces an unrooted tree by default. An unrooted tree shows the relationships between sequences, but does not indicate the direction of evolutionary divergence. If you don't have an appropriate root, you can still interpret the tree with that limitation.

You can root the tree by including an outgroup. An outgroup is a sequence or a group of sequences which is understood (independently of the tree you are constructing) to fall outside of the group you are investigating. The root is placed on a branch separating the outgroup from the ingroup. The outgroup should be sufficiently closely related that it can be aligned reliably, but sufficiently outside the ingroup that it provides a meaningful rooting point. Your outgroup could be a related group of sequences outside your group of interest, or it can be an homologous sequence from a taxonomically distinct lineage.

```
                    ┌── Gene A
                ┌───┤
                │   └── Gene B
        ┌───────┤
        │       │   ┌── Gene C
        │       └───┤
        │           └── Gene D
────────┤
        │
        └─────────── Outgroup
```
If a suitable outgroup is not available, you can **midpoint root** your tree. This places the root halfway between the longest path between two sequences in the tree. Midpoint rooting can be useful for visualisation, but the position of the root is an assumption rather than an independently supported evolutionary hypothesis. In particular, do not interpret the lineage closest to a midpoint root as the oldest or most ancestral lineage.

In `iTOL`, you can root a tree by selecting the appropriate branch or clade and use the rooting options. If an outgroup has been included, root the tree on the branch separating the outgroup from the ingroup.

## 6.b) Interpreting your tree
### Bootstrap support
Make sure to include bootstrap values in all phylogeny figures! Whether they are shown numerically or in other representation, I don't really mind.

Bootstraps are usually given as percentages and indicate how frequently a particular branch or relationship is recovered across bootstrap replicates. For example, a value of 95 means that the relationship was recovered in 95% of the bootstrap replicates. 

A poorly supported node means that the analysis does not provide strong evidence for that particular relationship.
If bootstraps are low, you may still be able to use the tree, but be very careful with your interpretation of low support nodes.

Rooting and branching support answer different questions. Branch support describes how consistently a particular relationship is recovered by the chosen support method, whereas rooting determines how the resulting tree is interpreted. Remember that **all** nodes can be rotated without changing the relationships represented - this is important for understanding your tree. The left-to-right or top-to-bottom order of branches does not, by itself, indicate evolutionary relationships.

### Branch lengths
Branch lengths represent the amount of evolutionary change inferred along each branch (**not** the passage of time). Very long branches can indicate rapidly evolving sequences, highly divergent sequences, truncated or incorrectly predicted proteins, or alignment problems. Long branches can make inference more difficult, and especially long branches should be investigated rather than interpreted as genuine evolutionary divergence. If your phylogeny has very long or very short branch lengths it is not necessarily a problem, and you can hide the branch lengths from the visualisation if it makes the tree impossible to read. This does not change the underlying tree or analysis, but the resulting figure will show the branching relationships without displaying the inferred amount of evolutionary change.

### Homologs vs orthologs vs paralogs
It is probably a good idea to cover this briefly here, though you can find much more comprehensive explanations elsewhere.

Homologs are sequences that share a common evolutionary origin. Orthologs are a type of homolog that diverged following a speciation event, whereas paralogs diverged following a gene duplication event. All nodes of a tree are either a duplication event, or a speciation event - however distinguishing these requires comparison with the species relationships and other evidence. 

The closest sequence to your gene in a tree is not necessarily its ortholog. Gene duplications can produce multiple paralogs within the same species, meaning that orthologous relationships can be more complicated than simply finding the most similar sequence. Lineage-specific duplication or proliferation is quite common and paralogous genes might have different functions within those lineages.

# Notes
## Common pitfalls
The most common causes of bad phylogenies are:
- Missing important taxa
- Including truncated or mis-annotated proteins
- Poor alignments
- Over-trimming or under-trimming

When a tree looks unusual, the problem is most likely input data rather than the phylogenetic algorithm.

## Other useful commands
Proteomes/genomes which are compressed may be in gunzip format (.gz), these can be unzipped with the following command:
```
gunzip file.gz
```
For extracting sequences by gene identifier:
```
seqtk subseq input.fasta name.list > output.fasta
```

## Notes on the use of AI
Large language models can be very helpful in troubleshooting bioinformatic pipelines, and I have used them for this purpose often. But, do not be overly reliant on them: they have some major limitations particularly when building large or complex pipelines. Examine and interpret your trees yourself and be skeptical of low support branches. 

## Further resources
- PhyloDig repository: automated database mining for comparative phylogenomics. https://github.com/davidjhoey/phylodig
- HMMER documentation. https://github.com/EddyRivasLab/hmmer, cite: doi.org/10.1371/journal.pcbi.1002195
- MAFFT documentation. https://mafft.cbrc.jp/alignment/software/manual/manual.html, cite: doi:10.1093/molbev/mst010
- PAL2NAL documentation. https://github.com/liaochenlanruo/PAL2NAL, cite: doi.org/10.1093/nar/gkl315
- IQ-TREE documentation. https://iqtree.github.io/doc/, cite: doi.org/10.1093/molbev/msaa015
- TrimAL documentation. https://trimal.readthedocs.io/en/latest/, cite: doi.org/10.1093/bioinformatics/btp348

## Referencing this page
Hoey DJ. plantphylo: Molecular Phylogenetics Guide. GitHub repository: https://github.com/davidjhoey/plantphylo
