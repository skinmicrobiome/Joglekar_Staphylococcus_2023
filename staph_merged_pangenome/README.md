Download this directory, update file locations in merge_pangenome.pl and run:

```perl -w merge_pangenomes.pl > output.txt```

or

```perl -w merge_pangenomes_stable.pl > output.txt```

The only difference in the above 2 commands is that the stable version has sorts that make the output reproducible. Without the sorts, unstable ordering of Perl hashes means merged cluster IDs will change

Files:

* Scap_pan_genome.emapper.annotations.tsv - eggnog annotation of scap.fa
* Sepi_pan_genome.emapper.annotations.tsv - eggnog annotation of sepi.fa
* Shom_pan_genome.emapper.annotations.tsv - eggnog annotation of shom.fa
* capitis_core_results - contains scap pangenome presence/absence file from panaroo
* epidermidis_core_results - contains sepi pangenome presence/absence file from panaroo
* hominis_core_results - contains shom pangenome presence/absence file from panaroo
* merge_pangenomes.pl - version used for manuscript
* merge_pangenomes_stable.pl - same as above but with sort added to species-specific loops. generates a stable output
* scap.fa - scap pangenome, 3187 genes
* sepi.fa - sepi pangenome, 4640 genes
* shom.fa - shom pangenome, 4849 genes
* solved.txt - 3-way complete graphs of species-level cluster representative (used to join species pangenomes). See RBBH code
* staph_merged_pangenome_v2_2023Feb07.tsv - output table from merge_pangenomes.pl
* staph_merged_pangenome_v2_2023Feb07.xlsx - formatted with a key
