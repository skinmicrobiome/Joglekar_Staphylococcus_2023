# greedy_pangenome_cover
The greedy_pangenome_cover.pl script selects genomes to maximize coverage of the pangenome (starting with core/softcore)

# example
```
# S. hominis
perl -w /data/Segrelab/bwbin/greedy_pangenome_cover.pl Shominis_gene_presence_absence_roary.csv Shom_DM21B05 Shom_DM18D10 Shom_DM21D06 Shom_DM12A09 > Shom_select.txt

# S. epidermidis
perl -w /data/Segrelab/bwbin/greedy_pangenome_cover.pl Sepidermidis_gene_presence_absence_roary.csv > Sepi_select.txt

# S. capitis
perl -w /data/Segrelab/bwbin/greedy_pangenome_cover.pl Scapitis_gene_presence_absence_roary.csv Scap_2C12_ad > Scap_select.txt
```
