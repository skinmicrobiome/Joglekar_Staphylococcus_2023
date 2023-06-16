module load seqkit
module load blast

perl -w reformat_pangenome.pl ../pangenome_staph_jan23/Sepi_Jan23/epidermidis_core_results/ 49 Epi > sepi.fa 2>sepi.err
perl -w reformat_pangenome.pl ../pangenome_staph_jan23/Shom_Jan23/hominis_core_results/ 55 Hom > shom.fa 2>shom.err
perl -w reformat_pangenome.pl ../pangenome_staph_jan23/Scap_Jan23/capitis_core_results/ 22 Cap > scap.fa 2>scap.err

seqkit translate -T 11 -M scap.fa --trim > scap.faa
seqkit translate -T 11 -M sepi.fa --trim > sepi.faa
seqkit translate -T 11 -M shom.fa --trim > shom.faa

makeblastdb -in scap.faa -dbtype prot -out scap_p
makeblastdb -in sepi.faa -dbtype prot -out sepi_p
makeblastdb -in shom.faa -dbtype prot -out shom_p

#run blasts
echo "swarm -f swarm_p.sh --module blast --time 10:00:00 -g 32"
