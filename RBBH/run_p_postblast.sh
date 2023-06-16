#edges
perl -w make_edges.pl 0.50 40 CvE.blastp EvC.blastp > edges.txt
perl -w make_edges.pl 0.50 40 EvH.blastp HvE.blastp >> edges.txt
perl -w make_edges.pl 0.50 40 HvC.blastp CvH.blastp >> edges.txt

perl -w make_graph_table.pl edges.txt > graphs.txt

perl -w graph_solver.pl graphs.txt > solved.txt

echo "protein"
echo "all"
grep -E '\sC\s' solved.txt | wc -l
echo "complete"
grep -E '\sC\s' solved.txt | grep -E '\=.+\=.+\=' | wc -l
echo "complete core"
grep -E '\sC\s' solved.txt | grep -E '\=.+\=.+\=' | grep -vE '\_[LO]\_' | wc -l
echo "missing an edge"
cat solved.txt | grep -E '\=.+\=.' | grep 'Hom_' | grep 'Cap_' | grep 'Epi_' | wc -l
