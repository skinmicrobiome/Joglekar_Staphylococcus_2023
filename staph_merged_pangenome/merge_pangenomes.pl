use strict;
use Data::Dumper;
use Bio::SeqIO;

#marker characters
our $missing="0"; #print where no value in p/a matrix (should be zero or NA to compute on)
our $nd="?";      #print for missing metadata (names, annotations)
my %info;

###############
# LOAD data
my ($solved,$solcat)=load_solved("solved.txt");
#print Dumper($solved);

# Each species gets their own reformatted->panaroo hash but these could be merged
# since the reformatted are unique
my $scap_lookup=load_alias("scap.fa");
my $sepi_lookup=load_alias("sepi.fa");
my $shom_lookup=load_alias("shom.fa");
#print Dumper($scap_lookup);

# Annotations stored by Panaroo name so we need separate hashes to prevent clashes (e.g., group_xxxx)
#6=COG_category,7=Description, 18=CAZy
my @alist=qw( 6 7 8 11 12 13 18 );
#my $scap_cog=load_nog(7,"/data/Segrelab/conlans/Staph_Diversity/results/pangenome_staph_jan23/emapper_web/Scap_pan_genome.emapper.annotations.tsv");
#my $sepi_cog=load_nog(7,"/data/Segrelab/conlans/Staph_Diversity/results/pangenome_staph_jan23/emapper_web/Sepi_pan_genome.emapper.annotations.tsv");
#my $shom_cog=load_nog(7,"/data/Segrelab/conlans/Staph_Diversity/results/pangenome_staph_jan23/emapper_web/Shom_pan_genome.emapper.annotations.tsv");
my $scap_cog=load_nog2("\t","/data/Segrelab/conlans/Staph_Diversity/results/pangenome_staph_jan23/emapper_web/Scap_pan_genome.emapper.annotations.tsv",@alist);
my $sepi_cog=load_nog2("\t","/data/Segrelab/conlans/Staph_Diversity/results/pangenome_staph_jan23/emapper_web/Sepi_pan_genome.emapper.annotations.tsv",@alist);
my $shom_cog=load_nog2("\t","/data/Segrelab/conlans/Staph_Diversity/results/pangenome_staph_jan23/emapper_web/Shom_pan_genome.emapper.annotations.tsv",@alist);

#print Dumper($scap_cog);

# Presence/Absence stored by Panaroo name so we need separate hashes to prevent clashes (e.g., group_xxxx)
#load presence-absence. $sxxx_genomes is a ref to the array of genome names
my ($scap_pa,$scap_genomes)=load_Rtab("/data/Segrelab/conlans/Staph_Diversity/results/pangenome_staph_jan23/Scap_Jan23/capitis_core_results/gene_presence_absence.Rtab");
my ($sepi_pa,$sepi_genomes)=load_Rtab("/data/Segrelab/conlans/Staph_Diversity/results/pangenome_staph_jan23/Sepi_Jan23/epidermidis_core_results/gene_presence_absence.Rtab");
my ($shom_pa,$shom_genomes)=load_Rtab("/data/Segrelab/conlans/Staph_Diversity/results/pangenome_staph_jan23/Shom_Jan23/hominis_core_results/gene_presence_absence.Rtab");
#print Dumper($scap_genomes);

###############
# Loop through Solved graph and write a merged presence/absence table
# Order: Scap,Sepi,Shom

#header
print join("\t","ID","Graph","subcat",@{$scap_genomes},@{$sepi_genomes},@{$shom_genomes});
print "\t".join("\t","Scap_clust",@alist,"Sepi_clust",@alist,"Shom_clust",@alist);
print "\n";

foreach my $id (sort {$a<=>$b} keys %{$solved})
  {
    my $gene_note="";        #hold panaroo genes/groups for graphs
    my $annotation_note="";  #hold eggnog annotation for graphs
    print "S".sprintf("%04d", $id)."\t".$solcat->{$id};

    #gather p/a data for each graph in the solved table
    my %dat; #reset for each line in solved
    foreach my $k (sort keys %{$solved->{$id}})
      {
         #$k is a gene/cluster/node in the graph; Ãsort forces node to be ordered by Species, then Partition
         #if a second node/gene is encountered for a genome, it is skipped
         #parse systematic/reformatted name:
         #Species_Partition_Index_Genelength
         #partitions: (C)ore, she(L)l, cl(O)ud
         my ($sp,$cat,$gid,$len);
         if ($k=~/([A-Z][a-z]+)\_([A-Z])\_(\d+)\_(\d+)/)
           {
             $sp=$1;
             $cat=$2;
             $gid=$3;
             $len=$4;
             if ($cat!~/^[CLO]$/){die "FATAL: $cat isn't a valid partition code\n"};
           }
         else {die "FATAL: can't parse $k\n"};

         ### CAP ###
         my $name;
         #For graphs, the presence absence, categories, etc... are assigned for
         #each species
         if ($sp eq "Cap" )
          {
            $name=$scap_lookup->{$k};
            if (!exists($scap_pa->{$name})){die "FATAL: $name not found in scap presence absence table"};
            if (!exists($dat{"Cap"}))
              {
                $dat{"Cap"}{"cat"}=$cat;
                $dat{"Cap"}{"pa"}=$scap_pa->{$name};
                $dat{"Cap"}{"name"}=$name;
                if (exists($scap_cog->{$name})){$dat{"Cap"}{"ann"}=$scap_cog->{$name}};
                $info{'cap in solved'}++;
              }
            else {print STDERR "INFO: Graph $id is ".$solcat->{$id}." and $k is ignored by precedence\n";
                  $info{'cap pseudo paralog'}++}
            delete($scap_pa->{$name}); #so we know which ones aren't covered by graphs
          }
         ### EPI ###
         elsif ($sp eq "Epi")
          {
            $name=$sepi_lookup->{$k};
            if (!exists($sepi_pa->{$name})){die "FATAL: $name not found in sepi presence absence table"};
            if (!exists($dat{"Epi"}))
              {
                $dat{"Epi"}{"cat"}=$cat;
                $dat{"Epi"}{"pa"}=$sepi_pa->{$name};
                $dat{"Epi"}{"name"}=$name;
                if (exists($sepi_cog->{$name})){$dat{"Epi"}{"ann"}=$sepi_cog->{$name}};
                $info{'epi in solved'}++;
              }
            else {print STDERR "INFO: Graph $id is ".$solcat->{$id}." and $k is ignored by precedence\n";
                  $info{'epi pseudo paralog'}++}
            delete($sepi_pa->{$name}); #so we know which ones aren't covered by graphs
          }
         ### HOM ###
         elsif ($sp eq "Hom")
          {
            $name=$shom_lookup->{$k};
            if (!exists($shom_pa->{$name})){die "FATAL: $name not found in shom presence absence table"};
            if (!exists($dat{"Hom"}))
              {
                $dat{"Hom"}{"cat"}=$cat;
                $dat{"Hom"}{"pa"}=$shom_pa->{$name};
                $dat{"Hom"}{"name"}=$name;
                if (exists($shom_cog->{$name})){$dat{"Hom"}{"ann"}=$shom_cog->{$name}};
                $info{'hom in solved'}++;
              }
            else {print STDERR "INFO: Graph $id is ".$solcat->{$id}." and $k is ignored by precedence\n";
                  $info{'hom pseudo paralog'}++}
            delete($shom_pa->{$name}); #so we know which ones aren't covered by graphs
          }
         else {die "$sp isn't a valid species\n"};
         #print "\t".join("\t",$k,$cat,$name);
      }

    #output line for this graph, using %dat. Write in alpha order by species
    #codes for individual pangenomes
    print "\t";
    if (exists($dat{"Cap"})){print $dat{"Cap"}{'cat'}} else {print $nd};
    if (exists($dat{"Epi"})){print $dat{"Epi"}{'cat'}} else {print $nd};
    if (exists($dat{"Hom"})){print $dat{"Hom"}{'cat'}} else {print $nd};

    #write presence/absence data
    # This could probably be a single loop over all genomes as below with singles
    foreach my $c (@{$scap_genomes})
      {
        if (exists($dat{"Cap"}) && exists($dat{"Cap"}{'pa'}{$c})){print "\t".$dat{"Cap"}{'pa'}{$c}}
        else {print "\t".$missing};
      }
    foreach my $c (@{$sepi_genomes})
      {
        if (exists($dat{"Epi"}) && exists($dat{"Epi"}{'pa'}{$c})){print "\t".$dat{"Epi"}{'pa'}{$c}}
        else {print "\t".$missing};
      }
    foreach my $c (@{$shom_genomes})
      {
        if (exists($dat{"Hom"}) && exists($dat{"Hom"}{'pa'}{$c})){print "\t".$dat{"Hom"}{'pa'}{$c}}
        else {print "\t".$missing};
      }

    #gene & annotation list
    print "\t";
    #if (exists($dat{"Cap"})){print $dat{"Cap"}{'ann'}.";"} else {print $nd.";"};
    #if (exists($dat{"Epi"})){print $dat{"Epi"}{'ann'}.";"} else {print $nd.";"};
    #if (exists($dat{"Hom"})){print $dat{"Hom"}{'ann'}.";"} else {print $nd.";"};
    if (exists($dat{"Cap"})){print $dat{"Cap"}{'name'}."\t"} else {print $nd."\t"};
    if (exists($dat{"Cap"}) && exists($dat{"Cap"}{'ann'})){print $dat{"Cap"}{'ann'}."\t"} 
      else {for (my $q=0;$q<scalar(@alist);$q++){print $nd."\t"}};

    if (exists($dat{"Epi"})){print $dat{"Epi"}{'name'}."\t"} else {print $nd."\t"};
    if (exists($dat{"Epi"}) && exists($dat{"Epi"}{'ann'})){print $dat{"Epi"}{'ann'}."\t"} 
       else {for (my $q=0;$q<scalar(@alist);$q++){print $nd."\t"}};

    if (exists($dat{"Hom"})){print $dat{"Hom"}{'name'}."\t"} else {print $nd."\t"};
    if (exists($dat{"Hom"}) && exists($dat{"Hom"}{'ann'})){print $dat{"Hom"}{'ann'}."\t"} 
       else {for (my $q=0;$q<scalar(@alist);$q++){print $nd."\t"}};
    print "\n";
  }

print STDERR scalar(keys(%{$scap_pa}))." Scap genes remain after graphs\n";
print STDERR scalar(keys(%{$sepi_pa}))." Sepi genes remain after graphs\n";
print STDERR scalar(keys(%{$shom_pa}))." Shom genes remain after graphs\n";

#print remaining that aren't in graphs
my @all_genomes=(@{$scap_genomes},@{$sepi_genomes},@{$shom_genomes});

### CAP ###
my $n=0;
foreach my $k (keys(%{$scap_pa}))
  {
    print join("\t","Sc".sprintf("%04d", $n),"NA",reverse_lookup_cat($k,$scap_lookup).$nd.$nd);
    $n++;
    foreach my $c (@all_genomes)
      {
        if (exists($scap_pa->{$k}->{$c})){print "\t".$scap_pa->{$k}->{$c}}
        else {print "\t".$missing};
      }
    print "\t".$k."\t";
    if (exists($scap_cog->{$k})){print $scap_cog->{$k}."\t"} else {for (my $q=0;$q<scalar(@alist);$q++){print $nd."\t"}};
    for (my $q=0;$q<=scalar(@alist);$q++){print $nd."\t"}; #pad for other species
    for (my $q=0;$q<=scalar(@alist);$q++){print $nd."\t"}; #pad for other species
    print "\n";
  }

### EPI ###
$n=0;
foreach my $k (keys(%{$sepi_pa}))
  {
    print join("\t","Se".sprintf("%04d", $n),"NA",$nd.reverse_lookup_cat($k,$sepi_lookup).$nd);
    $n++;
    foreach my $c (@all_genomes)
      {
        if (exists($sepi_pa->{$k}->{$c})){print "\t".$sepi_pa->{$k}->{$c}}
        else {print "\t".$missing};
      }
    print "\t";
    for (my $q=0;$q<=scalar(@alist);$q++){print $nd."\t"}; #pad for other species
    print $k."\t";
    if (exists($sepi_cog->{$k})){print $sepi_cog->{$k}."\t"} else {for (my $q=0;$q<scalar(@alist);$q++){print $nd."\t"}};
    for (my $q=0;$q<=scalar(@alist);$q++){print $nd."\t"}; #pad for other species
    print "\n";
  }

### HOM ###
$n=0;
foreach my $k (keys(%{$shom_pa}))
  {
    print join("\t","Sh".sprintf("%04d", $n),"NA",$nd.$nd.reverse_lookup_cat($k,$shom_lookup));
    $n++;
    foreach my $c (@all_genomes)
      {
        if (exists($shom_pa->{$k}->{$c})){print "\t".$shom_pa->{$k}->{$c}}
        else {print "\t".$missing};
      }
    print "\t";
    for (my $q=0;$q<=scalar(@alist);$q++){print $nd."\t"}; #pad for other species
    for (my $q=0;$q<=scalar(@alist);$q++){print $nd."\t"}; #pad for other species
    #had to change the tabs for Shom to make sure we didn't end up with a trailing tab
    print $k;
    if (exists($shom_cog->{$k})){print "\t".$shom_cog->{$k}} else {for (my $q=0;$q<scalar(@alist);$q++){print "\t".$nd}};
    print "\n";
  }

print STDERR Dumper(\%info);

### DONE ###

#########
# SUB

sub load_solved
  {
     #0       C       Cap_C_0001_1059=Epi_C_2710_1038,Epi_C_2710_1038=Hom_C_4010_1059,Hom_C_4010_1059=Cap_C_0001_1059
     my %dat;
     my %catret;
     my $n=0;
     open(INF,"<$_[0]") or die "can't open solved\n";
     while (my $newlin=<INF>)
      {
        chomp($newlin);
        my ($id,$cat,$list)=split(/\t/,$newlin);
        $catret{$id}=$cat;
        my @mem=split(/\,/,$list);
        $n++;

        #TODO: is this the place to add a mini-solver for degenerate graphs?

        foreach my $m (@mem)
          {
            my ($x,$y)=split(/\=/,$m);
            $dat{$id}{$x}=1;
            $dat{$id}{$y}=1;
          }
      }
     close(INF);
     print STDERR "$n clusters loaded into solved core\n";
     return(\%dat,\%catret);
  }

sub load_alias
  {
    my %dat;
    my $seqio_object = Bio::SeqIO->new(-file => $_[0]); 
    while (my $seq_object = $seqio_object->next_seq)
      {
        my $val=$seq_object->desc;
        $val=~s/\s+$//g;
        $dat{$seq_object->display_id()}=$val;
      }
    print STDERR scalar(keys(%dat))." aliases loaded from $_[0]\n";
    return(\%dat);
  }

sub load_nog
  {
    my %dat;
    my ($col,$file)=@_;
    open(INF,"<$file") or die "can't open eggnog\n";
    while (my $newlin=<INF>)
      {
        chomp($newlin);
        if ($newlin=~/^\#/ || $newlin=~/^$/){next};
        my @f=split(/\t/,$newlin);
        $dat{$f[0]}=$f[$col];
      }
    close(INF);
    print STDERR scalar(keys(%dat))." annotation $col loaded from $file\n";
    return(\%dat);
  }

sub load_nog2
  {
    #similar to load nog but grabs multiple fields
    my %dat;
    my ($delim,$file,@col)=@_;
    open(INF,"<$file") or die "can't open eggnog\n";
    while (my $newlin=<INF>)
      {
        chomp($newlin);
        if ($newlin=~/^\#/ || $newlin=~/^$/){next};
        my @f=split(/\t/,$newlin);
        foreach (my $i=0;$i<scalar(@col);$i++)
          {
            if ($i>0){$dat{$f[0]}.=$delim};
            $dat{$f[0]}.=$f[$col[$i]];
          }
      }
    close(INF);
    print STDERR scalar(keys(%dat))." annotation ".join(",",@col)." loaded from $file\n";
    return(\%dat);
  }

sub load_Rtab
  {
    my %dat;
    my ($file)=@_;
    open(INF,"<$file") or die "can't open eggnog\n";
    my $head=<INF>;
    chomp($head);
    my @genomes=split(/\t/,$head);
    my $qq=shift(@genomes);
    while (my $newlin=<INF>)
      {
        chomp($newlin);
        my @f=split(/\t/,$newlin);
        my $gene=shift(@f);
        for (my $i=0;$i<scalar(@genomes);$i++)
          {
            $dat{$gene}{$genomes[$i]}=$f[$i];
          }
      }
    close(INF);
    print STDERR scalar(@genomes)." genomes loaded from $file\n";
    print STDERR scalar(keys(%dat))." genes loaded from $file\n";
    return(\%dat,\@genomes);
  }

sub reverse_lookup_cat
  {
    my ($name,$lookup)=@_;
    #given a prokka name, lookup the standardized name and return prokka category
    foreach my $k (keys %{$lookup})
      {
        if ($lookup->{$k} eq $name)
          {
            if ($k=~/([A-Z][a-z]+)\_([A-Z])\_(\d+)\_(\d+)/)
             {
               return($2);
             }
            else {die "FATAL: error parsing $k during reverse lookup\n"};
          }

      }
    #Gene not in the pangenome reference but in the p/a table won't have a
    #lookup
    return($nd);
  }
