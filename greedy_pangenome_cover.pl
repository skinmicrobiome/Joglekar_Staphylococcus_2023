use strict;

#Reads a presence absence table from roary/panaroo
#Read in each genome as a string of 0/1 for presence/absence

my $skip_singles="yes";  #consider singletons? yes/no
                         #YES=totally ignore singles, yes=include in counts but not scoring
my $weighted="yes";      #weight genes based on how many genomes they appear in? yes/no
my $switch_scoring="no"; #switch from weighted to per-gene counting when shell covered
my $first_genome_col=15; #1-based position of first column of presence/absence data
my $max_rounds=20;       #how many genomes should we select?

my %gene_pa; #hold gene presence absence as hash of strings "111010101..."
my @genes;   #hold a list of gene names in order
my %subgrp;  #what type of gene is it (e.g., core) - indexed by gene name
my %prev;    #how many isolates have each gene - indexed by gene name

if (scalar(@ARGV)==0 || ! -e $ARGV[0])
   {print STDERR "please supply a roary P/A table, supplemented with a subgroup column\n";exit(0)};
my $csv=shift(@ARGV);
my @req;
if (scalar(@ARGV)>0){@req=@ARGV};
open(INF,"<$csv") or die "can't open matrix";

#get the genome list
my $header=<INF>;
chomp($header);
$header=~s/\r//g;
my @header=split(/\,/,$header);
for (my $q=0;$q<($first_genome_col-1);$q++){shift(@header)};

my $l=0;
while (my $newlin=<INF>)
  {
    $l++;
    chomp($newlin);
    $newlin =~ s/\r//g;
    my @f=split(/\,/,$newlin,-1);

    #right here we can skip singletons if we want...
    if ($skip_singles eq "YES" && $f[3]==1){next};

    my $w=scalar(@f);

    push(@genes,$f[0]);   #all genes
    #$subgrp{$f[0]}=$f[1]; #what kind they are (e.g., core); no longer supplied, calculate below
    $prev{$f[0]}=$f[3];   #prevalence, # of genomes they are in

    #calculate group
    my ($frac,$grp)=calc_roary_subgroup($f[3],scalar(@header));
    $subgrp{$f[0]}=$grp;

    #if ($grp_chk ne $f[1]){print "WARN: $f[0] has conflicting groups: ".join(" ",$f[1],$grp_chk,$frac)."\n"};

    #collect presence/absence from genomes
    for (my $i=0;$i<scalar(@header);$i++)
      {
         my $p=$i+($first_genome_col-1); #adjust for all the extra columns
         my $pre=$header[$i];        #expected gene prefix is genome name

         if ($p>=$w){die "FATAL: out of bounds $p >= $w at $l, $f[0]\n"};

         if ($f[$p] =~ /$pre/ || $f[$p] =~ /^[A-Z]{8,8}/){$gene_pa{$header[$i]}.="1"}
         elsif ($f[$p] =~ /refound/){$gene_pa{$header[$i]}.="1"}  #not sure what refound means, assume present
         elsif ($f[$i+14] eq ''    ){$gene_pa{$header[$i]}.="0"}
         else {die "FATAL: can't parse *".$f[$p]."* as a $pre gene presence/absence\n"};
      }
  }

if (scalar(@header) != scalar(keys(%gene_pa))){die "FATAL: header count doesn't match genome count\n"};
print scalar(keys(%gene_pa))." genomes found\n";
if ($skip_singles eq "yes")
  {
    #make clear the difference between yes and YES
    print length($gene_pa{$header[0]})." genes in pangenome (skip singles = $skip_singles, during scoring only)\n";
  }
else
  {
    print length($gene_pa{$header[0]})." genes in pangenome (skip singles = $skip_singles)\n";
  }

print "Selecting genomes... (weighted = $weighted), Switch scoring when shell covered? $switch_scoring\n";
if (scalar(@req)>0){print "Required genomes: ".join(",",@req)."\n"};

#Now, initialize a pangenome coverage vector and figure out the best way to
#merge them

my $coverage="0" x length($gene_pa{$header[0]});

my $r=0;
while ($r<$max_rounds)
  {
    my $r_flag="";  #just a flag to mark required genomes in output
    $r++;
    my $best_d=-1;      #best change in pangenome coverage
    my $best_p=-1;      #best delta prevalence (sum of isolates/gene)
    my $best_genome="";
    my $best_cov="";
    my $best_size=0;
    if (scalar(@req)>0)
      {
        #start with required genomes
        my $r1=shift(@req);
        if (!exists($gene_pa{$r1})){die "FATAL: $r1 is required but not found in matrix\n"};
        ($best_d,$best_p,$best_cov)=merge($coverage,$gene_pa{$r1});
        $best_size=size($best_cov);
        $best_genome=$r1;
        $r_flag="!";
      }
    else
      {
        #pick next best
        foreach my $k (sort keys %gene_pa)
          {
            my ($d,$p,$cov)=merge($coverage,$gene_pa{$k});
            if (($weighted eq "yes" && $p>$best_p) || ($weighted eq "no" && $d>$best_d))
              {
                #save genome resulting in best change in pangenome coverage
                $best_d=$d;
                $best_p=$p;
                $best_cov=$cov;
                $best_size=size($best_cov);
                $best_genome=$k;
              }
          }
      }

    print "round $r: $r_flag$best_genome added $best_d genes, resulting in a covered size of $best_size genes (";
    print sprintf("%.3f", ($best_size/length($gene_pa{$header[0]}))).")\n";
    $coverage=$best_cov; 

    #tally coverage of various gene classes
    my %dat;
    for (my $i=0;$i<scalar(@genes);$i++)
      {
        my $sg=$subgrp{$genes[$i]};
        $dat{$sg."_t"}++;
        if (substr($coverage,$i,1) eq "1"){$dat{$sg."_c"}++};
      }

    if (!exists($dat{"core_c"})){$dat{"core_c"}=0};
    if (!exists($dat{"softcore_c"})){$dat{"core_c"}=0};
    if (!exists($dat{"shell_c"})){$dat{"core_c"}=0};
    if (!exists($dat{"cloud_c"})){$dat{"core_c"}=0};

    print "\tcore\t".$dat{"core_c"}."/".$dat{"core_t"};
    print "\tsoftcore\t".$dat{"softcore_c"}."/".$dat{"softcore_t"};
    print "\tshell\t".$dat{"shell_c"}."/".$dat{"shell_t"};
    print "\tcloud\t".$dat{"cloud_c"}."/".$dat{"cloud_t"}."\n";

    if ($switch_scoring eq "yes" && $dat{"shell_c"}==$dat{"shell_t"}){$weighted="no"};

  }

sub calc_roary_subgroup
  {
    #based on the prevalence and total genomes, return a subgroup (Roary-like)
    my $in=$_[0];
    my $tot=$_[1];
    my $f=$in/$tot;
    if ($f>=0.99){return($f,'core')}
    elsif ($f>=0.95){return($f,'softcore')}
    elsif ($f>=0.15){return($f,'shell')}
    else {return($f,'cloud')};
  }

sub merge
  {
    my $s1=$_[0];
    my $s2=$_[1];
    my $d=0;      #delta in coverage
    my $p=0;      #delta in sum of prevalence
    for (my $i=0;$i<length($s1);$i++)
      {
         #update coverage vector s1 if coverage vector s2 is 1 at i
         if (substr($s1,$i,1) eq '0' && substr($s2,$i,1) eq '1')
           {
             #logic below allows singles to be skipped in scoring
             #but included in coverage counts
             if ($skip_singles =~ /yes/i && $prev{$genes[$i]}==1)
               {
                 #don't use in scoring
               }
             else
               {
                 $d++;
                 $p+=$prev{$genes[$i]};
               }
             substr($s1,$i,1,"1");
           }
      }
    return($d,$p,$s1);
  }

sub size
  {
    my $s=0;        #size of current pangenome
    my $s1=$_[0];
    for (my $i=0;$i<length($s1);$i++)
      {
         if (substr($s1,$i,1) eq '1'){$s++};
      }
    return($s);
  }

