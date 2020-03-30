#!/bin/perl
# 2020
# Alvaro Rodriguez del Rio (main code) and Juliane C. Dohm (some polishing)
#
# Script designed for visually represent syntenic regions based on orthology relations between
# genes in a single reference scaffold of one species and orthologous genes of another species
# which may appear in multiple scaffolds.

use GD::Simple;
use GD;
use List::MoreUtils qw(uniq);
use Getopt::Long;

# input files / parameters
GetOptions('o=s' => \$input, 'g1=s' => \$file_geneIDlist1, 'g2=s' => \$file_geneIDlist2, 'n=s' => \$suffix);
if ($input eq ''){
	die "No orthology file provided\n";
}

if ($file_geneIDlist1 eq '' || $file_geneIDlist2 eq ''){
	print "No gene files provided, only genes in $input will be displayed\n";
}

# open orthology file
open(INPUT, "< $input") or die "$!\n";
@input = <INPUT>;
close(INPUT);

# get reference scaffolds
foreach $line (@input){
	($ref, $query) = split(/\s/,$line) unless $line =~ m/^#/;
	($ref_scaf) = split(/__/,$ref);
	$ref_scafs{$ref_scaf} = 1;
}
$number_of_refScafs = keys %ref_scafs;


#processing each reference scaffold
foreach $refinput (sort keys %ref_scafs)
{
if ($refinput)
{
$grepinput = $refinput."__";
@orthology = grep(/^$grepinput/,@input);


# file opening
my @file_synteny_relations=@orthology;
chomp @file_synteny_relations;


# declaration of databases
our $colorControl=1; # controls the color assigned to each gene/gene pair
our %DB_geneRef_colour=(); # contains the color of each gene
our %DB_geneNoRef_orthologGenesRef=(); 
our %DB_NRscaffolds_genes=();
our %DB_Allscaff_genes=();

# creating database scaffold_gene{scaff}[gene1, gene2,..., genen] both for reference and non reference species
# if gene names of all the genes are provided in geneID lists, they are used for constructing the database

if ($file_geneIDlist1 eq "" || $file_geneIDlist2 eq ""){
  create_scaffold_DB_from_connections_file();
}
else{
  create_scaffold_DB_from_geneIDlist($file_geneIDlist1);
  create_scaffold_DB_from_geneIDlist($file_geneIDlist2);
}

#fill the databases 
my @orthologGenes=create_DBOrthologs(@file_synteny_relations);
create_dataBase_geneRef_color(@file_synteny_relations);
my @orderNRscaff = get_order_NR_scaff(@file_synteny_relations);
create_ortholog_scaffold_DB(@orthologGenes);
my $numberOfScafs=scalar keys(%DB_NRscaffolds_genes);


#geometry variables
our $img;
our $stepref;
our %lenScaff;
our $ymax;
our $ymaxRef;
our $yminRef;
if ($numberOfScafs<6){
  $ymax=1000;
  $ymaxRef=900;
  $yminRef=30;
}
else{
  $ymax=150*$numberOfScafs;
  $ymaxRef=$ymax-$ymax/10;
  $yminRef=$ymax-$ymax*9/10;
}
our $xmax=1000;
our $xmaxRef=300;
our $xminRef=100;
our $xminNoRef=600;
our $xMaxNoRef=800;
my $manage_NR_scaff_positions=0; #controls the position of the non reference scaffolds


# inizializing $img variable, where all the features of the image are saved during the program
$img = GD::Simple -> new($xmax,$ymax);
$img -> font('Arial');

# generating a vector of colors to be used for coloring genes
my @vectorColors;
for (my $i = 0;$i<10000;$i++){
  push(@vectorColors,$img->colorAllocate(rand()*1000,rand()*1000,rand()*1000));
}

# getting ref scaffold name
my $scaffoldRef=$refinput;

###########################
# draw the figure
###########################

#draw reference scaffold
draw_scaffold($scaffoldRef,$xminRef,$xmaxRef,$yminRef,$ymaxRef,0);
draw_guide_lines($scaffoldRef);

#draw non reference scaffolds
get_proportions_orthologs_scaffolds();#creates a DB with the relative size of the scaffold compared to the others
for (my $i = 0; $i < scalar @orderNRscaff; $i++){
  my $proportion=$lenScaff{@orderNRscaff[$i]};
  my $division=($ymax-15*$numberOfScafs)*$proportion;
  my $y1=$manage_NR_scaff_positions+15;
  #my $y2=$y1+$division;
  my $y2=$y1+$division-15/$numberOfScafs;	#space between last scaffold and bottom
  $manage_NR_scaff_positions=$y2;
  draw_scaffold(@orderNRscaff[$i],$xminNoRef,$xMaxNoRef,$y1,$y2,1);
}

#generating a file with the complete figure
my $out;
if ($suffix eq ''){
	open ($out, ">$refinput.png") or die "Could not write \"$refinput.png\"\n";
	print "Synteny image written to \"$refinput.png\".\n";
}
else {
	open ($out, ">$refinput.$suffix.png") or die "Could not wite \"$refinput.$suffix.png\"\n";
	print "Synteny image written to \"$refinput.$suffix.png\".\n"
	
}

binmode $out;
print $out $img->png;

}
}	#end of foreach loop

########################
#subroutines
########################

# for each reference gene, generate a database gene_color{gene}=color.
# orthologs have the same color.
# if a non reference gene has more than one ortholog in the reference scaffold,
# this database will contain all the colors of the orthologs.
sub create_dataBase_geneRef_color{
  my @file=@_;

  #go through the file to build a color database
  for (my $i=0;$i<scalar@file;$i++){
    my $gene=extract_info($file[$i],"g1");
    my $orthologousGene=extract_info($file[$i],"g2");
    $colorControl++;

    #filling the database
    if ($gene ne "-"){
      if ($orthologousGene ne "-"){#there is an ortholog gene for the gene in the refernece scaffold in this line
        if (!defined $DB_geneRef_colour{$gene}){#no color for gene in the reference scaffold
          push @{$DB_geneRef_colour{$gene}},$colorControl;
          push @{$DB_geneRef_colour{$orthologousGene}},$colorControl;
        }
        else{#there is already a color for that gene in the reference scaffold
          push @{$DB_geneRef_colour{$orthologousGene}},@{$DB_geneRef_colour{$gene}}[0];
        }
      }
      else{#this gene does not have an ortholog
        push @{$DB_geneRef_colour{$gene}},"black";
      }
    }
    else{#the line only contains a gene in the non reference scaffold
      push @{$DB_geneRef_colour{$orthologousGene}},"black";
    }
  }
}


#sub to order NR scaffolds so that they appear ordered in the image
sub get_order_NR_scaff{
  my @file = @_;
  my $prevScaffNR;
  my @order_scaff_NR;
  my @positions_ref;
  my %scaffoldNR_positionRef;
  my %scaffoldNR_numberGenesRef;
  for (my $i = 0; $i < scalar @file; $i++){
    my $geneNR = extract_info($file[$i],"g2");
    my $geneR = extract_info($file[$i],"g1");
    my $currentScaffNR=get_info_gene($geneNR,"c");
    if ($currentScaffNR eq ""){
      next;
    }
    my $coordinateGeneRef = get_info_gene($geneR,"b");
    push @{$scaffoldNR_positionRef{$currentScaffNR}},$coordinateGeneRef;#get all the ref coordinates in an hash{NRscaff}[pos_genes_ortholog_ref]
  }
  my %scaffoldNR_median;
  foreach my $i (keys %scaffoldNR_positionRef){ #calculate the median of orthologue gene positions in the ref scaffold
    my $median = int (scalar @{$scaffoldNR_positionRef{$i}})/2;
    $scaffoldNR_median{$i} = @{$scaffoldNR_positionRef{$i}}[$median];

  }
  @order_scaff_NR = sort { $scaffoldNR_median{$a} <=> $scaffoldNR_median{$b} } keys(%scaffoldNR_median);
  return @order_scaff_NR;
}

#draws the scale line for the reference scaffold
sub draw_guide_lines{
  my $scaffold=$_[0];
  my $length=get_scaffold_length($scaffold);
  $img->moveTo($xminRef-20,$yminRef);
  $img->lineTo($xminRef-20,$ymaxRef);
  my $divisionsScaffold=$length/10;
  my $divisionsLine=($ymaxRef-$yminRef)/10;
  my $currentLength=0;
  my $currenty=$yminRef;
  for(my $i=0;$i<=10;$i++){
    $img->moveTo($xminRef-20,$currenty);
    $img->lineTo($xminRef-25,$currenty);
    $img->moveTo($xminRef-95,$currenty+7);
    my $currentLengthPrint=int($currentLength/1000);
    $img->string("$currentLengthPrint kbp") if $currentLengthPrint == 0;
    $img->string("  $currentLengthPrint") if $currentLengthPrint > 0;
    $currentLength=$currentLength+$divisionsScaffold;
    $currenty=$currenty+$divisionsLine;
  }
}

# draws the scaffold, this subroutine also calls the draw_gene
# sobroutine for representing the genes in the scaffolds
sub draw_scaffold{
  my ($scaffold,$x1,$x2,$y1,$y2,$link)=@_;
  my $finish=get_scaffold_length($scaffold);
  my $step=($y2-$y1)/$finish;#divisions of the scaffold into little pieces which will be used for proportional placement of genes
  if ($link==0){#it's the reference genome
    $stepref=$step;
  }
  for ($i=0; $i < @{$DB_Allscaff_genes{$scaffold}}; ++$i)
  {
    draw_gene($DB_Allscaff_genes{$scaffold}[$i],$step,$x1,$x2,$y1,$y2,$link);
  }
  draw_big_rectangle($scaffold,$x1,$x2,$y1,$y2);
}

# creation of a Database $DB_NRscaffolds_genes{scaffold}[geneNoRef];
sub create_ortholog_scaffold_DB{
  my @orthologGenes=@_;
  my @U=uniq(@orthologGenes);#remove repeated genes
  my @sorted_Genes=sort{ lc($a) cmp lc($b) } @U;
  %DB_NRscaffolds_genes=();
  for (my $i=0;$i<scalar @sorted_Genes;$i++){
    my $scaffold=get_info_gene($sorted_Genes[$i],"c");
    push @{$DB_NRscaffolds_genes{$scaffold}},$sorted_Genes[$i];
  }
}

#creation of database scaffolds{scaff}[gene] including ref and noRef scaffolds
sub create_scaffold_DB_from_connections_file{
  foreach $line (@orthology)
  {
    my @line = split(/\s/,$line);
    my $scaffold1=get_info_gene($line[0],"c");
    my $gene1 = $line[0];
    my $scaffold2=get_info_gene($line[1],"c");
    my $gene2 = $line[1];
    if ($scaffold1 ne ""){
      push @{$DB_Allscaff_genes{$scaffold1}},$gene1;
    }
    if ($scaffold2 ne ""){
      push @{$DB_Allscaff_genes{$scaffold2}},$gene2;
    }
  }
}

sub create_scaffold_DB_from_geneIDlist{
  my $fileName=$_[0];
  open(INFILE,"<",$fileName) || die "# cannot read $fileName";

  while(my $line = <INFILE>){
  	chomp($line);
      my $gene = $line;
      my $scaffold=get_info_gene($line,"c");
      if ($scaffold ne ""){
        push @{$DB_Allscaff_genes{$scaffold}},$gene;
      }
    
  }
}

#creation of a database DBORtholog{geneNoRef}[....genesRef]
sub create_DBOrthologs{
  @orthologGenes=();
  my @file=@_;
  my $prevScaffNR;
  for (my $i=0;$i<scalar@file;$i++){
    my $gene=extract_info($file[$i],"g2");
    my $gene2=extract_info($file[$i],"g1");
    if ($gene ne "-"){
      my $orthologousGene=$gene2;
      if ($gene2 eq "-"){#no connection
        push @{$DB_geneNoRef_orthologGenesRef{$gene}},"-";
      }
      else{
        push @{$DB_geneNoRef_orthologGenesRef{$gene}},$orthologousGene;
      }
      #creation of a vector with all the genes
      push @orthologGenes,$gene;
    }
  }
  return @orthologGenes;
}

# links genes with lines when called by draw_genes,
# each link is created after drawing each gene in
# the noRef scaffold

sub link_orthologs{
  my ($orthologGene,$uppery2,$lowery2)=@_;
  for (my $i=0;$i<scalar@{$DB_geneNoRef_orthologGenesRef{$orthologGene}};$i++){#for each ortholog in the data base for that gene
    my $gene=$DB_geneNoRef_orthologGenesRef{$orthologGene}[$i];
    if ($gene eq "-"){# there is no ortholog connection
      next;
    }
    my $uperyGene=get_info_gene($gene,"b")*$stepref;#get coordinates from reference scaffold
    my $loweryGene=get_info_gene($gene,"e")*$stepref;
    $img->bgcolor($DB_geneRef_colour{$gene}[0]);
    $img->fgcolor($DB_geneRef_colour{$gene}[0]);
    $img->moveTo($xminNoRef,($uppery2+$lowery2)/2);
    $img->lineTo($xmaxRef,($uperyGene+$loweryGene)/2+$yminRef);
  }
}

#opens the file with the name given and stores it in a vector
sub open_file{
  my $fileName=$_[0];
  open (IN,$fileName);
  my @Lecture=<IN>;
  close IN;
    if (@Lecture eq""){
      die "\nfile couldnt be opened\n";
    }
  return @Lecture;
}

#gets a gene and extract the corresponding piece of information, specified by the second parameter
sub extract_info{
  my ($line,$option)=@_;
  $line=~/([^\s]*)\s+([^\s]*)/;
  my $info;
  if ($option eq "g1"){
      my $info=$1;
      return $info;
  }
  elsif ($option eq "g2"){
    my $info=$2;
    return $info;
  }
  return $info;
}

#draws the rectangle representing a scaffold
sub draw_big_rectangle{
  my ($scaffold,$x1,$x2,$y1,$y2)=@_;
  $img->bgcolor(undef);
  $img->fgcolor("black");
  $img->fontsize(16);
  $img->rectangle($x1,$y1,$x2,$y2);
  if ($x1 eq $xminRef){
    $img->moveTo($xminRef,$ymaxRef+50);
  }
  else{
    $img->moveTo($xMaxNoRef+10,($y1+$y2)/2);
  }
  my $print_length=int(get_scaffold_length($scaffold)/1000);
  $img->string("$scaffold\n($print_length"." kbp)");
}

# responsible for drawng genes in the scaffolds
sub draw_gene{
  my ($gene,$scale,$x1,$x2,$y1,$y2,$link)=@_;
  my $geneID=get_info_gene($gene,"n");
  my $beginning=get_info_gene($gene,"b");
  my $end=get_info_gene($gene,"e");
  my $uppery=$y1+$beginning*$scale;
  my $lowery=$y1+$end*$scale;
  my $def_color = $DB_geneRef_colour{$gene}[0];
  if ($def_color="black"){
    $img->bgcolor("black");
    $img->fgcolor("black");
    $img->rectangle($x1,$uppery,$x2,$lowery);
  }
  if ($link == 1){#for non ref scaffolds, colour half genes
    my $numberColors = scalar @{$DB_geneRef_colour{$gene}};
    my $newx1=$x1;
    my $newx2=$x2;
    for (my $i = 0; $i < $numberColors; $i++){
      $newx1=$x1+(($x2-$x1)/$numberColors)*$i;
      $newx2=$x1+(($x2-$x1)/$numberColors)*($i+1);
      $img->bgcolor($DB_geneRef_colour{$gene}[$i]);
      $img->fgcolor($DB_geneRef_colour{$gene}[$i]);
      $img->rectangle($newx1,$uppery,$newx2,$lowery);
    }
  }
  else{
    $img->bgcolor($DB_geneRef_colour{$gene}[0]);
    $img->fgcolor($DB_geneRef_colour{$gene}[0]);
    $img->rectangle($x1,$uppery,$x2,$lowery);
  }
  $img->bgcolor("white");
  $img->fgcolor("black");
  if ($link==1){
    link_orthologs($gene,$uppery,$lowery);
  }
}

# gets gene informatio
sub get_info_gene{
  my ($gene,$option)=@_;
  $gene=~/(.*)__(.*)__(.*)__(.*)/;
  if ($option eq "c"){#chromosome
    my $info=$1;
    return $info;
  }
  elsif ($option eq "n"){#name
    my $info=$4;
    return $info;
  }
  elsif ($option eq "b"){#beginning
    my $info=$2;
    return $info;
  }
  elsif ($option eq "e"){#end
    my $info=$3;
    return $info;
  }
}

#calculates scaffold length by getting the last coordinate fo the last gene in the scaffold
sub get_scaffold_length{
  my $scaffold=$_[0];
  my $defEnd = 0;
  my $end;
  for (my $i=0; $i < scalar @{$DB_Allscaff_genes{$scaffold}};$i++){
    $end=get_info_gene($DB_Allscaff_genes{$scaffold}[$i],"e");
    if ($end > $defEnd){
      $defEnd=$end;
    }
  }
  my $positionLastGene=scalar @{$DB_Allscaff_genes{$scaffold}}-1;
  return $defEnd;
}


#in case there are many scaffolds from the non reference species, the proportions can change
sub get_proportions_orthologs_scaffolds{
  my $totalLen=0;
  foreach my $i(keys %DB_NRscaffolds_genes){
    my $len=get_scaffold_length($i);
    $lenScaff{$i}=$len;
  }
  foreach my $i(keys %lenScaff){
    $totalLen=$totalLen+$lenScaff{$i};
  }
  foreach my $i(keys %lenScaff){
    $lenScaff{$i}=$lenScaff{$i}/$totalLen;
  }
}
