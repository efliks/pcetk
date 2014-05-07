#!/usr/bin/perl -w
#use lib qw(/home/essigke/svn/qmpb/trunk/qmpb);
#
# Use:
# setenv PERL5LIB /home/essigke/svn/qmpb/trunk/qmpb
# setenv MEADPATH /home/essigke/svn/mead/trunk/mead-2.2.3
use strict;
use pqr;
use Data::Dumper;

# Converts the multiflex input files (pqr,mgm,ogm,sites.del,sites.eps,*.st) to a qmpb input.

#############################################################
#
# Default values assumed for qmpb.in:
#
my %global = (
	      T => "300",
	      I => "0.1",
	      backfile => "back.pqr",
	      workdir => "qmpb",
	      meadpath => "$ENV{'MEADPATH'}/bin",
	      epsin1 => "1",
	      epsin2 => "4"
	     );
#
# epsin2 can be overwritten by the epsin command line parameter!
#
# Atoms included in model compound of previous and next residue:
my @prev_res = qw(C O);
#my @next_res = qw(N HN CA); # CHARMM nomenclature, should include HA
#my @next_res = qw(N H CA); # for lysozym example
#my @next_res = qw(N HN CA HA); # full CHARMM 
#my @next_PRO = qw(N CD HD1 HD2 CA HA);
#my @next_GLY = qw(N HN CA HA1 HA2);
my @next_res = qw(N H CA HA); # full CHARMM
my @next_PRO = qw(N CD HD1 HD2 CA HA);
my @next_GLY = qw(N H CA HA1 HA2);
my @aminoacid = qw(ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR HSP HSE HSD HIE HID);
my $conv = - 1.371783713;
#
#
##############################################################
my @ligands;

# Debug output - not used in final version
sub debug {
  my $line = shift;
  print STDERR $line,"\n";
}

# Reading of multiflex style .st-files
sub read_st {
  my $file = shift;
  my $st = shift;
  my %st;
  my $q1 = 0;
  my $q2 = 0;
  open (FH, "$file.st") || die "Can't open $file.st!";
  while (<FH>) {
    # No comments!!!
    my @l = split;
    next if ($#l < 0);
    debug "@l";
    my @t;
    unless (defined $st{Gmodel}) { # pKintr must be in first line
      #      push (@t, -$conv*$l[0]);
      @t = ($conv*$l[0], 0);
      $st{Gmodel} = \@t;
      debug $st{Gmodel}[0];
      next;
    }
    $st{st} = $file;
    $st{center} = $l[1] unless (defined $st{center}); # multiflex takes first atom in .st as center of interest
    $st{resname} = $l[0] unless (defined $st{resname});
    my @t1 = ('1', '0');
    $st{ligand}{proton} = \@t1;
    my @t3 = ('p', 'd');
    $st{label} = \@t3;
    my @t2 = ($l[2], $l[3]);
    $st{atom}{$l[0]}{$l[1]} = \@t2; #st: residue_name atom_name q* = charge
    $q1 += $l[2];
    $q2 += $l[3];
  }
  close (FH);
  die "Column p has a non-integer charge of $q1!" unless  (($q1 < sprintf("%.0f", $q1) + 0.001) && ($q1 > sprintf("%.0f", $q1) - 0.001));
  die "Column d has a non-integer charge of $q2!" unless  (($q2 < sprintf("%.0f", $q2) + 0.001) && ($q2 > sprintf("%.0f", $q2) - 0.001));
  return \%st;
}

# Reading new .est-files
sub read_est {
  my $file = shift; # residue name
  my $st = shift;
  my $no_column = -1;
  my %st;
  $st{center} = undef;
  open (FH, "$file.est") || die "Can't open $file.est!";
  my @q;
  while (<FH>) {
    next if (/^\s+#/); # skip lines with leading comment sign
    next if (/^\s+$/); # skip empty lines
    my @line = split "#";
    $_ = $line[0]; # take only part of line in front of comment sign
    my @tmp = split;
    next if ($#tmp < 1);
    my @l;
    # skip comments at end of line
    foreach my $element (0..$#tmp) {
      last if ($tmp[$element] =~ /#/);
      push (@l, $tmp[$element]);
    }
    debug "@l";
    if ($l[0] =~ /Gmodel/) {
      my @t = @l[1..$#l];
      $st{Gmodel} = \@t; # array slice with all array values
    } elsif ($line[0] =~ /label/) { # column label
      my @t = @l[1..$#l];
      $st{label} = \@t;
    } elsif ($l[0] =~ /center/) {
      $st{center} = $l[1];
    } elsif ($line[0] =~ /^\s*\S+\s+\d+\s+/) { # ligand
      my @t = @l[1..$#l];
      $st{ligand}{$l[0]} = \@t;
    } elsif ($line[0] =~ /^\s*\S+\s+\S+\s+[\d.-]+/) { # atom line
      my @t = @l[2..$#l];

      # Some extra code for error checking
      if ($#q < 0) {
	@q = @t;
      } elsif ($#q == $#t) {
	$q[$_] += $t[$_] foreach (0..$#t);
      } else {
	die "$file.est: Number of columns not equal!\n";
      }

      $st{atom}{$l[0]}{$l[1]} = \@t;
      $st{resname} = $l[0] unless (defined $st{resname});
    } else {
      die "$file.est: Unrecognized line $line[0]!\n";
    }
  }
  close (FH);
#  print Dumper %st;
  # Error checking...
  die "Gmodel not defined in $file.est!\n" unless (defined  $st{Gmodel});
  $no_column = $#{$st{Gmodel}}+1;

  debug ("$no_column instances for residue $file");

  die "Number of labels in residue $file expected to be $no_column, but got ", $#{$st{label}}+1, "!\n" unless ($#{$st{label}}+1 == $no_column);
  die "No ligand defined in $file.est!\n" unless (defined  $st{ligand});
  foreach (keys %{$st{ligand}}) {
    die "Number of ligands $_ in residue $file expected to be $no_column, but got ", $#{$st{ligand}{$_}}, "!\n" unless ($#{$st{ligand}{$_}}+1 == $no_column);
  }
  die "No atoms given in $file.est!\n" unless (defined $st{atom});
  my @res = keys %{$st{atom}};
  die "Only one residue can be given in $file.est!\n" unless ($#res == 0);
  foreach (keys %{$st{atom}{$res[0]}}) {
    die "Number of charges for atom $_ in residue $res[0] expected to be $no_column, but got ", $#{$st{atom}{$res[0]}{$_}}+1 ,"!\n" unless ($#{$st{atom}{$res[0]}{$_}}+1 == $no_column);
  }
  die "Center atom $st{atom} given in $file.st not part of site!\n" if ((defined $st{center}) && (!defined $st{atom}{$file}{$st{center}}));
  foreach my $column (0..$#{$st{label}}) {
    my $round = sprintf("%.1f", $q[$column]);
#    my $tmp = $q[$column];
#    debug "$tmp, round: $round ";
    print (STDERR "Check column $column for integer charge: ", $q[$column], " : $round\n");
    die "Column $st{label}[$column] has a non-integer charge of $q[$column]!" unless  (($q[$column] < sprintf("%.0f", $q[$column]) + 0.001) && ($q[$column] > sprintf("%.0f", $q[$column]) - 0.001));
  }
#  debug "\n";
  return \%st;
}



# Reading the multiflex .sites-file
# If .est-files are available they are read, else the .st-file is read
sub read_sites {
  my $file = shift;
  my $st = shift;
  my %site;
  open (FH1, $file) || die "Can't open $file!";
  my $sid = 0;
  while (<FH1>) {
    next if (/^\s+$/);
    my @l = split;
    debug "@l";
    unless (defined $$st{$l[1]}) {
      if (-e "$l[1].est") {
	$$st{$l[1]} = read_est ($l[1]);
      } elsif (-e "$l[1].st") {
	$$st{$l[1]} = read_st ($l[1]);
      } else {
	die "$l[1].est or $l[1].st not found!";
      }
    }
    my $chain = 0;
    $chain = $l[2] if ($#l > 1);
    $site{st}{$chain}{$l[0]} = $$st{$l[1]}; #site: chain, residue_number, pointer_to_st
    my $site_label = $$st{$l[1]}{resname}.$l[0]."-".$chain;
    $site{order_sid}{$site_label} = $sid;
    $sid++;
    foreach my $iid (0..$#{$$st{$l[1]}{label}}) {
      $site{order_iid}{$site_label}{$$st{$l[1]}{label}[$iid]} = $iid;
    }
  }
  close (FH1);
  return \%site;
}

# Read .ogm or .mgm file
sub read_grid {
  my $file = shift;
  my %grid;
  open (FH, $file) || die "Can't open $file!";
  while (<FH>) {
    my @l = split;
    push(@{$grid{center}}, $l[0]);
    push(@{$grid{point}}, $l[1]);
    push(@{$grid{space}}, $l[2]);
  }
  close (FH);
  return \%grid;
}

# Read .pqr-file
sub read_pqr {
  my $file = shift;
  return pqr->read($file);
}

sub write_back_pqr {
  my $pqr = shift;
  my $site = shift;
  my $file = shift;
  my $mypqr = pqr->clone($pqr);
  foreach my $chain (keys %$site) {
    foreach my $resno (keys %{$$site{$chain}}) {
      foreach my $atname (keys %{$$site{$chain}{$resno}{atom}{$$site{$chain}{$resno}{resname}}}) {
	debug "$resno, $atname";
#      next if ($atname eq "pKmodel" || $atname eq "center" ||$atname eq "resname" );
	$mypqr->delete_atom($resno, $atname, 0, $chain) || die "Atom $atname not found in residue $resno!\n";
      }
    }
  }
  $mypqr->write($file, 0, 0, 0, 0);
}

sub write_site {
#  print "write site: @_\n";
  my $name = shift;
  my $site = shift;
  my $res_no = shift;
  my $pqr = shift;
  my $inst = shift; # instance number
  my $chain = shift;
  my $conf = 0;
  mkdir $chain unless (-d $chain);
  open(FH, ">$chain/$name") || die "Can't open file $chain/$name for writing!\n";
#  die "DEB: not defined ($chain, $res_no): $site\n";
  my $resname = $$site{$chain}{$res_no}{resname};
  foreach my $atom (keys %{$$site{$chain}{$res_no}{atom}{$resname}}) {
    print FH $pqr->print_atom($res_no, $atom, $$site{$chain}{$res_no}{atom}{$resname}{$atom}[$inst], 0, 0, $conf, $chain);
  }
  close (FH);
}

sub write_model {
  my $name = shift;
  my $site = shift;
  my $res_no = shift;
  my $pqr = shift;
  my $proto = shift;
  my $chain = shift;
  my $resname = $$site{$chain}{$res_no}{resname};
  my $isAA = 0;
  foreach my $AA (@aminoacid) {
    if ($AA eq $resname) {
      $isAA = 1;
      last;
    }
  }
  mkdir $chain unless (-d $chain);
  open(FH, ">$chain/$name") || die "Can't open file $chain/$name for writing!\n";
  if ($isAA) {
    foreach (@prev_res) {
      print FH $pqr->print_atom($res_no-1, $_, undef, 0, 0, 0, $chain) if ($pqr->atom_exists($res_no-1, $_, 0, $chain, 0, 0));
    }
  }
  foreach my $atom_name (@{$pqr->get_atoms_of_residue($res_no, 0, $chain, 0, 0)}) { #atom name
#    debug "$_";
#    print "Model 1: ", Dumper $site;
    unless (exists ($$site{$chain}{$res_no}{atom}{$resname}{$atom_name})) {
      print FH $pqr->print_atom($res_no, $atom_name, undef, 0, 0, 0, $chain);
#      print "Model 2: ", Dumper $site;
    } else {
      print FH $pqr->print_atom($res_no, $atom_name, $$site{$chain}{$res_no}{atom}{$resname}{$atom_name}[$proto], 0, 0, 0, $chain);
      #      print "Model 3: ", Dumper $site;
    }
  }
  if ($isAA) {
    my $next_resname = $pqr->get_residue_name($res_no+1, 0, $chain, 0, 0);
    if ((defined $next_resname) && ($next_resname eq "GLY")) {
      foreach (@next_GLY) {
	print FH $pqr->print_atom($res_no+1, $_, undef, 0, 0, 0, $chain) if ($pqr->atom_exists($res_no+1, $_, 0, $chain, 0, 0));
      }
    } elsif ((defined $next_resname) && ($next_resname eq "PRO")) {
      foreach (@next_PRO) {
	print FH $pqr->print_atom($res_no+1, $_, undef, 0, 0, 0, $chain) if ($pqr->atom_exists($res_no+1, $_, 0, $chain, 0, 0));
      }
    } else {
      foreach (@next_res) {
	print FH $pqr->print_atom($res_no+1, $_, undef, 0, 0, 0, $chain) if ($pqr->atom_exists($res_no+1, $_, 0, $chain, 0, 0));
      }
    }
  }
  close (FH);
}

sub print_site {
  my $fullsite = shift;
  my $site = $$fullsite{st};
  my $res_no = shift;
  my $chain = shift;
  my $pqr = shift;
#  print Dumper $site;
  print "\n#####################\n\n";
  foreach my $inst (0..$#{$$site{$chain}{$res_no}{Gmodel}}) {
    my $inst_label;
    my $model_label;
    if (defined $$site{$chain}{$res_no}{label}[$inst]) {
      $inst_label = $$site{$chain}{$res_no}{label}[$inst];
      $model_label = "model_".$$site{$chain}{$res_no}{label}[$inst];
    } else { # should never be used!
      $inst_label = "inst_{$inst}";
      $model_label = "model_{$inst}";
    }
    my $site_label = $$site{$chain}{$res_no}{resname}.$res_no."-".$chain;
    print "MMsite ", $site_label, " $inst_label\n";
    my $name = "site_".$$site{$chain}{$res_no}{resname}."_".$res_no."_".$inst_label.".pqr";
#    print "Site: $site\n";
    write_site ($name, $site, $res_no, $pqr, $inst, $chain);
    print "  file=$chain/$name\n  ref=$model_label\n  Gmm=0\n  sid=".$$fullsite{order_sid}{$site_label}."\n  iid=".$$fullsite{order_iid}{$site_label}{$inst_label}."\n";
    if (defined $$site{$chain}{$res_no}{center}) {
      my @c = $pqr->get_atom_coor($res_no, $$site{$chain}{$res_no}{center}, 0, $chain, 0, 0);
      my $center = "$c[0], $c[1], $c[2]";
      print "  center=$center\n  nohomo\n\n";
    } else {
      print "  nohomo\n\n";
    }
    print "Modelsite ", $site_label, "  $model_label\n";
    $name = "model_".$$site{$chain}{$res_no}{resname}."_".$res_no."_".$inst_label.".pqr";
    write_model ($name, $site, $res_no, $pqr, $inst, $chain);
    my @N;
    foreach my $lig_name (@ligands) {
      if (defined $$site{$chain}{$res_no}{ligand}{$lig_name}) {
	push (@N, $$site{$chain}{$res_no}{ligand}{$lig_name}[$inst]);
      } else {
	push (@N, 0);
      }
    }
    print "  file=$chain/$name\n  ref=$inst_label\n  Gmodel= ",  $$site{$chain}{$res_no}{Gmodel}[$inst], "\n  sid=".$$fullsite{order_sid}{$site_label}."\n  iid=".$$fullsite{order_iid}{$site_label}{$inst_label}."\n  eps=$global{epsin2}\n  N=@N\n\n";
  }
}


my $file = $ARGV[0];
die "Usage: $0 molname > qmpb.in\n" if ($#ARGV != 0);
my %st;
my $fullsite = read_sites ("$file.sites", \%st);
my $site = $$fullsite{st};
my %seen = ();
# Get unique list of ligand names for all sites
foreach my $chain (keys %$site) {
  foreach my $resno (keys %{$$site{$chain}}) {
    foreach my $ligand (keys %{$$site{$chain}{$resno}{ligand}}) {
      push (@ligands, $ligand) unless $seen{$ligand}++;
    }
  }
}
my $pqr = read_pqr ("$file.pqr");
my $ogm = read_grid ("$file.ogm");
my $mgm = read_grid ("$file.mgm");
write_back_pqr ($pqr, $site, $global{backfile});

print "# Simulate Mutiflex by QMPB\n";
print "# Input generated by $0\n\n";
print "$_ = $global{$_}\n" foreach (sort keys %global);
print "MGMcenter = @{$$mgm{center}}\n";
print "MGMpoints = @{$$mgm{point}}\n";
print "MGMspace = @{$$mgm{space}}\n";
print "OGMcenter = @{$$ogm{center}}\n";
print "OGMpoints = @{$$ogm{point}}\n";
print "OGMspace = @{$$ogm{space}}\n";
print "Ligand_Labels = @ligands\n";
print "\n\n";
foreach my $chain (sort keys %$site) {
  foreach my $resno (sort keys %{$$site{$chain}}) {
    print_site ($fullsite, $resno, $chain, $pqr);
  }
}
