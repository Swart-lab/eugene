#!/usr/bin/perl
use strict;

if (scalar(@ARGV) < 2)
{
	print "Error: need two GFF3 files (just CDS or exon lines): first the reference, second the prediction\n";
	exit;
}

my @a_reel = `cat $ARGV[0]` ;
my @a_pred = `cat $ARGV[1]` ;
my $nb_reel = scalar(@a_reel);
my $nb_pred = scalar(@a_pred);
my $cmd = "/usr/local/bioinfo/bin/intersectBed -a ".$ARGV[0]." -b ".$ARGV[1]." -f 1 -r -s -u";
my @a_good_pred = ` $cmd`;
my $nb_good_pred = scalar( @a_good_pred);

my $spec = $nb_good_pred/$nb_pred*100;
my $sens = $nb_good_pred/$nb_reel*100;

print "Specifity:   $spec %. \nSensibility: $sens %.\n";

#my $res = system($cmd);

#print $res;

