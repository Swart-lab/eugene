#!/usr/bin/perl

use strict;
use IO::File;
use Carp;


# >../qd_comp_matrix.pl Output/Output_Tester TestTrace//Output_Tester


my $file1 = $ARGV[0];
my $file2 = $ARGV[1];
my @a_col_to_ignore = (19,20); # doivnet etre ordonne

my $f_1 = new IO::File($file1) or Carp::carp "cant open $file1\n";
my $f_2 = new IO::File($file2) or Carp::carp "cant open $file2\n";
my $nb_line =0;
my $line2;

while (my $line1 = <$f_1>)
{
	$nb_line++;
	$line2 = <$f_2>;
	$line1 =~ s/^ +//;$line2 =~ s/^ +//;
	$line1 =~ s/ +/\t/g;
	$line2 =~ s/ +/\t/g;
	my @a_1 = split(/\t/, $line1);
	my @a_2 = split(/\t/, $line2);

	my $cpt = 0;
	foreach my $field1 (@a_1)
	{
		foreach my $no_col (@a_col_to_ignore)
		{
			if ($cpt eq $no_col)
			{
				$cpt++;
			}
		}
 
		if ($field1 != $a_2[$cpt])
		{
			print "Diff line $nb_line, col $cpt in second file: $field1 != ".$a_2[$cpt]."  \n";
			exit ;
		}
		$cpt++;
	}
}
print "Identical files\n";
