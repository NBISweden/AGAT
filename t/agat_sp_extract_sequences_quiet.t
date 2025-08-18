use strict;
use warnings;
use Test::More tests => 1;
use File::Temp qw(tempdir);
use File::Copy qw(copy);
use Cwd qw(abs_path);

my $tmpdir = tempdir(CLEANUP => 1);
my $repo   = abs_path('.');

copy("$repo/share/agat_config.yaml", "$tmpdir/agat_config.yaml");
system("perl -i -pe 's/^verbose: 1/verbose: 0/' $tmpdir/agat_config.yaml");
system("perl -i -pe 's/^progress_bar: true/progress_bar: false/' $tmpdir/agat_config.yaml");
system("perl -i -pe 's/^log: true/log: false/' $tmpdir/agat_config.yaml");

my $cmd = "$repo/bin/agat_sp_extract_sequences.pl --gff $repo/t/scripts_output/in/1.gff --fasta $repo/t/scripts_output/in/1.fa --config $tmpdir/agat_config.yaml -o $tmpdir/out.fa";
my $output = `$cmd 2>&1`;

is($output, '', 'agat_sp_extract_sequences.pl is quiet with verbose=0');
