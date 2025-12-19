#!/usr/bin/env perl
use strict;
use warnings;

my $pod_file = 'AGAT_END_sp.pod';
my $bin_dir  = '../bin';

die "File $pod_file not found\n" unless -f $pod_file;
die "Directory $bin_dir not found\n" unless -d $bin_dir;

# Read AGAT_END.pod content
open my $pod_fh, '<', $pod_file or die "Cannot open $pod_file: $!";
local $/;
my $end_pod = <$pod_fh>;
close $pod_fh;
$end_pod .= "\n" unless $end_pod =~ /\n\z/;

# Process each file in bin/
opendir my $dh, $bin_dir or die "Cannot open $bin_dir: $!";
while (my $file = readdir $dh) {
    next if $file =~ /^\./; # skip . and ..
    # Only scripts starting with agat_sp
    next if $file =~ /^agat_sq/;
    next if $file =~ /^agat_convert_bed/;
    next if $file =~ /^agat_convert_embl/;
    next if $file =~ /^agat_convert_genscan/;
    next if $file =~ /^agat_convert_mfannot/;
    next if $file =~ /^agat_convert_minimap/;
    my $path = "$bin_dir/$file";
    next unless -f $path;

    # Read file content
    open my $fh, '<', $path or die "Cannot open $path: $!";
    local $/;
    my $content = <$fh>;
    close $fh;

    # Only replace if =head1 FEEDBACK exists
    next unless $content =~ /^=head1 SHARED OPTIONS/m;

    # Remove =head1 FEEDBACK and everything after, replace with pod content
    $content =~ s/^=head1 SHARED OPTIONS.*\z/$end_pod/sm;

    # Write back
    open my $out, '>', $path or die "Cannot write $path: $!";
    print $out $content;
    close $out;

    print "Updated: $path\n";
}
closedir $dh;

print "Done.\n";