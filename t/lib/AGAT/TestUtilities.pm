package AGAT::TestUtilities;
use strict;
use warnings;
use Exporter 'import';
use Test::TempDir::Tiny qw(tempdir);
use File::chdir;
use File::Copy qw(copy);
use File::Spec;
use Test::More;

our @EXPORT = qw(setup_tempdir check_diff script_prefix);
my @DIRS;    # keep temp dirs alive until program end

sub setup_tempdir {
    my $dir = tempdir();
    copy('share/agat_config.yaml', File::Spec->catfile($dir, 'agat_config.yaml'));
    $CWD = $dir;
    push @DIRS, $dir;
    return $dir;
}

sub check_diff {
    my ( $got, $expected, $label, $opts ) = @_;
    $opts //= '';
    my $diff_output = qx(diff $opts $got $expected 2>&1);
    my $exit_code   = $? >> 8;
    diag("Diff output:\n$diff_output") if $exit_code != 0;
    ok( $exit_code == 0, $label );
}

sub script_prefix {
    return (
        $ENV{HARNESS_PERL_SWITCHES}
          && $ENV{HARNESS_PERL_SWITCHES} =~ m/Devel::Cover/
      ) ? 'perl -MDevel::Cover ' : '';
}

BEGIN {
    setup_tempdir();
}

1;
