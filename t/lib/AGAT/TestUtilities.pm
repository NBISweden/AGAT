package AGAT::TestUtilities;
use strict;
use warnings;
use Exporter 'import';
use Test::TempDir::Tiny qw(tempdir);
use File::chdir;
use File::Copy qw(copy);
use File::Spec;
use File::Temp qw(tempfile);
use Test::More;

our @EXPORT =
  qw(setup_tempdir check_diff script_prefix check_quiet_run);
my @DIRS;    # keep temp dirs alive until program end

sub setup_tempdir {
    my $dir = tempdir();
    copy('share/agat_config.yaml', File::Spec->catfile($dir, 'agat_config.yaml'));
    copy('share/feature_levels.yaml', File::Spec->catfile($dir, 'feature_levels.yaml'));
    $CWD = $dir;
    push @DIRS, $dir;
    return $dir;
}

sub check_diff {
    my ( $got, $expected, $label, $opts ) = @_;
    $opts //= '';
    # Ignore dynamic footer and usage/help lines
    my $diff_opts = join(' ',
        '-b',
        # any line containing the runtime banner produced by file_text_line
        q{-I '.*Job done in.*'},
        # ignore command line echo and start time
        q{-I '^command :'},
        q{-I '^date :'},
        # ignore pure empty lines that may differ
        q{-I '^$'}
    );
    $diff_opts .= " $opts" if $opts ne '';
    my $diff_output = qx(diff $diff_opts $got $expected 2>&1);
    my $exit_code   = $? >> 8;
    diag("Diff output:\n$diff_output") if $exit_code != 0;
    ok( $exit_code == 0, $label );
}

sub script_prefix {
    # Build a perl prefix for launching repo scripts from tests.
    # - Honors HARNESS_PERL_SWITCHES (e.g., -MDevel::Cover, -I t/lib)
    # - Adds -I <lib> only if a non-system lib path is present in \@INC
    #   (i.e., tests were started with -I lib or use lib 'lib'), and
    #   no -I is already present in HARNESS_PERL_SWITCHES.
    # Local test call must be like: perl -I $PWD/lib t/thetest.t

    my $switches = $ENV{HARNESS_PERL_SWITCHES} // '';
    my $has_cover = ($switches =~ /Devel::Cover/);
    my $has_I_in_switches = ($switches =~ /\s-I\S+/);

    my @lib_to_inject = ();
    if (!$has_I_in_switches) {
        # Mirror a user-provided include dir from \@INC (prefer the first
        # non-system path that contains our AGAT module).
        for my $p (@INC) {
            next unless defined $p && length $p;
            next if $p eq '.';
            # skip common system prefixes
            next if $p =~ m{^/(?:usr|lib|etc)/};
            # if this path contains our local AGAT, use it
            my $agat_pm = File::Spec->catfile($p, 'AGAT');
            if (-e $agat_pm) {
                push @lib_to_inject, $p;
            }
        }
    }
    
    my @parts;
    if ($has_cover || $switches ne '' || @lib_to_inject) {
        push @parts, 'perl';
        push @parts, '-MDevel::Cover' if $has_cover;
        push @parts, $switches if $switches ne '';
        foreach my $lib_to_inject (@lib_to_inject) {
            push @parts, '-I', $lib_to_inject if defined $lib_to_inject;
        }
        return join(' ', @parts) . ' ';
    }
    return '';
}

sub check_quiet_run {
    my ($cmd) = @_;

    my $stdout = File::Spec->catfile( $CWD, 'stdout.txt' );
    my $stderr = File::Spec->catfile( $CWD, 'stderr.txt' );
    my $exit = system("$cmd --verbose 0 1>$stdout 2>$stderr");
    if ( -s $stdout ) {
        open my $out_fh, '<', $stdout or die "Cannot open $stdout: $!";
        local $/;
        my $out = <$out_fh>;
        close $out_fh;
        diag("$stdout:\n$out");
    }
    if ( -s $stderr ) {
        open my $err_fh, '<', $stderr or die "Cannot open $stderr: $!";
        local $/;
        my $err = <$err_fh>;
        close $err_fh;
        diag("$stderr:\n$err");
    }
    ok( -z $stdout, 'stdout is empty' );
    ok( -z $stderr, 'stderr is empty' );
    return $exit;
}

BEGIN {
    setup_tempdir();
}

1;
