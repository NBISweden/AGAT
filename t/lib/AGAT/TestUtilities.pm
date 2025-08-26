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
  qw(setup_tempdir check_diff script_prefix check_quiet_run check_quiet_and_normal_run check_console_output);
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
    my $diff_opts = q{-b -I '^Job done in' -I '^usage:'};
    $diff_opts .= " $opts" if $opts ne '';
    my $diff_output = qx(diff $diff_opts $got $expected 2>&1);
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

sub check_quiet_run {
    my ($cmd) = @_;
    my $stdout = File::Spec->catfile( $CWD, 'stdout.txt' );
    my $stderr = File::Spec->catfile( $CWD, 'stderr.txt' );
    my $exit = system("$cmd --quiet 1>$stdout 2>$stderr");
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

sub _sh_escape {
    my ($v) = @_;
    return "''" if !defined $v || $v eq '';
    $v =~ s/'/'"'"'/g;    # POSIX-safe single-quote escaping
    return "'$v'";
}

sub _build_cmd {
    my ( $script, $args_in, $outtmp ) = @_;
    my @parts = ($script);

    # allow a hashref or an arrayref of hashrefs so that the same option
    # (e.g. --gff) can be supplied multiple times
    my @args_sets;
    if ( ref $args_in eq 'ARRAY' ) {
        @args_sets = @$args_in;
    }
    else {
        @args_sets = ($args_in // {});
    }

    for my $args_hr (@args_sets) {
        my %args = %{ $args_hr // {} };
        delete @args{qw/o output/};    # enforce our -o
        for my $k ( sort keys %args ) {
            my $v    = $args{$k};
            my $flag = ( length($k) == 1 ) ? "-$k" : "--$k";    # allow single-letter flags
            if ( !defined $v || $v eq '' || $v eq 1 ) {          # valueless option
                push @parts, $flag;
            }
            else {
                push @parts, $flag, _sh_escape($v);
            }
        }
    }

    push @parts, "-o", _sh_escape($outtmp);
    return join( ' ', @parts );
}

sub check_quiet_and_normal_run {
    my ( $script, $args, $stdout_expected, $results, $out_suffixes ) = @_;
    die "need script"           unless defined $script;
    die "need expected stdout"   unless defined $stdout_expected;

    my @results;
    if ( defined $results ) {
        @results = ref $results eq 'ARRAY' ? @{$results} : ($results);
    }

    my @suffixes;
    if ( defined $out_suffixes ) {
        @suffixes =
          ref $out_suffixes eq 'ARRAY' ? @{$out_suffixes} : ($out_suffixes);
    }
    push @suffixes, (undef) x ( @results - @suffixes );

    # run in quiet mode first
    my $dir    = setup_tempdir();
    my $outtmp = File::Spec->catfile( $dir, 'tmp.gff' );
    my $cmd    = _build_cmd( $script, $args, $outtmp );
    ok( check_quiet_run($cmd) == 0, "quiet run $script" );

    # normal mode
    $dir       = setup_tempdir();
    $outtmp    = File::Spec->catfile( $dir, 'tmp.gff' );
    my $outprefix = File::Spec->catfile( $dir, 'tmp' );
    $cmd       = _build_cmd( $script, $args, $outtmp );
    check_console_output( $cmd, $stdout_expected );

    for my $i ( 0 .. $#results ) {
        my $suffix   = $suffixes[$i];
        my $outfile =
          defined $suffix && $suffix ne ''
          ? $outprefix . $suffix
          : $outtmp;
        check_diff( $outfile, $results[$i], "output $script" );
    }
}

sub check_console_output {
    my ( $cmd, $stdout_fixture ) = @_;
    my $stdout_tmp = File::Spec->catfile( $CWD, 'stdout.txt' );
    my $exit       = system("$cmd --no-progressbar 1>$stdout_tmp 2>&1");
    unless ( -e $stdout_fixture ) {
        diag("$stdout_fixture fixture does not exist");
        return $exit;
    }

    my @ignore_starts   = ( '=> Using standard', 'Using ', 'usage:', 'Reading' , 'Parse file ', 'Parsing ',
                            'Feature discarded by applying the test (see', 'ARG  ' 
                          );
    my @ignore_contains = (
        'AGAT/so.obo',
        '(AGAT) - Version: ',
        'This script is being run',
        'Bioperl location being used:',
        'Operating system being used:',
        ' file parsed',
        'IDs checked and fixed.',
        'Result available in ',
        'done in ',
        '/2025', '/2026', '/2027',      # replace with robust date filtering  
        'Parsing Finished',
        'Compute statistics',
        'Look at the fasta database',
        'load FUNCTIONAL information',
        'Writing result...',
        'End of script.'
    );

    my $filter = sub {
        my ($file) = @_;
        open my $fh, '<', $file or die "Cannot open $file: $!";
        my @lines = <$fh>;
        close $fh;
        return [
            grep {
                my $line = $_;
                !( grep { $line =~ /^\Q$_/ } @ignore_starts )
                  && !( grep { index( $line, $_ ) != -1 } @ignore_contains )
            } @lines
        ];
    };

    my $got      = $filter->($stdout_tmp);
    my $expected = $filter->($stdout_fixture);

    my ( $fh_got, $file_got )             = tempfile();
    my ( $fh_expected, $file_expected )   = tempfile();
    print {$fh_got} @$got;
    print {$fh_expected} @$expected;
    close $fh_got;
    close $fh_expected;

    my $diff_output = qx(diff -b $file_got $file_expected 2>&1);
    my $exit_code   = $? >> 8;
    diag("$stdout_fixture diff:\n$diff_output") if $exit_code != 0;
    ok( $exit_code == 0, "$stdout_fixture matches expected" );

    return $exit;
}

BEGIN {
    setup_tempdir();
}

1;
