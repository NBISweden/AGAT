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
    my ( $script, $args_hr, $outtmp ) = @_;
    my @parts = ($script);
    my %args  = %{ $args_hr // {} };

    delete @args{qw/o output/};    # enforce our -o
    for my $k ( sort keys %args ) {
        my $v    = $args{$k};
        my $flag = ( length($k) == 1 ) ? "-$k" : "--$k";
        if ( !defined $v || $v eq '' || $v eq 1 ) {
            push @parts, $flag;
        }
        else {
            push @parts, $flag, _sh_escape($v);
        }
    }
    push @parts, "-o", _sh_escape($outtmp);
    return join( ' ', @parts );
}

sub check_quiet_and_normal_run {
    my ( $script, $args_hr, $result, $out_suffix ) = @_;
    die "need script" unless defined $script;
    die "need result" unless defined $result;

    for my $mode (qw/quiet console/) {
        my $dir       = setup_tempdir();
        my $outtmp    = File::Spec->catfile( $dir, 'tmp.gff' );
        my $outprefix = File::Spec->catfile( $dir, 'tmp' );
        my $cmd       = _build_cmd( $script, $args_hr, $outtmp );

        if ( $mode eq 'quiet' ) {
            check_quiet_run( $cmd, $result );
        }
        else {
            check_console_output( $cmd, $result );
        }

        if ($out_suffix) {
            check_diff( $outprefix . $out_suffix, $result,
                "$mode output $script" );
        }
        else {
            check_diff( $outtmp, $result, "$mode output $script" );
        }
    }
}

sub check_console_output {
    my ( $cmd, $result ) = @_;
    my $stdout_fixture = $result . '.stdout';
    my $stdout_tmp     = File::Spec->catfile( $CWD, 'stdout.txt' );
    my $exit           = system("$cmd --no-progressbar 1>$stdout_tmp 2>&1");
    unless ( -e $stdout_fixture ) {
        diag("$stdout_fixture fixture does not exist");
        return $exit;
    }

    my @ignore_starts   = ( '=> Using standard', 'Using standard', 'usage:', 'Reading' );
    my @ignore_contains = (
        'AGAT/so.obo',
        'This script is being run',
        'Bioperl location being used:',
        'Operating system being used:',
        'done in '
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

    my $diff_output = qx(diff $file_got $file_expected 2>&1);
    my $exit_code   = $? >> 8;
    diag("$stdout_fixture diff:\n$diff_output") if $exit_code != 0;
    ok( $exit_code == 0, "$stdout_fixture matches expected" );

    return $exit;
}

BEGIN {
    setup_tempdir();
}

1;
