use 5.006;
use strict;
use warnings FATAL => 'all';
use ExtUtils::MakeMaker;
use File::ShareDir::Install;
use Cwd qw(abs_path);

install_share dist => 'share';

# ensure temporary test directories are isolated
$ENV{HARNESS_PERL_SWITCHES} = join ' ', '-I' . abs_path('t/lib') . ' -MAGAT::TestUtilities', ($ENV{HARNESS_PERL_SWITCHES} // '');

# ------------------- CREATE LIST OF EXE FILE ----------------------------------

# define function to retrieve list of exe files
sub get_list_exe{
  my @list_exe;

  opendir (DIR, "bin") or die $!;
  while (my $file = readdir(DIR)) {
     # Use a regular expression to ignore files beginning with a period
     next if ($file =~ m/^\./);

     #add exe file
     push @list_exe, "bin/".$file;
  }
  closedir(DIR);

  return \@list_exe;
}
my $ref_list_exe = get_list_exe();


# ------------------- CREATE PREREQ_PM Hash -----------------------------------
# As it is used twice, I factorze it here, to avoid de-synchronization

my %prereq_hash = ( "Bio::DB::Fasta" => 0,
                    "Bio::DB::Taxonomy" => 0,
                    "Bio::OntologyIO::obo" => 0,
                    "Bio::Ontology::OntologyEngineI" => 0,
                    "Bio::Seq" => 0,
                    "Bio::SeqIO" => 0,
                    "Bio::Tools::CodonTable" => 0,
                    "Bio::Tools::GFF" => 0,
                    "Carp" => 0,
                    "Clone" => 0,
                    "Cwd" => 0,
                    "Exporter" => 0,
                    "File::Basename" => 0,
                    "File::Copy" => 0,
		    "File::Glob" => 0,
                    "File::Share" => 0,
	  	    "File::Spec" => 0,
                    "Getopt::Long" => 0,
                    "IO::File" => 0,
                    "IPC::Open2" => 0,
                    "YAML" => 0,
                    "LWP::UserAgent" => 0,
                    "LWP::Protocol::https" => 0,
                    "List::MoreUtils" => 0,
                    "Moose" => 0,
                    "POSIX" => 0,
                    "Pod::Usage" => 0,
                    "Scalar::Util" => 0,
                    "Sort::Naturally" => 0,
                    "Term::ProgressBar" => 0,
                    "Time::Piece" => 0,
                    "Time::Seconds" => 0,
                    "Try::Tiny" => 0,
                    "URI::Escape" => 0,
                    "strict" => 0,
                    "warnings" => 0
                  );

# ------------------- CREATE the WriteMakefile hash ----------------------------
# will be used to create the WriteMakefile object

my %WriteMakefileArgs = (
  NAME            => 'AGAT',
  AUTHOR          => 'Jacques Dainat <jacques.dainat@nbis.se>',
  VERSION_FROM    => 'lib/AGAT/AGAT.pm',
  ABSTRACT        => 'Module to deal comprehensively with GFF and GTF format',
  LICENSE         => 'gpl_3',
  PREREQ_PM       => \%prereq_hash, # give a ref to the hash
  MIN_PERL_VERSION => '5.006',
  BUILD_REQUIRES  => {"ExtUtils::MakeMaker" => 6.31,
                      "File::ShareDir::Install" => 0
                    },
  TEST_REQUIRES   => {
                     "Test::More"        => 0,
                     "Test::TempDir::Tiny" => 0,
                     "File::chdir"       => 0,
                    },
  EXE_FILES       => $ref_list_exe,
  META_MERGE      => {
      'meta-spec' => { version => 2 },
      resources   => {
          repository => {
              type => 'git',
              url  => 'https://github.com/NBISweden/AGAT.git',
              web  => 'https://github.com/NBISweden/AGAT',
          },
          bugtracker => {web => 'https://github.com/NBISweden/AGAT/issues'},
          homepage   => 'https://nbis.se',
      },
  },
  dist  => { COMPRESS => 'gzip -9f', SUFFIX => 'gz', },
  clean => { FILES => 'AGAT-*' },
);

my %FallbackPrereqs = %prereq_hash;


unless ( eval { ExtUtils::MakeMaker->VERSION(6.63_03) } ) {
  delete $WriteMakefileArgs{TEST_REQUIRES};
  delete $WriteMakefileArgs{BUILD_REQUIRES};
  $WriteMakefileArgs{PREREQ_PM} = \%FallbackPrereqs;
}

delete $WriteMakefileArgs{CONFIGURE_REQUIRES}
  unless eval { ExtUtils::MakeMaker->VERSION(6.52) };

if ( $^O eq 'MSWin32' ) {
	$WriteMakefileArgs{PREREQ_PM}{'Win32'} = $FallbackPrereqs{'Win32'} = '0';
}

WriteMakefile(%WriteMakefileArgs);

#Mandatory otherwise the data files from the share folder will not be copied.
package MY;
use File::ShareDir::Install 'postamble';
