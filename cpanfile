requires 'perl', '5.36.0';
requires 'Bio::DB::Fasta';
requires 'Bio::DB::Taxonomy';
requires 'Bio::OntologyIO::obo';
requires 'Bio::Ontology::OntologyEngineI';
requires 'Bio::Seq';
requires 'Bio::SeqIO';
requires 'Bio::Tools::CodonTable';
requires 'Bio::Tools::GFF';
requires 'Carp';
requires 'Clone';
requires 'Cwd';
requires 'DB_File';
requires 'Exporter';
requires 'ExtUtils::InstallPaths';
requires 'File::Basename';
requires 'File::Copy';
requires 'File::Glob';
requires 'File::Share';
requires 'File::ShareDir';
requires 'File::Spec';
requires 'Getopt::Long';
requires 'IO::File';
requires 'IO::Uncompress::Gunzip';
requires 'IPC::Open2';
requires 'YAML';
requires 'LWP::UserAgent';
requires 'LWP::Protocol::https';
requires 'List::MoreUtils';
requires 'Moose';
requires 'POSIX';
requires 'Pod::Usage';
requires 'Scalar::Util';
requires 'Sort::Naturally';
requires 'Term::ProgressBar';
requires 'Time::Piece';
requires 'Time::Seconds';
requires 'Try::Tiny';
requires 'URI::Escape';

recommends 'Parallel::ForkManager';

test_requires 'Capture::Tiny';
test_requires 'File::ShareDir';
test_requires 'Test::Differences';
test_requires 'Test::More', '0.98';
test_requires 'Test::Exception';
test_requires 'Test::TempDir::Tiny';
test_requires 'File::chdir';
test_requires 'Path::Tiny';


on 'develop' => sub {
  requires 'Perl::Critic';
  requires 'Devel::Cover';
  requires 'Devel::Cover::Report::Coveralls';
};
