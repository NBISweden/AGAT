use 5.006;
use strict;
use warnings FATAL => 'all';
use ExtUtils::MakeMaker;
use File::ShareDir::Install;

install_share dist => 'share';

my %WriteMakefileArgs = (
  NAME            => 'Agat',
  AUTHOR          => 'Jacques Dainat <jacques.dainat@nbis.se>',
  VERSION_FROM    => 'lib/Agat/Omniscient.pm',
  ABSTRACT        => 'Module to deal comprehensively with GFF and GTF format',
  LICENSE         => 'gpl_3',
  PREREQ_PM       => {
                      "Bio::DB::Fasta" => 0,
                      "Bio::DB::Taxonomy" => 0,
                      "Bio::OntologyIO::obo" => 0,
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
                      "File::Share" => 0,
                      "Getopt::Long" => 0,
                      "IO::File" => 0,
                      "IPC::Open2" => 0,
                      "JSON" => 0,
                      "LWP::UserAgent" => 0,
                      "List::MoreUtils" => 0,
                      "POSIX" => 0,
                      "Pod::Usage" => 0,
                      "Scalar::Util" => 0,
                      "Sort::Naturally" => 0,
                      "Statistics::R" => 0,
                      "Time::Piece" => 0,
                      "Time::Seconds" => 0,
                      "Try::Tiny" => 0,
                      "URI::Escape" => 0,
                      "strict" => 0,
                      "warnings" => 0
              },
  TEST_REQUIRES   => {},
  MIN_PERL_VERSION => '5.006',
  BUILD_REQUIRES  => {"ExtUtils::MakeMaker" => 6.31,
                      "File::ShareDir::Install" => 0
                    },
  EXE_FILES => [
          "bin/gff3_sp_Prokka_inferNameFromAttributes.pl",
          "bin/gff3_sp_add_introns.pl",
          "bin/gff3_sp_add_start_and_stop.pl",
          "bin/gff3_sp_alignment_output_style.pl",
          "bin/gff3_sp_clipN_seqExtremities_and_fixCoordinates.pl",
          "bin/gff3_sp_complement_annotations.pl",
          "bin/gff3_sp_ensembl_output_style.pl",
          "bin/gff3_sp_extract_attributes.pl",
          "bin/gff3_sp_extract_sequences.pl",
          "bin/gff3_sp_filter_by_ORF_size.pl",
          "bin/gff3_sp_filter_by_locus_distance.pl",
          "bin/gff3_sp_filter_by_mrnaBlastValue.pl",
          "bin/gff3_sp_filter_incomplete_gene_coding_models.pl",
          "bin/gff3_sp_fix_cds_frame.pl",
          "bin/gff3_sp_fix_features_locations_duplicated.pl",
          "bin/gff3_sp_fix_fusion.pl",
          "bin/gff3_sp_fix_longest_ORF.pl",
          "bin/gff3_sp_fix_overlaping_genes.pl",
          "bin/gff3_sp_fix_small_exon_from_extremities.pl",
          "bin/gff3_sp_flag_short_intron.pl",
          "bin/gff3_sp_functional_statistics.pl",
          "bin/gff3_sp_infer_function_colocalisation.pl",
          "bin/gff3_sp_keep_longest_isoform.pl",
          "bin/gff3_sp_list_short_introns.pl",
          "bin/gff3_sp_load_function_from_protein_align.pl",
          "bin/gff3_sp_manage_IDs.pl",
          "bin/gff3_sp_manage_UTRs.pl",
          "bin/gff3_sp_manage_attributes.pl",
          "bin/gff3_sp_manage_functional_annotation.pl",
          "bin/gff3_sp_manage_introns.pl",
          "bin/gff3_sp_merge_annotations.pl",
          "bin/gff3_sp_splitByLevel2Feature.pl",
          "bin/gff3_sp_statistics.pl",
          "bin/gff3_sp_to_tabulated.pl",
          "bin/gff3_sp_webApollo_compliant.pl",
          "bin/gff3_sq_add_hash_tag.pl",
          "bin/gff3_sq_add_locus_tag.pl",
          "bin/gff3_sq_keep_annotation_from_fastaSeq.pl",
          "bin/gff3_sq_list_attributes.pl",
          "bin/gff3_sq_manage_ID.pl",
          "bin/gff3_sq_mask.pl",
          "bin/gff3_sq_remove_redundant_entries.pl",
          "bin/gff3_sq_repeats_analyzer.pl",
          "bin/gff3_sq_rfam_analyzer.pl",
          "bin/gff3_sq_split.pl",
          "bin/gff3_sq_stat_basic.pl",
          "bin/gff3_sp_gxf_to_gff3.pl"
  ],
  META_MERGE    => {
      'meta-spec' => { version => 2 },
      resources => {
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
  clean => { FILES => 'Agat-*' },
);

my %FallbackPrereqs = (
    "Bio::DB::Fasta" => 0,
    "Bio::DB::Taxonomy" => 0,
    "Bio::OntologyIO::obo" => 0,
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
    "Getopt::Long" => 0,
    "IO::File" => 0,
    "IPC::Open2" => 0,
    "JSON" => 0,
    "LWP::UserAgent" => 0,
    "List::MoreUtils" => 0,
    "POSIX" => 0,
    "Pod::Usage" => 0,
    "Scalar::Util" => 0,
    "Sort::Naturally" => 0,
    "Statistics::R" => 0,
    "Time::Piece" => 0,
    "Time::Seconds" => 0,
    "Try::Tiny" => 0,
    "URI::Escape" => 0,
    "strict" => 0,
    "warnings" => 0
);


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
