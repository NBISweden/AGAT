#!/usr/bin/perl -w

package AGAT::PlotR;

use strict;
use warnings;
use Exporter;
use Bio::Tools::GFF;
use Statistics::R;
use AGAT::Utilities;

our @ISA         = qw(Exporter);
our @EXPORT   = qw(execute_R_command rcc_density_one_row_per_file rcc_plot_from_list);
sub import {
  AGAT::PlotR->export_to_level(1, @_); # to be able to load the EXPORT functions when direct call; (normal case)
  AGAT::PlotR->export_to_level(2, @_); # to be able to load the EXPORT functions when called from one level up;
}

=head1 SYNOPSIS

=head1 DESCRIPTION

=head1 AUTHOR

   Jacques Dainat - jacques.dainat@nbis.se

=cut

################## GENERAL #######################
sub execute_R_command {
  my ($command)=@_;

  my $R = Statistics::R->new() or die "Problem with R : $!\n";

#R command
$R->send(
      qq`

      $command

      dev.off()`
          );

# Close the bridge
$R->stopR();
}

################## INPUT FROM LIST #######################

#rcc = R create Command
sub rcc_plot_from_list{
  #@ $list = Ref to a list of numerical value
  #@ $outputPDF = String output name
  #@ $xName = String
  #@ $yaxis = Integer Y axis
  #@ $title = String title
  #@ $plotType = String "Histograme" or other
  my ($list, $legend, $plotType, $xName, $yaxis, $title, $outputPDF)=@_;

  _checkPlotType($plotType);

  my $R_command;
  my $nbColor=1;
  my $opt_breaks = 100;
  # CREATE HEART OF R CODE
  my $list_clean = join ',', @{$list};

  ##################################
  # create main part of R command  #
  ##################################


  # write it for each file
  if($plotType eq "histogram" or $plotType eq "all"){ # write R histogram command
    $R_command.=_histogram_from_list($list_clean, $legend, $xName, $yaxis, $nbColor, $opt_breaks);
  }
  elsif($plotType eq "density" or $plotType eq "all"){ # write R density command
    $R_command.=_density_from_list($list_clean, $legend, $xName, $yaxis, $nbColor, $opt_breaks);
  }


  #add header
    my $final_R_command='#create output
                  pdf("'.$outputPDF.'")
                  # create an empty vector
                  listlegend <- c();';

    #add heart
    $final_R_command.=$R_command;

    #add footer
  $final_R_command .='# Add Title
                title(main="'.$title.'")';
  if($legend){
                #Add Legend
                $final_R_command .= 'legend("topright", col=(1:'.$nbColor.'), lty=1, c(listlegend))';
  }

  return $final_R_command;
}


sub _density_from_list{
  my ($list,$legend,$xName,$yAxisValue,$colorNb,$opt_breaks)=@_;

  my $command='

      listValues<-as.matrix(c('.$list.'))
      legendInfo=paste("'.$legend.'")
      listlegend<-c(listlegend,legendInfo)
      plot(density(listValues),xlim=c(0,1), ylim=c(0,'.$yAxisValue.'), col='.$colorNb.', xlab="'.$xName.'", main="")
      ';

      return $command
}

sub _histogram_from_list{
  my ($list,$legend,$xName,$yAxisValue,$colorNb,$opt_breaks)=@_;

  my $command='

      listValues<-as.matrix(c('.$list.'))
      legendInfo=paste("'.$legend.'")
      listlegend<-c(listlegend,legendInfo)

      a<-seq(min(listValues), max(listValues), '.$opt_breaks.')
      a[length(a)+1]<-max(listValues)

      hist(listValues, breaks=a, col='.$colorNb.', xlab="'.$xName.'", freq = TRUE, main="")
      ';

      return $command
}


################## INPUT FROM FILE #######################

#rcc = R create Command
sub rcc_density_one_row_per_file{
	#@ $fileList = Ref to a list of string tuple file-legend
	#@ $outputPDF = String output name
	#@ $xName = String
	#@ $yaxis = Integer Y axis
	#@ $title = String title
	#@ $plotType = String "Histograme" or other
	my ($fileList, $plotType, $xName, $yaxis, $title, $outputPDF)=@_;

	_checkPlotType($plotType);

	my $R_command;
	my $nbFile=0;

	# CREATE HEART OF R CODE
	foreach my $tuple (@{$fileList}){
		my ($file, $legend)=@$tuple;
   		$nbFile++;
   		my $ref_stream = IO::File->new();
   		$ref_stream->open( $file, 'r' ) or
   		croak(
      		sprintf( "Can not open '%s' for reading: %s", $file, $! ) );

   		##################################
   		# create main part of R command  #
   		##################################
   		if($nbFile > 1){ # only one plot to do
   		  $R_command.='
                 par(new=TRUE)';
   		}

   		# write it for each file
   		if($plotType eq "histogram"){ # write R histogram command
   			$R_command.=_histogram_one_row_in_file($file, $legend, $xName, $yaxis, $nbFile, $outputPDF);
   		}
   		elsif($plotType eq "density"){ # write R density command
   			$R_command.=_density_one_row_in_file($file, $legend, $xName, $yaxis, $nbFile, $outputPDF);
   		}
	}

	#add header
  	my $final_R_command='#create output
                  pdf("'.$outputPDF.'")
                  # create an empty vector
                  listlegend <- c();';

    #add heart
  	$final_R_command.=$R_command;

  	#add footer
 	$final_R_command .='# Add Title
                title(main="'.$title.'")

                #Add Legend
                legend("topright", col=(1:'.$nbFile.'), lty=1, c(listlegend))';

	return $final_R_command;
}


sub _density_one_row_in_file{
  my ($pathIn,$legend,$xName,$yAxisValue,$colorNb,$outputPDF)=@_;

  my $command='

      listValues=as.matrix(read.table("'.$pathIn.'", sep="\t", he=F))
      legendInfo=paste("'.$legend.'")
      listlegend<-c(listlegend,legendInfo)
      plot(density(listValues),xlim=c(0,1), ylim=c(0,'.$yAxisValue.'), col='.$colorNb.', xlab="'.$xName.'", main="")
      ';

      return $command
}

sub _histogram_one_row_in_file{
  my ($pathIn,$legend,$xName,$yAxisValue,$colorNb,$outputPDF)=@_;

  my $command='

      listValues=as.matrix(read.table("'.$pathIn.'", sep="\t", he=F))
      legendInfo=paste("'.$legend.'")
      listlegend<-c(listlegend,legendInfo)

      hist(listValues, breaks=seq(0,100,5), col='.$colorNb.', xlab="'.$xName.'", main="")
      ';

      return $command
}

sub _checkPlotType{
	my ($plotType) = @_;

	my @listPlotOK=("histogram","density");
	my $testOK="false";

	foreach my $type (@listPlotOK){
		if($type eq $plotType){
			$testOK="true";
			last;
		}
	}

	if($testOK eq "false"){
		print "This plot type <$plotType> doesn't exist or has not been yet implemented!\n";exit;
	}


}

1;
