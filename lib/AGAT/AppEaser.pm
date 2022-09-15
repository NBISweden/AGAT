#!/usr/bin/perl -w

package AGAT::AppEaser;
use v5.24;
use warnings;
use experimental qw< signatures >;
no warnings qw< experimental::signatures >;
{ our $VERSION = '0.011' }

use Exporter 'import';
our @EXPORT_OK = qw< d run >;

sub add_auto_commands ($application) {
   my $commands = $application->{commands};
   $commands->{help} //= {
      name                     => 'help',
      supports                 => ['help'],
      help                     => 'print a help message',
      description              => 'print help for (sub)command',
      'allow-residual-options' => 0,
      leaf                     => 1,
      execute                  => \&stock_help,
   };
   $commands->{commands} //= {
      name                     => 'commands',
      supports                 => ['commands'],
      help                     => 'list sub-commands',
      description              => 'Print list of supported sub-commands',
      'allow-residual-options' => 0,
      leaf                     => 1,
      execute                  => \&stock_commands,
   };
   return $application;
} ## end sub add_auto_commands ($application)

sub collect ($self, $spec, $args) {
   my @sequence;
   my $config = {};
   my @residual_args;

   my $merger = merger($self, $spec);

   for my $source_spec (sources($self, $spec, $args)) {
      my ($src, $src_cnf) =
        'ARRAY' eq ref $source_spec
        ? $source_spec->@*
        : ($source_spec, {});
      $src = $self->{factory}->($src, 'collect');    # "resolve"
      $src_cnf = {$spec->%*, $src_cnf->%*, config => $config};
      my ($slice, $residual_args) = $src->($self, $src_cnf, $args);
      push @residual_args, $residual_args->@* if defined $residual_args;
      push @sequence, $slice;
      $config = $merger->(@sequence);
   } ## end for my $source_spec (sources...)

   return ($config, \@residual_args);
} ## end sub collect

sub collect_options ($self, $spec, $args) {
   my $factory = $self->{factory};
   my $collect = $spec->{collect}
     // $self->{application}{configuration}{collect} // \&collect;
   my $collector = $factory->($collect, 'collect');    # "resolve"
   (my $config, $args) = $collector->($self, $spec, $args);
   push $self->{configs}->@*, $config;
   return $args;
} ## end sub collect_options

sub commandline_help ($getopt, $shortbool) {
   my @retval;

   my ($mode, $type, $desttype, $min, $max, $default);
   if (substr($getopt, -1, 1) eq '!') {
      $type = 'bool';
      substr $getopt, -1, 1, '';
      push @retval, 'boolean option';
   }
	 elsif (substr($getopt, -2, 1) eq '!-') {
			$type = 'boolshort';
			substr $getopt, -2, 1, '';
			push @retval, 'boolean option';
	 }
   elsif (substr($getopt, -1, 1) eq '+') {
      $mode = 'increment';
      substr $getopt, -1, 1, '';
      push @retval,
        'incremental option (adds 1 every time it is provided)';
   } ## end elsif (substr($getopt, -1...))
   elsif (
      $getopt =~ s<(
         [:=])    # 1 mode
         ([siof]) # 2 type
         ([@%])?  # 3 desttype
         (?:
            \{
               (\d*)? # 4 min
               ,?
               (\d*)? # 5 max
            \}
         )? \z><>mxs
     )
   {
      $mode     = $1 eq '=' ? 'mandatory' : 'optional';
      $type     = $2;
      $desttype = $3;
      $min      = $4;
      $max      = $5;
      if (defined $min) {
         $mode = $min ? 'optional' : 'required';
      }
      $type = {
         s => 'string',
         i => 'integer',
         o => 'perl-extended-integer',
         f => 'float',
      }->{$type};
      my $line = "$mode $type option";
      $line .= ", at least $min times" if defined($min) && $min > 1;
      $line .= ", no more than $max times"
        if defined($max) && length($max);
      $line .= ", list valued" if defined($desttype) && $desttype eq '@';
      push @retval, $line;
   } ## end elsif ($getopt =~ s<(          )? \z><>mxs)
   elsif ($getopt =~ s<: (\d+) ([@%])? \z><>mxs) {
      $mode     = 'optional';
      $type     = 'i';
      $default  = $1;
      $desttype = $2;
      my $line = "optional integer, defaults to $default";
      $line .= ", list valued" if defined($desttype) && $desttype eq '@';
      push @retval, $line;
   } ## end elsif ($getopt =~ s<: (\d+) ([@%])? \z><>mxs)
   elsif ($getopt =~ s<:+ ([@%])? \z><>mxs) {
      $mode     = 'optional';
      $type     = 'i';
      $default  = 'increment';
      $desttype = $1;
      my $line = "optional integer, current value incremented if omitted";
      $line .= ", list valued" if defined($desttype) && $desttype eq '@';
      push @retval, $line;
   } ## end elsif ($getopt =~ s<:+ ([@%])? \z><>mxs)

   my @alternatives = split /\|/, $getopt;
   if ($type eq 'bool') {
		 if ($shortbool){
			 push @retval, map {
					if   (length($_) == 1) { "-$_" }
					else                   { "--$_" }
			 } @alternatives;
		 }
		 else {
      push @retval, map {
         if   (length($_) == 1) { "-$_" }
         else                   { "--$_ | --no-$_" }
      } @alternatives;
		 }
   } ## end if ($type eq 'bool')
   elsif ($mode eq 'optional') {
      push @retval, map {
         if   (length($_) == 1) { "-$_ [<value>]" }
         else                   { "--$_ [<value>]" }
      } @alternatives;
   } ## end elsif ($mode eq 'optional')
   else {
      push @retval, map {
         if   (length($_) == 1) { "-$_ <value>" }
         else                   { "--$_ <value>" }
      } @alternatives;
   } ## end else [ if ($type eq 'bool') ]

   return @retval;
} ## end sub commandline_help ($getopt)

sub commit_configuration ($self, $spec, $args) {
   my $commit = $spec->{commit} // return;
   $self->{factory}->($commit, 'commit')->($self, $spec, $args);
}

sub d (@stuff) {
   no warnings;
   require Data::Dumper;
   local $Data::Dumper::Indent = 1;
   warn Data::Dumper::Dumper(@stuff % 2 ? \@stuff : {@stuff});
} ## end sub d (@stuff)

sub default_getopt_config ($self, $spec) {
   my @r = qw< gnu_getopt >;
   push @r, qw< require_order pass_through >
      if has_children($self, $spec);
   push @r, qw< pass_through > if  $spec->{'allow-residual-options'};
   return \@r;
}

sub execute ($self, $args) {
   my $command    = $self->{trail}[-1][0];
   my $executable = fetch_spec_for($self, $command)->{execute}
     or die "no executable for '$command'\n";
   $executable = $self->{factory}->($executable, 'execute');    # "resolve"
   my $config = $self->{configs}[-1] // {};
   return $executable->($self, $config, $args);
} ## end sub execute

sub fetch_subcommand_default ($self, $spec) {
   my $acfg = $self->{application}{configuration};
   my $child = exists($spec->{'default-child'}) ? $spec->{'default-child'}
      : exists($acfg->{'default-child'}) ? $acfg->{'default-child'}
      : get_child($self, $spec, 'help'); # help is last resort
   return ($child, $child) if defined $child && length $child;
   return;
}

sub fetch_subcommand ($self, $spec, $args) {
   my ($subc, $alias) = fetch_subcommand_wh($self, $spec, $args)
      or return;
   my $r = ref $subc;
   if ($r eq 'HASH') {
      $subc = $spec->{children}[$subc->{index}]
         if scalar(keys $subc->%*) == 1 && defined $subc->{index};
      $r = ref $subc;
      return ($subc, $subc->{supports}[0]) if $r eq 'HASH';
      $alias = $subc;
   }
   die "invalid sub-command (ref to $r)" if $r;
   return ($subc, $alias);
}

sub fetch_subcommand_wh ($self, $spec, $args) {
   # if there's a dispatch, use that to figure out where to go next
   # **this** might even overcome having children at all!
   for my $cfg ($spec, $self->{application}{configuration}) {
      next unless exists $cfg->{dispatch};
      my $sub = $self->{factory}->($cfg->{dispatch}, 'dispatch');
      defined(my $child = $sub->($self, $spec, $args)) or return;
      return ($child, $child);
   }

   # regular course here, no point in going forth without children
   return unless has_children($self, $spec);

   # use defaults if there's no argument to investigate
   return fetch_subcommand_default($self, $spec) unless $args->@*;

   # try to get a child from the first argument
   if (my $child = get_child($self, $spec, $args->[0])) {
      return ($child, shift $args->@*); # consumed arg name
   }

   # the first argument didn't help, but we might want to fallback
   for my $cfg ($spec, $self->{application}{configuration}) {
      if (exists $cfg->{fallback}) { # executable
         defined(my $fb = $cfg->{fallback}) or return;
         my $sub = $self->{factory}->($fb, 'fallback'); # "resolve"
         defined(my $child = $sub->($self, $spec, $args)) or return;
         return ($child, $child);
      }
      if (exists $spec->{'fallback-to'}) {
         defined(my $fbto = $spec->{'fallback-to'}) or return;
         return ($fbto, $fbto);
      }
      return fetch_subcommand_default($self, $spec)
         if $cfg->{'fallback-to-default'};
   }

   # no fallback at this point... it's an error, build a message and die!
   my @names = map { $_->[1] } $self->{trail}->@*;
   shift @names;    # remove first one
   my $path = join '/', @names, $args->[0]; # $args->[0] was the candidate
   die "cannot find sub-command '$path'\n";
} ## end sub fetch_subcommand_wh

sub generate_factory ($c) {
   my $w = \&stock_factory;    # default factory
   $w = stock_factory($c->{create}, 'factory', $c) if defined $c->{create};
   return sub ($e, $d = '') { $w->($e, $d, $c) };
}

sub get_child ($self, $spec, $name) {
   for my $child (get_children($self, $spec)) {
      my $command = fetch_spec_for($self, $child);
      next
        unless grep { $_ eq $name }
        ($command->{supports} //= [$child])->@*;
      return $child;
   } ## end for my $child (get_children...)
   return;
} ## end sub get_child

sub get_children ($self, $spec) {
   return if $spec->{leaf};
   return if exists($spec->{children}) && !$spec->{children};
   my @children = ($spec->{children} // [])->@*;

   # set auto-leaves as 1 by default, new in 0.007002
   $self->{application}{configuration}{'auto-leaves'} = 1
      unless exists $self->{application}{configuration}{'auto-leaves'};

   return
     if $self->{application}{configuration}{'auto-leaves'}
     && @children == 0;    # no auto-children for leaves under auto-leaves
   my @auto =
     exists $self->{application}{configuration}{'auto-children'}
     ? (($self->{application}{configuration}{'auto-children'} // [])->@*)
     : (qw< help commands >);
   if (exists $spec->{'no-auto'}) {
      if (ref $spec->{'no-auto'}) {
         my %no = map { $_ => 1 } $spec->{'no-auto'}->@*;
         @auto = grep { !$no{$_} } @auto;
      }
      elsif ($spec->{'no-auto'} eq '*') {
         @auto = ();
      }
      else {
         die "invalid no-auto, array or '*' are allowed\n";
      }
   } ## end if (exists $spec->{'no-auto'...})
   return (@children, @auto);
} ## end sub get_children

# traverse a whole @$list of sub-commands from $start. This is used to
# list "commands" at a certain sub-level or show help
sub get_descendant ($self, $start, $list) {
   my $target = $start;
   my $path;
   for my $desc ($list->@*) {
      $path = defined($path) ? "$path/$desc" : $desc;
      my $command = fetch_spec_for($self, $target)
        or die "cannot find sub-command '$path'\n";
      defined($target = get_child($self, $command, $desc))
        or die "cannot find sub-command '$path'\n";
   } ## end for my $desc ($list->@*)

   # check that this last is associated to a real command
   return $target if fetch_spec_for($self, $target);
   die "cannot find sub-command '$path'\n";
} ## end sub get_descendant

sub has_children ($self, $spec) { get_children($self, $spec) ? 1 : 0 }

sub hash_merge {
   my (%retval, %is_overridable);
   for my $href (@_) {
      for my $src_key (keys $href->%*) {
         my $dst_key = $src_key;
         my $this_overridable;
         if ($dst_key =~ m{\A //= (.*) \z}mxs) { # overridable
            $dst_key = $1;
            $is_overridable{$dst_key} = 1 unless exists $retval{$dst_key};
            $this_overridable = 1;
         }
         $retval{$dst_key} = $href->{$src_key}
            if $is_overridable{$dst_key} || ! exists($retval{$dst_key});
         $is_overridable{$dst_key} = 0 unless $this_overridable;
      }
   }
   return \%retval;
   # was a simple: return {map { $_->%* } reverse @_};
}

sub list_commands ($self, $children) {
   my $retval = '';
   open my $fh, '>', \$retval;
   for my $child ($children->@*) {
      my $command = fetch_spec_for($self, $child);
      my $help    = $command->{help};
      my @aliases = ($command->{supports} // [$child])->@*;
      next unless @aliases;
      printf {$fh} "%15s: %s\n", shift(@aliases), $help;
      printf {$fh} "%15s  (also as: %s)\n", '', join ', ', @aliases
        if @aliases;
   } ## end for my $child ($children...)
   close $fh;
   return $retval;
} ## end sub list_commands

sub load_application ($application) {
   return $application if 'HASH' eq ref $application;

   my $text;
   if ('SCALAR' eq ref $application) {
      $text = $$application;
   }
   else {
      my $fh =
        'GLOB' eq ref $application
        ? $application
        : do {
         open my $fh, '<:encoding(UTF-8)', $application
           or die "cannot open '$application'\n";
         $fh;
        };
      local $/;    # slurp mode
      $text = <$fh>;
      close $fh;
   } ## end else [ if ('SCALAR' eq ref $application)]

   return eval {
      require JSON::PP;
      JSON::PP::decode_json($text);
   } // eval { eval $text; } // die "cannot load application\n";
} ## end sub load_application ($application)

sub merger ($self, $spec = {}) {
   my $merger = $spec->{merge}
     // $self->{application}{configuration}{merge} // \&hash_merge;
   return $self->{factory}->($merger, 'merge');    # "resolve"
}

sub env_namer ($self, $cspec) {
   my $namenv = $cspec->{namenv}
     // $self->{application}{configuration}{namenv} // \&stock_NamEnv;
   $namenv = $self->{factory}->($namenv, 'namenv'); # "resolve"
   return sub ($ospec) { $namenv->($self, $cspec, $ospec) };
} ## end sub name_for_option ($o)

sub name_for_option ($o) {
   return $o->{name} if defined $o->{name};
   return $1 if defined $o->{getopt} && $o->{getopt} =~ m{\A(\w+)}mxs;
   return lc $o->{environment}
      if defined $o->{environment} && $o->{environment} ne '1';
   return '~~~';
} ## end sub name_for_option ($o)

sub params_validate ($self, $spec, $args) {
   my $validator = $spec->{validate}
     // $self->{application}{configuration}{validate} // return;
   require Params::Validate;
   Params::Validate::validate($self->{configs}[-1]->%*, $validator);
} ## end sub params_validate

sub print_commands ($self, $target) {
   my $command = fetch_spec_for($self, $target);
   my $fh =
     $self->{application}{configuration}{'help-on-stderr'}
     ? \*STDERR
     : \*STDOUT;
   if (my @children = get_children($self, $command)) {
      print {$fh} list_commands($self, \@children);
   }
   else {
      print {$fh} "no sub-commands\n";
   }
}

sub print_help ($self, $target) {
   my $command = fetch_spec_for($self, $target);
   my $enamr   = env_namer($self, $command);
   my $fh =
     $self->{application}{configuration}{'help-on-stderr'}
     ? \*STDERR
     : \*STDOUT;

   print {$fh} $command->{help}, "\n\n";

   if (defined(my $description = $command->{description})) {
      $description =~ s{\A\s+|\s+\z}{}gmxs;    # trim
      $description =~ s{^}{    }gmxs;          # add some indentation
      print {$fh} "Description:\n$description\n\n";
   }

   printf {$fh} "Can be called as: %s\n\n", join ', ',
     $command->{supports}->@*
     if $command->{supports};

   my $options = $command->{options} // [];
   if ($options->@*) {
      print {$fh} "Options:\n";
      my $n = 0; # count the option
      for my $option ($options->@*) {
         print {$fh} "\n" if $n++;

         printf {$fh} "%15s: %s\n", name_for_option($option),
           $option->{help} // '';

         if (exists $option->{getopt}) {
            my $shortbool;  if (exists $option->{shortbool}){$shortbool=1;};
            my @lines = commandline_help($option->{getopt}, $shortbool);
            printf {$fh} "%15s  command-line: %s\n", '', shift(@lines);
            printf {$fh} "%15s                %s\n", '', $_ for @lines;
         }

         if (defined(my $env_name = $enamr->($option))) {
            printf {$fh} "%15s  environment : %s\n", '', $env_name;
         }

         printf {$fh} "%15s  default     : %s\n", '',
           $option->{default} // '*undef*'
           if exists $option->{default};
      } ## end for my $option ($options...)
      print {$fh} "\n";
   } ## end if ($options->@*)
   else {
      print {$fh} "This command has no options.\n\n";
   }

   if (my @children = get_children($self, $command)) {
      print {$fh} "Sub commands:\n", list_commands($self, \@children),
        "\n";
   }
   else {
      print {$fh} "no sub-commands\n\n";
   }
}

sub stock_SpecFromHash ($s, $cmd) {
   return $cmd if ref($cmd) eq 'HASH';
   return $s->{application}{commands}{$cmd} // undef;
}

sub stock_SpecFromHashOrModule ($s, $cmd) {
   return $cmd if ref($cmd) eq 'HASH';
   return $s->{application}{commands}{$cmd}
      //= $s->{factory}->($cmd, 'spec')->();
}

sub fetch_spec_for ($self, $command) {
   my $fetcher = $self->{application}{configuration}{specfetch}
      // \&stock_SpecFromHash;
   return $self->{factory}->($fetcher, 'specfetch')->($self, $command);
}

sub run ($application, $args) {
   $application = add_auto_commands(load_application($application));
   my $self = {
      application => $application,
      configs     => [],
      factory     => generate_factory($application->{factory} // {}),
      helpers     => {
         'print-commands' => \&print_commands,
         'print-help'     => \&print_help,
      },
      trail       => [['MAIN', $application->{commands}{MAIN}{name}]],
   };

   while ('necessary') {
      my $command = $self->{trail}[-1][0];
      my $spec    = fetch_spec_for($self, $command)
        or die "no definition for '$command'\n";

      $args = collect_options($self, $spec, $args);
      validate_configuration($self, $spec, $args);
      commit_configuration($self, $spec, $args);

      my ($subc, $alias) = fetch_subcommand($self, $spec, $args) or last;
      push $self->{trail}->@*, [$subc, $alias];
   } ## end while ('necessary')

   return execute($self, $args) // 0;
} ## end sub run

sub slurp ($file, $mode = '<:encoding(UTF-8)') {
   open my $fh, $mode, $file or die "open('$file'): $!\n";
   local $/;
   return <$fh>;
}

sub sources ($self, $spec, $args) {
   my $s = $spec->{sources}
     // $self->{application}{configuration}{sources}
     // \&stock_DefaultSources;
   $s = $self->{factory}->($s, 'sources')->() if 'ARRAY' ne ref $s;
   return $s->@*;
} ## end sub sources

sub stock_CmdLine ($self, $spec, $args) {
   my @args = $args->@*;
   my $goc  = $spec->{'getopt-config'}
     // default_getopt_config($self, $spec);
   require Getopt::Long;
   Getopt::Long::Configure('default', $goc->@*);

   my %option_for;
   my @specs = map {
      my $go = $_->{getopt};
      ref($go) eq 'ARRAY'
        ? ($go->[0] => sub { $go->[1]->(\%option_for, @_) })
        : $go;
     }
     grep { exists $_->{getopt} } ($spec->{options} // [])->@*;
   Getopt::Long::GetOptionsFromArray(\@args, \%option_for, @specs)
     or die "bailing out\n";

   # Check if we want to forbid the residual @args to start with a '-'
   my $strict = !$spec->{'allow-residual-options'};
   if ($strict && @args && $args[0] =~ m{\A -}mxs) {
      Getopt::Long::Configure('default', 'gnu_getopt');
      Getopt::Long::GetOptionsFromArray(\@args, {});
      die "bailing out\n";
   }

   return (\%option_for, \@args);
} ## end sub stock_CmdLine

sub stock_JsonFileFromConfig ($self, $spec, $args) {
   my $key = $spec->{'config-option'} // 'config';
   return {} if !exists($spec->{config}{$key});
   require JSON::PP;
   return JSON::PP::decode_json(slurp($spec->{config}{$key}));
} ## end sub stock_JsonFileFromConfig

sub stock_JsonFiles ($self, $spec, @ignore) {
   return merger($self, $spec)->(
      map {
         require JSON::PP;
         JSON::PP::decode_json(slurp($_));
        }
        grep { -e $_ } ($spec->{'config-files'} // [])->@*
   );
} ## end sub stock_JsonFiles

sub stock_Default ($self, $spec, @ignore) {
   return {
      map { '//=' . name_for_option($_) => $_->{default} }
      grep { exists $_->{default} } ($spec->{options} // [])->@*
   };
} ## end sub stock_Default

sub stock_Environment ($self, $spec, @ignore) {
   my $enamr = env_namer($self, $spec);
   return {
      map {
         my $en = $enamr->($_); # name of environment variable
         defined($en) && exists($ENV{$en})
            ? (name_for_option($_) => $ENV{$en}) : ();
        } ($spec->{options} // [])->@*
   };
} ## end sub stock_Environment

sub stock_NamEnv ($self, $cspec, $ospec) {
   my $aek = 'auto-environment';
   my $autoenv = exists $cspec->{$aek} ? $cspec->{$aek}
      : $self->{application}{configuration}{$aek} // undef;
   my $env = exists $ospec->{environment} ? $ospec->{environment}
      : $autoenv ? 1 : undef;
   return $env unless ($env // '') eq '1';
   my $appname = $self->{application}{configuration}{name} // '';
   my $optname = name_for_option($ospec);
   return uc(join '_', $appname, $optname);
}

sub stock_Parent ($self, $spec, @ignore) { $self->{configs}[-1] // {} }

sub stock_commands ($self, $config, $args) {
   die "this command does not support arguments\n" if $args->@*;
   my $target = get_descendant($self, $self->{trail}[-2][0], $args);
   print_commands($self, $target);
   return 0;
} ## end sub stock_commands

sub stock_factory ($executable, $default_subname = '', $opts = {}) {
   state $factory = sub ($executable, $default_subname) {
      my @prefixes =
          !defined $opts->{prefixes}       ? ()
        : 'ARRAY' eq ref $opts->{prefixes} ? ($opts->{prefixes}->@*)
        :                                    ($opts->{prefixes});
      push @prefixes, {'+' => 'AGAT::AppEaser#stock_'};
    SEARCH:
      for my $expansion_for (@prefixes) {
         for my $p (keys $expansion_for->%*) {
            next if $p ne substr $executable, 0, length $p;
            substr $executable, 0, length $p, $expansion_for->{$p};
            last SEARCH;
         }
      } ## end SEARCH: for my $expansion_for (...)

      # if it *still* "starts" with '=', it's "inline" Perl code
      return eval $executable if $executable =~ s{\A \s* = \s* }{}mxs;

      my ($package, $sname) = split m{\#}mxs, $executable;
      $sname = $default_subname unless defined $sname && length $sname;

      # first try to see if the sub is already available in $package
      if (my $s = $package->can($sname)) { return $s }

      # otherwise force loading of $package and retry
      (my $path = "$package.pm") =~ s{::}{/}gmxs;
      require $path;
      if (my $s = $package->can($sname)) { return $s }

      die "no '$sname' in '$package'\n";
   };
   state $cache = {};

   my $args;
   ($executable, $args) = ($executable->{executable}, $executable)
     if 'HASH' eq ref $executable;
   $executable = $cache->{$executable . ' ' . $default_subname} //=
     $factory->($executable, $default_subname)
     if 'CODE' ne ref $executable;
   return $executable unless $args;
   return sub { $executable->($args, @_) };
} ## end sub stock_factory

sub stock_help ($self, $config, $args) {
   print_help($self, get_descendant($self, $self->{trail}[-2][0], $args));
   return 0;
} ## end sub stock_help

sub stock_DefaultSources { [qw< +Default +CmdLine +Environment +Parent >] }

sub stock_SourcesWithFiles {
   [
      qw< +Default +CmdLine +Environment +Parent
         +JsonFileFromConfig +JsonFiles
        >
   ]
} ## end sub stock_SourcesWithFiles

sub validate_configuration ($self, $spec, $args) {
   my $from_spec = $spec->{validate};
   my $from_self = $self->{application}{configuration}{validate};
   my $validator;
   if (defined $from_spec && 'HASH' ne ref $from_spec) {
      $validator = $self->{factory}->($from_spec, 'validate');
   }
   elsif (defined $from_self && 'HASH' ne ref $from_self) {
      $validator = $self->{factory}->($from_self, 'validate');
   }
   else {    # use stock one
      $validator = \&params_validate;
   }
   $validator->($self, $spec, $args);
} ## end sub validate_configuration

exit run(
   $ENV{APPEASER} // {
      commands => {
         MAIN => {
            name        => 'main app',
            help        => 'this is the main app',
            description => 'Yes, this really is the main app',
            options     => [
               {
                  name        => 'foo',
                  description => 'option foo!',
                  getopt      => 'foo|f=s',
                  environment => 'FOO',
                  default     => 'bar',
               },
            ],
            execute => sub ($global, $conf, $args) {
               my $foo = $conf->{foo};
               say "Hello, $foo!";
               return 0;
            },
            'default-child' => '',    # run execute by default
         },
      },
   },
   [@ARGV]
) unless caller;

1;
