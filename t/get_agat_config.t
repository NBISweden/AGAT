use strict;
use warnings;
use Test::More tests => 1;
use AGAT::AGAT;
use AGAT::Config;

# Ensure no local config file exists
my $config_file = 'agat_config.yaml';
unlink $config_file;

# Copy the default config locally
expose_config_file();

# Set verbose to 0
system(q{perl -i -pe 's/^verbose: 1/verbose: 0/' agat_config.yaml});

# Capture output from get_agat_config
my $output = '';
{
    open my $fh, '>', \$output;
    local *STDOUT = $fh;
    get_agat_config();
}

is($output, '', 'get_agat_config is quiet when verbose=0');

unlink $config_file;
