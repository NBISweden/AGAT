use strict;
use warnings;
use Test::More tests => 2;
use AGAT::AGAT;
use AGAT::Config;

# Ensure no local config file exists
my $config_file = 'agat_config.yaml';
unlink $config_file;

# Copy the default config locally
expose_config_file();

# Set verbose to 0 in config file
system(q{perl -i -pe 's/^verbose: 1/verbose: 0/' agat_config.yaml});

# Capture output from get_agat_config
my $output = '';
{
    open my $fh, '>', \$output;
    local *STDOUT = $fh;
    get_agat_config();
}

is($output, '', 'get_agat_config is quiet when verbose=0');

# Reset config to verbose 1
system(q{perl -i -pe 's/^verbose: 0/verbose: 1/' agat_config.yaml});
$output = '';
{
    open my $fh, '>', \$output;
    local *STDOUT = $fh;
    get_agat_config({ verbose => 0 });
}
is($output, '', 'CLI verbose=0 suppresses header');

unlink $config_file;
