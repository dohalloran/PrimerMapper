# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl PrimerMapper.t'

#########################

# change 'tests => 3' to 'tests => last_test_to_print';

use strict;
use warnings;

use Test::More tests => 3;
BEGIN { use_ok('PrimerMapper'), use_ok('PrimerMapperGraphics'), use_ok('PrimerMapperUtilities') };

#########################

