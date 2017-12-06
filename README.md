## PrimerMapper version 3.0
Author: Damien O'Halloran, The George Washington University

[![GitHub issues](https://img.shields.io/github/issues/dohalloran/PrimerMapper.svg)](https://github.com/dohalloran/PrimerMapper/issues)

## Installation
1. Download and extract the primerview.zip file  
`tar -xzvf PrimerMapper.zip`  
2. The extracted dir will be called PRIMERVIEW  
```cmd
  cd PrimerMapper  
  perl Makefile.PL  
  make  
  make test  
  make install 
```  
## Getting Started  
1. You must have `primer3` in your PATH  
[Click here to get primer3](https://sourceforge.net/projects/primer3/) 
2. Must have following [BioPerl](https://github.com/bioperl) modules:  
Bio::SeqIO    
Bio::Graphics  
Bio::SeqFeature::Generic  
WARNING: the option 'clean_up' deletes the '.txt' extension files generated from cwd   

## Usage 
Run as follows:  
```perl
  use strict;
  use warnings;
  use PRIMERVIEW;

  my $in_file = "test_seqs.fasta";

  my $tmp = PRIMERVIEW->new();
 
   $tmp->load_selections(  
      in_file         => $in_file, 
      single_view     => "1",   
      batch_view      => "1",      
      clean_up        => "1"   
   ); 
   
   $tmp->run_primerview();  
``` 

## Contributing
All contributions are welcome.

## Testing

PrimerMapper was successfully tested on:

- [x] Microsoft Windows 7 Enterprise ver.6.1
- [x] Linux Mint 64-bit ver.18

## Support
If you have any problem or suggestion please open an issue [here](https://github.com/dohalloran/PrimerMapper/issues).

## License 
GNU GENERAL PUBLIC LICENSE





