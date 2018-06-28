## PrimerMapper version 3.0
Author: Damien O'Halloran, The George Washington University

[![GitHub issues](https://img.shields.io/github/issues/dohalloran/PrimerMapper.svg)](https://github.com/dohalloran/PrimerMapper/issues)

## Installation
1. Clone the PrimerMapper repo  
`git clone https://github.com/dohalloran/PrimerMapper.git`  
2. The extracted dir will be called PrimerMapper  
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
or on Linux `sudo apt-get install primer3` 
2. Must have following [BioPerl](https://github.com/bioperl) modules:  
Bio::SeqIO    
Bio::Graphics  
Bio::SeqFeature::Generic  
WARNING: the option 'clean_up' deletes the '.txt' extension files generated from cwd   

## Usage 
From the commandline, type:  
```cmd
  perl PrimerMapper_driver.pl
``` 

## Contributing
All contributions are welcome.

## Testing

PrimerMapper was successfully tested on:

- [x] Linux Mint 64-bit ver.18

## Support
If you have any problem or suggestion please open an issue [here](https://github.com/dohalloran/PrimerMapper/issues).

## License 
GNU GENERAL PUBLIC LICENSE





