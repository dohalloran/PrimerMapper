## PrimerMapper version 3.2
Author: Damien O'Halloran, The George Washington University

[![GitHub issues](https://img.shields.io/github/issues/dohalloran/PrimerMapper.svg)](https://github.com/dohalloran/PrimerMapper/issues)

## What PrimerMapper Does

PrimerMapper is a GUI-based tool for PCR primer design and visualization. It supports:

1. Traditional primer design against one or more FASTA sequences
2. SNP-flanking primer design using IUPAC degenerate base notation
3. Allele-specific primer design for all ten IUPAC ambiguity codes
4. Graphical output showing primer positions across each sequence and a concatenated single-view
5. Remote sequence retrieval from NCBI GenBank and dbSNP
6. Remote BLAST of designed primers against multiple NCBI databases

## Dependencies

### Required external tool
- **primer3** — the primer design engine that PrimerMapper wraps

### Required Perl modules
- `Tkx` — GUI framework
- `LWP::Simple` — remote sequence retrieval from NCBI
- `Bio::SeqIO` — FASTA parsing (BioPerl)
- `Bio::Graphics` — graphical primer maps (BioPerl)
- `Bio::SeqFeature::Generic` — sequence feature objects (BioPerl)

## Installation

### Linux (Debian / Ubuntu / Mint)

```bash
# 1. Install primer3 and Tk
sudo apt-get install primer3 tk-dev

# 2. Install Perl modules
cpan Tkx LWP::Simple Bio::SeqIO Bio::Graphics

# Alternatively, use the BioPerl apt package:
sudo apt-get install libbioperl-perl
```

### macOS (Homebrew)

```bash
brew install perl primer3 gd tcl-tk
cpan App::cpanminus
cd PrimerMapper
cpanm --notest Tkx
cpanm --installdeps .
```

If your shell still uses Apple's system Perl after installing Homebrew Perl, restart Terminal and check `which perl`.
On Apple Silicon Macs, Homebrew Perl is usually under `/opt/homebrew/bin/perl`; on Intel Macs it is usually under `/usr/local/bin/perl`.

### Windows (Strawberry Perl)

1. Install [Strawberry Perl](https://strawberryperl.com/)
2. Install [primer3](https://sourceforge.net/projects/primer3/) and add it to your PATH
3. From the Strawberry Perl CPAN shell:

```cmd
cpan Tkx LWP::Simple Bio::SeqIO Bio::Graphics
```

### Install PrimerMapper itself

```bash
git clone https://github.com/dohalloran/PrimerMapper.git
cd PrimerMapper/PrimerMapper
perl Makefile.PL
make
make test
make install
```

## Running PrimerMapper

From the top-level directory:

```bash
perl PrimerMapper_driver.pl
```

This opens the GUI. No command-line arguments are needed.

### macOS launch command used for the dependency-fixed build

If you installed PrimerMapper with Homebrew Perl and the isolated local library
used in the dependency-fixed package, launch it like this:

```bash
conda deactivate 2>/dev/null || true

export PM_PERL="$HOME/perl5-primermapper"
export BREW_PERL="$(brew --prefix perl)/bin/perl"
eval "$("$BREW_PERL" -I"$PM_PERL/lib/perl5" -Mlocal::lib="$PM_PERL")"

cd ~/Desktop/PrimerMapper-fixed
"$BREW_PERL" PrimerMapper_driver.pl
```

If the GUI opens but graphics fail with a GD/PNG error, confirm PNG support:

```bash
"$BREW_PERL" -MGD -e 'print "GD version: $GD::VERSION\n"; print GD::Image->can("png") ? "PNG OK\n" : "PNG MISSING\n";'
```

The expected result is `PNG OK`.

## Usage

### Standard primer design
1. Click **Load File** to select a multi- or single-entry FASTA file, or enter accession numbers (comma-separated, e.g. `NM_001307521.1,FJ455617.1`) in the text box and click **Get Sequences** to fetch directly from NCBI.
2. Adjust primer design parameters in the **PRIMER DESIGN** panel (product size range, primer length, GC%, Tm, GC clamp, salt concentration).
3. Click **Design Primers and Graphics** to run primer3 and generate output files.

Regular FASTA input must start with a descriptor line beginning with `>`:

```text
>example_sequence
ATGCGTACGTACGTACGTACGTACGTACGTAC
```

### SNP primer design
1. Click **Load SNP File** to select a FASTA file where the SNP position is encoded as a degenerate IUPAC base (R, Y, S, W, K, M, B, D, H, or V), or enter rs numbers and click **Get SNP Sequences**.
2. Set the SNP design parameters (flanking window sizes, primer length, GC%, Tm, minimum distance from SNP).
3. Click **Design SNP Primers and Graphics**.

SNP FASTA input must also start with `>`, and the header must contain a
parseable SNP coordinate in the form `pos=N.len`. The base at that coordinate
must be an IUPAC ambiguity code:

```text
>example_snp pos=25.len=80
ATGCGTACGTACGTACGTACGTACRACGTACGTACGTACGTACGTACGTAC
```

In this example, position 25 contains `R`, meaning A or G. PrimerMapper will now
show a message if the SNP coordinate is missing, outside the sequence, or does
not point to an ambiguity code.

### Input validation and user messages

PrimerMapper now checks common problems before drawing graphics:

- Missing or unreadable FASTA/SNP FASTA files
- FASTA files that do not begin with `>`
- Invalid product-size, primer-length, GC%, Tm, salt, or GC-clamp settings
- Minimum product size larger than maximum product size
- SNP files missing `pos=N.len`
- SNP coordinates outside the sequence
- SNP coordinate not matching an IUPAC ambiguity code
- SNP sequences shorter than the 100 bp SNP product-size minimum
- Cases where primer3 runs successfully but returns no primer pairs

When one of these problems occurs, the GUI should now show a message explaining
what to change instead of printing a Perl/Tkx stack trace.

### BLAST
1. Paste your primer sequences in FASTA format into the BLAST text box.
2. Select a database (Nucleotide, Human Genomic, etc.).
3. Click **BLAST**. Results are written to a text file in the working directory.

### Clean-up
Click **Clean-up** to remove temporary files generated during a run. Output graphics (`.png`/`.gif`) and TSV files are preserved; only intermediate `.txt` working files are deleted.

## Output Files

| File | Description |
|---|---|
| `GRAPHIC_<id>.png` | Per-sequence primer position map |
| `single_view.png` | Concatenated view of all sequences and primers |
| `primers.tsv` | Primer sequences with coordinates |
| `SNP_primers.tsv` | SNP primer sequences with coordinates |
| `Allele_specific_primers.tsv` | Allele-specific oligos with Tm and GC% |
| `blast_primers.fasta` | Primer FASTA submitted to BLAST |

## Changes in version 3.1

- Fixed `calcTm` reading from `$_` instead of the `$sequence` argument (incorrect Tm values for all allele-specific primers)
- Fixed NCBI URLs from `http://` to `https://` (sequence retrieval was silently failing)
- Fixed filehandle not closed before `system()` call to primer3 (potential truncated input on Windows)
- Fixed typo in temp filenames (`specificty` → `specificity`; trailing space in SNP filename)
- Fixed output graphics being opened in append mode — re-runs no longer corrupt PNG files
- Fixed `clean_up` being too aggressive (previously deleted all `.txt` in cwd, now targets only PrimerMapper temp files)
- Fixed allele-specific primer loop opening output file on every kmer iteration (file descriptor leak)
- Fixed `"Reverse:"` label incorrectly applied to forward primers for IUPAC type D
- Fixed `graphics_all_primers_SNP` crash when a sequence produced zero primers
- Fixed `PrimerMapper_driver.pl` to use `FindBin` so it runs from any directory
- Added input validation dialogs (missing file, missing BLAST database, empty BLAST text box)
- Source SNP FASTA file is no longer deleted after preprocessing

## Changes in version 3.2

- Added GUI error messages for invalid regular primer-design settings, including product-size range and primer-length range
- Added GUI error messages when regular FASTA or SNP FASTA input is missing, empty, unreadable, or not FASTA format
- Added SNP-specific validation for `pos=N.len`, SNP coordinate bounds, and IUPAC ambiguity code at the declared SNP coordinate
- Added a clear message when primer3 returns no regular or SNP primer pairs, including primer3 explanation lines when available
- Guarded regular and SNP graphics generation so missing `single_viewer.txt` / `single_viewer_SNP.txt` no longer causes a crash
- Added NCBI retrieval validation so accession/rs lookup failures do not write invalid FASTA files silently
- Added README launch instructions for the Homebrew Perl/local-lib macOS setup

## Testing

PrimerMapper has been tested on:

- [x] Linux Mint 64-bit ver. 18
- [x] Linux Mint 64-bit ver. 21+

## Contributing
All contributions are welcome.

## Support
If you have any problem or suggestion please open an issue [here](https://github.com/dohalloran/PrimerMapper/issues).

## License
GNU GENERAL PUBLIC LICENSE
