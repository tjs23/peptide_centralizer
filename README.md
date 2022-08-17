# peptide_centralizer
Extract protein sub-sequences around identified sites of interest, e.g. phosphorylation sites.

## Command Line Usage
```
centralise_peptides.py [-h] [-o OUT_FILE_ROOT] [-a ALIGN_PATH]
                       [-na NUM_ORTHO_ALIGN] [-w HALF_WIDTH]
                       [-p PAD_CHAR] [-ci COLUMN] [-cp COLUMN]
                       CSV_FILE
```
### Positional Arguments
```
   CSV_FILE             Input CSV file containing mass spec data, must contain columns for protein ID with site location
                        and for peptide sub-sequence.
```
### Optional Arguments
```
  -h, --help            Show this help message and exit

  -o OUT_FILE_ROOT      Output file path prefix (i.e. no file extension) for writing output CSV and FASTA format files.
                        If unspecified a default will be based on the input CSV file.

  -a ALIGN_PATH         Optional alignment file path, in concatenated FASTA format, as output by "align_orthologues.py".

  -na NUM_ORTHO_ALIGN   Number of aligned ortholgue sequences (from alignment file) to use for each peptide; to expand
                        range of sequences used. Default: 5

  -w HALF_WIDTH         Number of amino acid residues to add to both sides the central site to make each subsequence.
                        Default: 7

  -p PAD_CHAR           Padding/extension character to pad the start/end of protein sequences. Default: X

  -ci COLUMN            Column in the CSV file containing UniProt protein IDs and site of interest in the format
                        {PROT_ID}_{AA}{RES_NUM}___{COUNT}, e.g. "Q9Z2D1_S631___1". Column identifier may be a number,
                        starting from 1, of a letter starting from A. Default 1

  -cp COLUMNN           Column in the CSV file containing MS identified peptide subsequence(s). Column identifier may
                        be a number, starting from 1, of a letter starting from A. Default 5
```
