# CRISPRLeaderExtractor

Author: Murat Buyukyoruk

    CRISPRLeaderExtractor help:

This script is developed to fetch flank sequences of gene/CRISPR loci of interest by using the fasta file with provided position and strand information. 

SeqIO and Seq packages from Bio is required to fetch flank sequences. Additionally, tqdm is required to provide a progress bar since some multifasta files can contain long and many sequences.

Syntax:

    python CRISPRLeaderExtractor.py -i demo.fasta -o demo_gene_flanks.fasta -d demo_CRISPRDetect.txt

    OR

    python CRISPRLeaderExtractor.py -i demo.fasta -o demo_gene_flanks.fasta -d demo_CRISPRDetect.txt -f 200 -c Y -r Y -s I-E

CRISPRLeaderExtractor dependencies:

Bio module, SeqIO and Seq available in this package     refer to https://biopython.org/wiki/Download
	
tqdm                                                    refer to https://pypi.org/project/tqdm/

Input Paramaters (REQUIRED):
----------------------------
	-i/--input		FASTA			Specify a fasta file. FASTA file requires headers starting with accession number. (i.e. >NZ_CP006019 [fullname])

	-o/--output		Output file	        Specify a output file name that should contain fetched sequences.

	-d/--data		Dataframe		Specify a list of accession (Accession only). Each accession should be included in a new line (i.e. generated with Excel spreadsheet). Script works with or without '>' symbol before the accession.

Parameters [optional]:
----------------------
	-f/--flank		200         This is the default length of flanks that is fetched.

	-c/--circular           N           This is the default option that is assuming the sequence is not circular. Type "Y" insted if you know it is cirgular genome.

	-r/--repeat	        Y           Include first repeat of not? (Default: Y).

	-s/--subtype	        I-F         Specify CRISPR subtype assigned by CRISPRDetect.

Basic Options:
--------------
	-h/--help		HELP			Shows this help text and exits the run.

Output header will contain array no, original accession number, positions and fullname (if included in original fasta), observed stand.
