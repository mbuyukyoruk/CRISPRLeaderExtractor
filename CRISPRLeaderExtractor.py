
import argparse
import os
import sys
import subprocess
import re
import textwrap

orig_stdout = sys.stdout

try:
    from Bio import SeqIO
except ImportError as e:
    print("SeqIO module is not installed! Please install SeqIO and try again.")
    sys.exit()

try:
    from Bio.Seq import Seq
except ImportError as e:
    print("Seq module is not installed! Please install Seq and try again.")
    sys.exit()

try:
    import tqdm
except ImportError as e:
    print("tqdm module is not installed! Please install tqdm and try again.")
    sys.exit()

parser = argparse.ArgumentParser(prog='python CRISPRLeaderExtractor.py',
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 epilog=textwrap.dedent('''\

      	Author: Murat Buyukyoruk

        CRISPRLeaderExtractor help:

This script is developed to fetch flank sequences of gene/CRISPR loci of interest by using the fasta file with provided position and strand information. 

SeqIO and Seq packages from Bio is required to fetch flank sequences. Additionally, tqdm is required to provide a progress bar since some multifasta files can contain long and many sequences.

Syntax:

    python CRISPRLeaderExtractor.py -i demo.fasta -o demo_gene_flanks.fasta -d flank_info_dataframe.txt

    OR

    python CRISPRLeaderExtractor.py -i demo.fasta -o demo_gene_flanks.fasta -d demo_flank_info_dataframe.txt -f 200 -c Y -r Y -s I-E

CRISPRLeaderExtractor dependencies:

	Bio module, SeqIO and Seq available in this package     refer to https://biopython.org/wiki/Download
	
	tqdm                                                    refer to https://pypi.org/project/tqdm/

Input Paramaters (REQUIRED):
----------------------------
	-i/--input		FASTA			Specify a fasta file. FASTA file requires headers starting with accession number. (i.e. >NZ_CP006019 [fullname])

	-o/--output		Output file	    Specify a output file name that should contain fetched sequences.

	-d/--data		Dataframe		Specify a list of accession (Accession only). Each accession should be included in a new line (i.e. generated with Excel spreadsheet). Script works with or without '>' symbol before the accession.

Parameters [optional]:
----------------------
	-f/--flank		200			    This is the default length of flanks that is fetched.

	-c/--circular	N			    This is the default option that is assuming the sequence is not circular. Type "Y" insted if you know it is cirgular genome.

	-r/--repeat	    Y			    Include first repeat of not? (Default: Y).

	-s/--subtype	I-F			    Specify CRISPR subtype assigned by CRISPRDetect.

Basic Options:
--------------
	-h/--help		HELP			Shows this help text and exits the run.

Output header will contain array no, original accession number, positions and fullname (if included in original fasta), observed stand.

      	'''))
parser.add_argument('-i', '--input', required=True, type=str, dest='filename',
                    help='Specify a fastafile to fetch regions from.\n')
parser.add_argument('-o', '--output', required=True, dest='out',
                    help='Specify a output file name.\n')
parser.add_argument('-d', '--data', required=True, dest='data',
                    help='Specify a dataframe with accession, prot_accession/array_no, start, stop, strand info in that order. Leave NA if not known.\n')
parser.add_argument('-f', '--flank', type=int, required=False, default=200, dest='flank',
                    help='Specify length of flanking sequence to fetch.\n')
parser.add_argument('-c', '--circular', type=str, required=False, default="N", dest='circular',
                    help='Circular genome? Y/N (Default: N). Split circluar and linear sequences if you have a mixture.\n')
parser.add_argument('-r', '--repeat', type=str, required=False, default="Y", dest='repeat',
                    help='Get repeat? Y/N (Default: Y)\n')
parser.add_argument('-s', '--subtype', type=str, required=False, default="All", dest='subtype',
                    help='Specify subtype (i.e. IE or I-E) or type All for all subtypes (Default: All)\n')

results = parser.parse_args()
filename = results.filename
out = results.out
data = results.data
flank = results.flank
circ = results.circular
repeat = results.repeat
subtype_search = results.subtype

if circ in ("Y", "N", "y", "n"):
    pass
else:
    print("""\nWrong argument! Type "Y" for Yes(cicular) or "N" for No (linear).\n""")

    os.system('ORF_flank_fetch.py -h')
    sys.exit()

# Generate summary Dataframe of CRISPRDetect output

df_out = f'{data}.df'

os.system("> " + df_out)

f = open(df_out, 'a')
sys.stdout = f

print(
    "Array_no\tAccession\tName\tStart\tStop\tStrand\tSubtype\tRepeat_seq\tScore")

proc = subprocess.Popen("grep -c '>' " + data, shell=True, stdout=subprocess.PIPE, text=True)
length = int(proc.communicate()[0].split('\n')[0])

with tqdm.tqdm(range(length)) as pbar_sum:
    pbar_sum.set_description('Creating summary file...' + data)
    with open(data, 'r') as file:
        for line in file:
            line_arr = line.split()
            if (len(line_arr) != 0):
                if line_arr[0] == "Array":
                    array_no = line_arr[1]
                if "# Array family :" in line:
                    subtype = line.split("# Array family : ")[1].split('\n')[0]
                if "# Summary:" in line:
                    pbar_sum.update()
                    arr = line.split("ID_START_STOP_DIR: ")[1].split(';')
                    name = arr[0]
                    fullname = name.split("-")[0]
                    acc = name.split()[0]
                    low_bound = name.split("-")[-3]
                    high_bound = name.split("-")[-2]
                    ori = name.split("-")[-1]
                    dr_seq = arr[1].split(":")[1]
                    score = arr[-2].split(":")[1]

                    print(
                        acc + '_Array_' + array_no + '\t' + acc + '\t' + fullname + '\t' + low_bound + '\t' + high_bound + '\t' + ori + '\t' + subtype + '\t' + dr_seq + "\t" + score)

sys.stdout = orig_stdout

id_list = []
seq_list= []

os.system('> ' + out)

proc = subprocess.Popen("grep -c '>' " + filename, shell=True, stdout=subprocess.PIPE, text=True)
length = int(proc.communicate()[0].split('\n')[0])

with tqdm.tqdm(range(length), desc='Importing...') as pbar:
    for record in SeqIO.parse(filename, "fasta"):
        pbar.update()
        id_list.append(record.id)
        seq_list.append(record.seq)

proc = subprocess.Popen("wc -l < " + df_out, shell=True, stdout=subprocess.PIPE, text=True)
length = int(proc.communicate()[0].split('\n')[0])

with tqdm.tqdm(range(length), desc='Fetching...') as pbar:
    with open(df_out, 'r') as file:
        f = open(out, 'a')
        sys.stdout = f
        for line in file:
            pbar.update()
            if "Accession" in line:            
                arr = line.split('\t')
                index_Array_no = arr.index("Array_no")
                index_Accession = arr.index("Accession")
                index_Name = arr.index("Name")
                index_Start = arr.index("Start")
                index_Stop = arr.index("Stop")
                index_Strand = arr.index("Strand")
                index_Subtype = arr.index("Subtype")
                index_Repeat_seq = arr.index("Repeat_seq")
                index_Score = arr.index("Score\n")

            if "Accession" not in line:
                array = line.split('\t')[index_Array_no]
                acc = line.split('\t')[index_Accession]
                name = line.split('\t')[index_Name]
                start = line.split('\t')[index_Start]
                stop = line.split('\t')[index_Stop]
                strand = line.split('\t')[index_Strand]
                subtype = line.split('\t')[index_Subtype].replace('"','')
                if subtype != 'NA':
                    if "Matched known repeat from this family" in subtype:
                        subtype = subtype.split()[0]
                    else:
                        subtype = 'NA'

                if repeat == "Y":
                    repeat_seq = line.split('\t')[index_Repeat_seq]
                elif repeat == "N":
                    repeat_seq_ori = line.split('\t')[index_Repeat_seq]
                    repeat_seq = ""

                score = line.split('\t')[index_Score].replace("\n","")

                try:
                    ind = id_list.index(acc)
                except:
                    sys.stdout = orig_stdout
                    print('\n\n' + acc + ' is missing in the fasta file. Please make sure that you are using the correct fasta.\n')
                    sys.exit()

                seq_raw = seq_list[ind]
                if subtype.lower().replace('-','') == subtype_search.lower().replace('-','') or subtype_search.lower().replace('-','') == 'all':
                    if strand == "F":
                        if int(start) - flank > 0:
                            print('>' + array + " | " + start + "-" + stop + " | " + strand + " | " + subtype + " | " + score + " | " + acc + " | " + name + " | " +str(int(start) - flank) + '-' + str(int(start) - 1 + len(repeat_seq)))
                            seq = str(Seq(str(seq_raw[int(start) - 1 - (flank):int(start) - 1 + len(repeat_seq)])))
                            print(re.sub("(.{60})", "\\1\n", seq, 0, re.DOTALL))
                        elif int(start) - flank <= 0:
                            if circ == 'N':
                                print('>' + array + " | " + start + "-" + stop + " | " + strand + " | " + subtype + " | " + score + " | " + acc + " | " + name + " | "  + str(1) + '-' + str(int(start) - 1 + len(repeat_seq)))
                                seq = str(Seq(str(seq_raw[0:int(start) - 1 + len(repeat_seq)])))
                                print(re.sub("(.{60})", "\\1\n", seq, 0, re.DOTALL))

                            elif circ == 'Y':
                                remaining = (int(start) - flank) * (-1)
                                print('>' + array + " | " + start + "-" + stop + " | " + strand + " | " + subtype + " | " + score + " | " + acc + " | " + name + " | "  + str(len(seq_raw)-remaining) + '-' + str(len(seq_raw)) + '~' + str(1) + '-' + str(int(start) - 1 + len(repeat_seq)))
                                seq_part = str(Seq(str(seq_raw[0:int(start) - 1+ len(repeat_seq)])))
                                seq_round = str(Seq(str(seq_raw[len(seq_raw) - ((flank + 1) - int(start)):len(seq_raw)])))
                                seq = seq_round + seq_part
                                print(re.sub("(.{60})", "\\1\n", seq, 0, re.DOTALL))

                    elif strand == "R":
                        if int(stop) + flank <= len(seq_raw):
                            print('>' + array + " | " + start + "-" + stop + " | " + strand + " | " + subtype + " | " + score + " | " + acc + " | " + name + " | "  + str(int(stop) - len(repeat_seq)) + '-' + str(int(stop) + flank -1))
                            seq = str(Seq(str(seq_raw[int(stop)-len(repeat_seq)-1:int(stop)-1 + flank])).reverse_complement())
                            print(re.sub("(.{60})", "\\1\n", seq, 0, re.DOTALL))

                        elif int(stop) + flank > len(seq_raw):
                            if circ == 'N':
                                print('>' + array + " | " + start + "-" + stop + " | " + strand + " | " + subtype + " | " + score + " | " + acc + " | " + name + " | "  + str(int(stop) - len(repeat_seq)) + '-' + str(len(seq_raw)))
                                seq = str(Seq(str(seq_raw[int(stop)-len(repeat_seq)-1:len(seq_raw)])).reverse_complement())
                                print(re.sub("(.{60})", "\\1\n", seq, 0, re.DOTALL))

                            elif circ == 'Y':
                                remaining = int(stop) + flank - len(seq_raw)-1
                                print('>' + array + " | " + start + "-" + stop + " | " + strand + " | " + subtype + " | " + score + " | " + acc + " | " + name + " | "  + str(int(stop) - len(repeat_seq)) + '-' + str(len(seq_raw)) + '~' + str(1) + '-' + str(remaining))
                                seq_part = str(Seq(str(seq_raw[int(stop) - len(repeat_seq) - 1:len(seq_raw)])).reverse_complement())
                                seq_round = str(Seq(str(seq_raw[0:remaining])).reverse_complement())
                                seq = seq_round + seq_part
                                print(re.sub("(.{60})", "\\1\n", seq, 0, re.DOTALL))
                    if strand == 'NA':
                        if int(start) - flank > 0:
                            print('>' + array + " | " + start + "-" + stop + " | " + strand + "_F | " + subtype + " | " + score + " | " + acc + " | " + name + " | "  + str(int(start) - flank) + '-' + str(int(start) - 1 + len(repeat_seq)))
                            seq = str(Seq(str(seq_raw[int(start) - 1 - (flank):int(start) - 1 + len(repeat_seq)])))
                            print(re.sub("(.{60})", "\\1\n", seq, 0, re.DOTALL))

                        elif int(start) - flank <= 0:
                            if circ == 'N':
                                print('>' + array + " | " + start + "-" + stop + " | " + strand + "_F | " + subtype + " | " + score + " | " + acc + " | " + name + " | "  + str(1) + '-' + str(int(start) - 1 + len(repeat_seq)))
                                seq = str(Seq(str(seq_raw[0:int(start) - 1 + len(repeat_seq)])))
                                print(re.sub("(.{60})", "\\1\n", seq, 0, re.DOTALL))

                            elif circ == 'Y':
                                remaining = (int(start) - flank) * (-1)
                                print('>' + array + " | " + start + "-" + stop + " | " + strand + "_F | " + subtype + " | " + score + " | " + acc + " | " + name + " | " + str(len(seq_raw)-remaining) + '-' + str(len(seq_raw)) + '~' + str(1) + '-' + str(int(start) - 1 + len(repeat_seq)))
                                seq_part = str(Seq(str(seq_raw[0:int(start) - 1+ len(repeat_seq)])))
                                seq_round = str(Seq(str(seq_raw[len(seq_raw) - ((flank + 1) - int(start)):len(seq_raw)])))
                                seq = seq_round + seq_part
                                print(re.sub("(.{60})", "\\1\n", seq, 0, re.DOTALL))

                        if int(stop) + flank <= len(seq_raw):
                            print('>' + array + " | " + start + "-" + stop + " | " + strand + "_R | " + subtype + " | " + score + " | " + acc + " | " + name + " | "  + str(int(stop) - len(repeat_seq)) + '-' + str(int(stop) + flank -1))
                            seq = str(Seq(str(seq_raw[int(stop)-len(repeat_seq)-1:int(stop)-1 + flank])).reverse_complement())
                            print(re.sub("(.{60})", "\\1\n", seq, 0, re.DOTALL))

                        elif int(stop) + flank > len(seq_raw):
                            if circ == 'N':
                                print('>' + array + " | " + start + "-" + stop + " | " + strand + "_R | " + subtype + " | " + score + " | " + acc + " | " + name + " | "  + str(int(stop) - len(repeat_seq)) + '-' + str(len(seq_raw)))
                                seq = str(Seq(str(seq_raw[int(stop)-len(repeat_seq)-1:len(seq_raw)])).reverse_complement())
                                print(re.sub("(.{60})", "\\1\n", seq, 0, re.DOTALL))

                            elif circ == 'Y':
                                remaining = int(stop) + flank - len(seq_raw)-1
                                print('>' + array + " | " + start + "-" + stop + " | " + strand + "_R | " + subtype + " | " + score + " | " + acc + " | " + name + " | "  + str(int(stop) - len(repeat_seq)) + '-' + str(len(seq_raw)) + '~' + str(1) + '-' + str(remaining))
                                seq_part = str(Seq(str(seq_raw[int(stop) - len(repeat_seq) - 1:len(seq_raw)])).reverse_complement())
                                seq_round = str(Seq(str(seq_raw[0:remaining])).reverse_complement())
                                seq = seq_round + seq_part
                                print(re.sub("(.{60})", "\\1\n", seq, 0, re.DOTALL))
                else:
                    continue