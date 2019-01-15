#!/home/linking/.pyenv/versions/3.5.1/bin/python
import re
import os
import argparse
import time
from Bio import SeqIO


# crispr.py
# Version 0.1
# 12.28.18
# by linxingchen


localtime = time.asctime(time.localtime(time.time()))
pwd = os.getcwd()

parser = argparse.ArgumentParser(description="This script gets all spacers from both paired reads, " +
                                             "when providing scaffolds, proteins predicted from scaffolds" +
                                             " and scaffolds mapping sam file.")
requiredNamed = parser.add_argument_group('required named arguments')
requiredNamed.add_argument("-f", "--fasta", type=str, help="The scaffolds file", required=True)
requiredNamed.add_argument("-s", "--sam", type=str, help="Scaffolds mapping sam", required=True)
requiredNamed.add_argument("-p", "--protein", type=str, help="Predicted proteins of scaffolds in AA format",
                           required=True)
parser.add_argument("-q", "--query", type=str, help="Scaffolds for spacer targeting search, will search the fasta "
                                                    "file if this is not provided")
parser.add_argument("-g", "--gene", type=str, help="Predicted genes in DNA or AA format of query scaffolds, "
                                                   "provide only when query is given")
args = parser.parse_args()


# Open scaffolds file
fasta_dict = {}

with open(args.fasta, 'r') as f:
    for record in SeqIO.parse(f, "fasta"):
        header = str(record.id).strip()
        seq = str(record.seq)
        fasta_dict[header] = seq


# Log file
log = open("{0}.crispr.log".format(args.fasta), 'w+')
print("The start time is: {0}".format(localtime), flush=True, file=log)
print('\n', end='', flush=True, file=log)

# Check if output directory exists or not
if os.path.exists("{0}_crispr".format(args.fasta)):
    print('Warning: Output directory exists, please rename it!')
    exit()
else:
    os.mkdir("{0}_crispr".format(args.fasta))

# Open protein file and run hmm for cas protein
with open(args.protein, 'r') as p:
    os.system("/home/davidbur/scripts/hmmsearchTable {0} "
              "/data10/Rifle/202_genomes/hmm/TIGRFAMs_15.0_HMM/TIGR00287.HMM 6 --cut_tc "
              "> {0}.crispr.Cas1.hmm.txt".format(args.protein))

print("Step 1: hmm search for cas proteins is finished.", flush=True, file=log)
print('\n', end='', flush=True, file=log)


# Parse hmm for scaffolds with cas protein(s)
cas = []
cas_hmm = open("{0}.crispr.Cas1.hmm.txt".format(args.protein), 'r')
for line in cas_hmm:
    if line.startswith('targetName') is False:
        line = line.strip().split()
        line[0] = line[0].split('_')
        scaf = '_'.join(line[0][0:int(len(line[0]))-1])
        cas.append(scaf)


# Create a fasta file to store scaffolds with cas protein(s)
cas_scaf = open("{0}.crispr.with_cas_proteins".format(args.fasta), 'w')
for scaffold in set(cas):
    print(">" + scaffold, file=cas_scaf)
    print(fasta_dict[scaffold], file=cas_scaf)

cas_scaf.close()

print("Step 2: scaffolds with cas proteins have been extracted.", flush=True, file=log)
print('\n', end='', flush=True, file=log)


# Run minced for scaffolds
os.system("minced -minSL 17 -gffFull {0} >{0}.crispr.minced.gffFull.txt".format(args.fasta))

print("Step 3: minced prediction of repeat and spacers on scaffolds is finished.", flush=True, file=log)
print('\n', end='', flush=True, file=log)


# Get direct repeat sequences from minced results
scaffold2repeat_output = open("{0}.crispr.scaffold2repeat.txt".format(args.fasta), 'w+')

minced_scaffold = list()
scaffold2repeat = {}

minced = open("{0}.crispr.minced.gffFull.txt".format(args.fasta), 'r')
for line in minced.readlines():
    line = line.strip().split()
    if line[0].startswith('#') is False and line[8][-3:] == 'DR2':
        if line[0] not in list(scaffold2repeat.keys()):
            print(line[0] + '\t' + fasta_dict[line[0]][int(line[3])-1:int(line[4])], file=scaffold2repeat_output)
            scaffold2repeat[line[0]] = fasta_dict[line[0]][int(line[3])-1:int(line[4])]
            minced_scaffold.append(line[0])
        elif line[0] in list(scaffold2repeat.keys()):
            if fasta_dict[line[0]][int(line[3])-1:int(line[4])] == scaffold2repeat[line[0]]:
                pass
            else:
                print(line[0] + '\t' + fasta_dict[line[0]][int(line[3])-1:int(line[4])], file=scaffold2repeat_output)
                minced_scaffold.append(line[0])


minced.close()
scaffold2repeat_output.close()
print("Step 4: direct repeat sequences have been extracted.", flush=True, file=log)
print('\n', end='', flush=True, file=log)

# For those scaffolds with two or more direct repeat sequences reported by minced, save them to an independent
# file and their spacers will not be saved to the output *spacers.fasta file (see 'single_array_scaffold' below).
multiple_DR = open("{0}.crispr.multiple_potential_arrays.txt".format(args.fasta), 'w+')
multiple_dr = list()

for item in set(minced_scaffold):
    if minced_scaffold.count(item) > 1:
        multiple_dr.append(item)

scaffold2repeat_output = open("{0}.crispr.scaffold2repeat.txt".format(args.fasta), 'r')
for line in scaffold2repeat_output:
    line = line.strip().split()
    if line[0] in multiple_dr:
        print(line[0] + '\t' + line[1], file=multiple_DR)

print("Step 5: scaffolds may have multiple potential arrays (if presented) have been extracted", flush=True, file=log)
print('\n', end='', flush=True, file=log)


# Function to check if a nucleotide sequence contains other element excluding A, T, C, G.
def check_dna(dna):
    bp = list()
    for letter in dna:
        bp.append(letter)
        return set(dna) <= {'A', 'T', 'C', 'G', 'N'}


# Function to get the reverse complementary sequence of a read.
def complement_reverse(dna):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    bases = list(dna)
    bases = [complement[base] for base in bases]
    return ''.join(bases)[::-1]


# The 1st dict to save scaffold mapped reads, the 2nd dict saves the spacers from mapped.reads. These information
# will be retrieved from the provided sam file.
scaffold2reads = {}
scaffold2spacers = {}

for scaffold in minced_scaffold:
    scaffold2reads[scaffold] = list()
    scaffold2spacers[scaffold] = list()

# This puts all mapped reads into the list of each scaffold in the dict of scaffold2reads.
no_mapped_reads = open("{0}.crispr.without.mapping.txt".format(args.fasta), 'w+')
sam_scaffold = list()

with open(args.sam, 'r') as sam:
    for line in sam.readlines():
        if line.startswith('@') is False:
            line = line.strip().split('\t')
            sam_scaffold.append(line[2])
            if line[2] in minced_scaffold and check_dna(line[9]):
                scaffold2reads[line[2]].append(line[9])

for scaffold in minced_scaffold:
    if scaffold not in sam_scaffold:
        print(scaffold, file=no_mapped_reads)

sam.close()
print("Step 6: scaffolds have no mapped have been listed in '*.without.mapping.txt'.", flush=True, file=log)
print('\n', end='', flush=True, file=log)

# Create a file to save the summary information of all spacers.
summary_txt = open("{0}.crispr.spacers.summary.txt".format(args.fasta), 'w+')
print("Scaffold_ID" + "\t" + "Repeat_Sequence" + '\t' + "Total_Spacers_Count" + '\t' + 'Unique_Spacers_Count', file=
      summary_txt)


# This finds all spacers covered by two repeats by mapping half of the direct repeat sequence,
# and put all spacers into the list of each scaffold in the dict of scaffold2spacers.
# for scaffold in scaffold2repeat.keys():
spacers_fasta = open("{0}.crispr.spacers.fasta".format(args.fasta), 'w+')
spacers_txt = open("{0}.crispr.spacers.txt".format(args.fasta), 'w+')

single_array_scaffold = list(set(minced_scaffold))
single_array_scaffold.sort()

for item in multiple_dr:
    single_array_scaffold.remove(item)


# Those spacers on the scaffolds should also be included, in case some of them not on any of the mapped reads,
# this could be possible if the repeat and spacer are very long.
for scaffold in single_array_scaffold:
    scaffold2reads[scaffold].append(fasta_dict[scaffold])

    if len(scaffold2repeat[scaffold]) % 2 == 1:
        length = len(scaffold2repeat[scaffold]) + 1
    else:
        length = len(scaffold2repeat[scaffold])

    repeat_first_half = scaffold2repeat[scaffold][:int(length/2)]
    repeat_second_half = scaffold2repeat[scaffold][int(length/2):]

    for sequence in scaffold2reads[scaffold]:
        if scaffold2repeat[scaffold] in sequence:
            for spacer in re.findall(repeat_second_half + "(.+?)" + repeat_first_half, sequence):
                if 51 > len(spacer) > 15 and repeat_first_half not in spacer and repeat_second_half not in spacer:
                    scaffold2spacers[scaffold].append(spacer)
        else:
            for spacer in re.findall(repeat_second_half + "(.+?)" + repeat_first_half, complement_reverse(sequence)):
                if 51 > len(spacer) > 15 and repeat_first_half not in spacer and repeat_second_half not in spacer:
                    scaffold2spacers[scaffold].append(spacer)

# This prints the unique spacer sequences into fasta format.
    for unique_spacer in list(set(scaffold2spacers[scaffold])):
        print(">" + scaffold + "_spacer_", end='', file=spacers_fasta)
        print(list(set(scaffold2spacers[scaffold])).index(unique_spacer), file=spacers_fasta)
        print(unique_spacer, file=spacers_fasta)

# This prints the total counts of each unique spacer in all mapped reads and scaffolds
    for unique_spacer in set(scaffold2spacers[scaffold]):
        print(scaffold + "_spacer_" + str(list(set(scaffold2spacers[scaffold])).index(unique_spacer)) + "\t" +
              unique_spacer + "\t" + str(len(unique_spacer)) + "\t" +
              str(scaffold2spacers[scaffold].count(unique_spacer)), file=spacers_txt)

    print(scaffold + '\t' + scaffold2repeat[scaffold] + '\t' + str(len(scaffold2spacers[scaffold])) + '\t' +
          str(len(set(scaffold2spacers[scaffold]))), file=summary_txt)

spacers_fasta.close()
spacers_txt.close()
summary_txt.close()

print("Step 7: all spacers from scaffolds and paired-reads have been extracted.", flush=True, file=log)
print('\n', end='', flush=True, file=log)

# Get the scaffolds without cas protein(s) and/or repeat
no_cas_scaf = open("{0}.crispr.no.cas.protein.no.repeat.fasta".format(args.fasta), 'w')
for scaffold in fasta_dict.keys():
    if scaffold not in set(minced_scaffold + cas):
        print(">" + scaffold, file=no_cas_scaf)
        print(fasta_dict[scaffold], file=no_cas_scaf)

print("Step 8: scaffolds without cas protein and/or repeat have been extracted for spacer targeting search.",
      flush=True, file=log)
print('\n', end='', flush=True, file=log)

no_cas_scaf.close()


# Build the blast database for spacer target search
os.system("formatdb -i {0}.crispr.spacers.fasta -p f".format(args.fasta))

# with open(args.query, 'r') as q:
if args.query:
    os.system("blastn -task blastn-short -db {0}.crispr.spacers.fasta -query {1} -out {1}.blastn.spacers.fasta "
              "-evalue 1e-3 -max_target_seqs 50 -num_threads 6 -perc_identity 90 -outfmt 6".
              format(args.fasta, args.query))
else:
    os.system("blastn -task blastn-short -db {0}.crispr.spacers.fasta "
              "-query {0}.crispr.no.cas.protein.no.repeat.fasta "
              "-out {0}.crispr.no.cas.protein.no.repeat.fasta.blastn.spacers.out "
              "-evalue 1e-3 -max_target_seqs 50 -num_threads 6 -perc_identity 90 -outfmt 6".format(args.fasta))

print("Step 9: blastn search for spacer target is finished.", flush=True, file=log)
print('\n', end='', flush=True, file=log)

# Parse blastn results
spacer_target = open("{0}.crispr.spacers.targets.parsed.txt".format(args.fasta), 'w')
print("Targeted_scaffold" + '\t' + "Scaffold_length" + '\t' + 'Targeted_gene' + '\t' + "Targeting_spacer" +
      '\t' + "Spacer_sequence" + '\t' + "Spacer_length" + '\t' + '\t' + "Aligned sequence" + '\t' +
      "Aligned_length" + '\t' + "Aligned_mismatches", file=spacer_target)

# gene and location dictionary
gene_id = []
gene_loc = []

if args.query:
    with open(args.gene, 'r') as g:
        for line in g:
            if line.startswith('>'):
                line = line.strip().split(' # ')
                gene_id.append(line[0][1:])
                gene_loc.append([line[1], line[2]])
else:
    with open(args.protein, 'r') as p:
        for line in p:
            if line.startswith('>'):
                line = line.strip().split(' # ')
                gene_id.append(line[0][1:])
                gene_loc.append([line[1], line[2]])

gene2loc = dict(zip(gene_id, gene_loc))


def get_gene_id(scaffold_id, align_start, align_end):
    for gene in gene_id:
        if '_'.join(gene.split('_')[0:-1]) == scaffold_id \
                and int(gene2loc[gene][0]) <= int(align_start) and int(gene2loc[gene][1]) >= int(align_end):
            return gene
        else:
            continue
    return "gene-intergenic/intergenic"


# spacer dictionary
spacer_id = []
spacer_seq = []

spacers_fasta = open("{0}.crispr.spacers.fasta".format(args.fasta), 'r')
for line in spacers_fasta:
    line = line.strip().split()
    if line[0].startswith('>'):
        spacer_id.append(line[0][1:])
    if line[0].startswith('>') is False:
        spacer_seq.append(line[0])

spacer_dict = dict(zip(spacer_id, spacer_seq))

# query scaffold dictionary if provided
query_dict = {}

if args.query:
    with open(args.query, 'r') as q:
        for record in SeqIO.parse(q, "fasta"):
            header = str(record.id).strip()
            seq = str(record.seq)
            query_dict[header] = seq


# Remove those targeted scaffold with repeat the same as that of the targeting spacer,
# only keep those spacers aligned at least 90% of length.
if args.query:
    blastn = open("{0}.blastn.spacers.fasta".format(args.query), 'r')
    for line in blastn:
        line = line.strip().split()
        print(line)
        if scaffold2repeat[line[1].split('_spacer')[0]] not in query_dict[line[0]] and \
                int(line[3])/len(spacer_dict[line[1]]) >= 0.9:
            print(line[0], '\t', len(query_dict[line[0]]), '\t', get_gene_id(line[0], line[6], line[7]),
                  '\t', line[1], '\t', spacer_dict[line[1]], '\t', len(spacer_dict[line[1]]), '\t',
                  query_dict[line[0]][int(line[6])-1:int(line[7])], '\t', line[3], '\t',
                  round(float(line[3])*(100 - float(line[2]))/100), file=spacer_target)
else:
    blastn = open("{0}.crispr.no.cas.protein.no.repeat.fasta.blastn.spacers.out".format(args.fasta), 'r')
    for line in blastn:
        line = line.strip().split()
        if scaffold2repeat[line[1].split('_spacer')[0]] not in fasta_dict[line[0]] and \
                int(line[3])/len(spacer_dict[line[1]]) >= 0.9:
            print(line[0], '\t', len(fasta_dict[line[0]]), '\t', get_gene_id(line[0], line[6], line[7]),
                  '\t', line[1], '\t', spacer_dict[line[1]], '\t', len(spacer_dict[line[1]]), '\t',
                  fasta_dict[line[0]][int(line[6])-1:int(line[7])], '\t', line[3], '\t',
                  round(float(line[3])*(100-float(line[2]))/100), file=spacer_target)

print("Step 10: blastn table parsing is finished", flush=True, file=log)
print('\n', end='', flush=True, file=log)

if args.query:
    os.system("mv {0}.blastn.spacers.fasta {1}_crispr".format(args.query, args.fasta))
else:
    pass

os.system("rm formatdb.log {0}.crispr.spacers.fasta.n*".format(args.fasta))
os.system("mv {0}.crispr.* *.Cas1.hmm.txt {0}_crispr".format(args.fasta))

print("All done! Please check files in {0}_crispr".format(args.fasta), flush=True, file=log)

