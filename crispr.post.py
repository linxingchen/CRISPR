##!/home/linking/.pyenv/versions/3.5.1/bin/python

import os
import argparse
import time
from Bio import SeqIO
from collections import defaultdict


# crispr.post.py
# Version 0.1
# 01.01.19
# by linxingchen


localtime = time.asctime(time.localtime(time.time()))
print("The start time is: {0}".format(localtime), flush=True)

pwd = os.getcwd()

parser = argparse.ArgumentParser(description="This script searches provided CRISPR scaffolds against provided spacers, "
                                             "spacers could be obtained using another script called crispr.py.")
requiredNamed = parser.add_argument_group('required named arguments')
requiredNamed.add_argument("-q", "--query", type=str, help="Scaffolds for spacer targeting search.")
requiredNamed.add_argument("-g", "--gene", type=str, help="Predicted genes in DNA or AA format of query scaffolds.")
requiredNamed.add_argument("-s", "--spacer", type=str, help="CRISPR spacers obtained by crispr.py or other tools.")
parser.add_argument("-p", "--process", type=int, default=6, help="Number of processes for blastn (default=6).")
args = parser.parse_args()


spacer_dict = defaultdict(str)
with open(args.spacer, 'r') as s:
    for line in SeqIO.parse(s, "fasta"):
        header = str(line.id).strip()
        seq = str(line.seq).strip()
        spacer_dict[header] = seq

os.system("formatdb -i {0} -p f".format(args.spacer))


query_dict = defaultdict(str)
with open(args.query, 'r') as q:
    for line in SeqIO.parse(q, "fasta"):
        header = str(line.id).strip()
        seq = str(line.seq).strip()
        query_dict[header] = seq

name_query = args.query.name.strip().split('/')[-1]

processes = args.process

os.system("blastn -task blastn-short -db {0} -query {1} -out {2}.blastn.spacers.out -evalue 1e-3 -max_target_seqs 50 "
          "-num_threads {3} -perc_identity 90 -outfmt 6".format(args.spacer, args.query, name_query, processes))


# gene and location dictionary
gene2loc = defaultdict(list)

with open(args.gene, 'r') as g:
    for line in g:
        if line.startswith('>'):
            line = line.strip().split(' # ')
            gene2loc[line[0][1:]] = [line[1], line[2]]


def get_gene_id(scaffold_id, align_start, align_end):
    for gene in gene2loc.keys():
        if '_'.join(gene.split('_')[0:-1]) == scaffold_id \
                and int(gene2loc[gene][0]) <= int(align_start) and int(gene2loc[gene][1]) >= int(align_end):
            return gene
        else:
            continue
    return "gene-intergenic/intergenic"


# Parse blastn results
spacer_target = open("{0}.blastn.spacers.parsed.txt".format(name_query), 'w')
print("Targeted_scaffold" + '\t' + "Scaffold_length" + '\t' + 'Targeted_gene' + '\t' + "Targeting_spacer" +
      '\t' + "Spacer_sequence" + '\t' + "Spacer_length" + '\t' + '\t' + "Aligned sequence" + '\t' +
      "Aligned_length" + '\t' + "Aligned_mismatches", file=spacer_target)


# Remove those targeted scaffold with repeat the same as that of the targeting spacer,
# only keep those spacers aligned at least 90% of length.

blastn = open("{0}.blastn.spacers.out".format(name_query), 'r')
for line in blastn:
    line = line.strip().split()
    if int(line[3])/len(spacer_dict[line[1]]) >= 0.9:
        print(line[0], '\t', len(query_dict[line[0]]), '\t', get_gene_id(line[0], line[6], line[7]),
              '\t', line[1], '\t', spacer_dict[line[1]], '\t', len(spacer_dict[line[1]]), '\t',
              query_dict[line[0]][int(line[6])-1:int(line[7])], '\t', line[3], '\t',
              round(float(line[3])*(100 - float(line[2]))/100), file=spacer_target)


os.system("rm formatdb.log {0}.n*".format(args.spacer))

print("All done! Please check {0}.blastn.spacers.parsed.txt".format(name_query))
