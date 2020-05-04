from Bio import SeqIO
from datetime import datetime
import os
import csv
from optparse import OptionParser
from Bio.Alphabet import generic_dna, generic_protein
start_time = datetime.now()

"""
This script is assuming that the picrust genomes are saved in the same folder as the script (as picrust_genomes.fasta)
and that all hidden markov models are saved separately in a folder named 'hmms' (also inside the same folder) and that :
conda is installed
hmm is installed: conda install -c biocore hmmer in terminal
bio is installed: conda install biopython

It also adds the results of the hmms on to the end of the ko file currently, but if you want this to be different then just
provide a .csv file that contains a list of the assembly names and that should work fine, too (it will output .txt still)
"""


picrust_seqs = 'picrust_genomes.fasta'
hmms = os.listdir(os.getcwd()+'/hmms/')
ko = 'ko_empty.csv'

try:
    os.mkdir('hmms_out')
except:
    didnt_mkdir = True

"""
for hmm in hmms:
    os.system('nhmmer hmms/'+hmm+' '+picrust_seqs+' > hmms_out/'+hmm[:-4]+'.out ')
"""

with open(ko, 'rU') as f:
    ko_data = []
    for row in csv.reader(f):
        ko_data.append(row)


hmms_out = os.listdir(os.getcwd()+'/hmms_out')
main_dir = os.getcwd()
for hmm in hmms_out:
    ko_data[0].append(hmm[:-4])
    included_genomes = []
    with open(main_dir+'/hmms_out/'+hmm, 'rU') as f:
        contents = f.read()
    row, rows = '', []
    for a in range(len(contents)-1):
        if contents[a:a+1] == '\n':
            if row == '  ------ inclusion threshold ------':
                break
            rows.append(row)
            row = ''
        else:
            row += contents[a]
    after_start = False
    other_count = 0
    for r in range(len(rows)):
        if after_start:
            block = 0
            this_genome = ''
            for b in range(1, len(rows[r])):
                if rows[r][b-1] == ' ' and rows[r][b] != ' ':
                    block += 1
                if block == 4 and rows[r][b] != ' ':
                    this_genome += rows[r][b]
            if this_genome != '':
                included_genomes.append(this_genome)
        count = 0
        for a in range(len(rows[r])):
            if rows[r][a] == '-':
                count += 1
            if count > 40:
                after_start = True
                continue
    for a in range(len(included_genomes)):
        if included_genomes[a][-11:] == 'Description':
            included_genomes[a] = included_genomes[a][:-11]
    for b in range(len(ko_data)):
        if ko_data[b][0] in included_genomes or ko_data[b][0][:-8] in included_genomes:
            num = included_genomes.count(ko_data[b][0])
            ko_data[b].append(str(num))
        else:
            if b > 0:
                ko_data[b].append('0')
with open(ko[:-4]+'_new.csv', 'w') as f:
    writer = csv.writer(f)
    for row in ko_data:
        writer.writerow(row)
    
print(datetime.now() - start_time)
    
