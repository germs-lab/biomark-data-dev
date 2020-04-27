
#Usage: python biomark_count.py [input file] [# of samples] [# of genes] [# of replicates per gene] > [name of output file]
#Make sure the same primers are grouped within the column before running script

# Written for  274 samples x 92 Genes - Phil April 2020

#Example of input


import sys

gene_list = []
last_gene = ''
sample_list = []
counter = 0
count_list = []
gene_list = []
count_list_counter = 0
big_count_list = []
all_genes = []
all_samples = []
d_genes = {}
d_all_samples = {}

for n, line in enumerate(open(sys.argv[1])):
    if n ==  0:  #Header filae for biomark is 17 lines
        continue
    else:
        dat = line.rstrip().split(',')
        sample_name = dat[0] #modify as neeeded
        ct = str(dat[18]) #modify as needed
        gene = dat[5] #modify as needed
        all_genes.append(gene)
        all_samples.append(sample_name)

all_genes_set = set(all_genes)
all_samples_set = set(all_samples)
for i in all_samples_set:
    d_all_samples[i] = 0

for i in all_genes_set:
    d_genes[i] = d_all_samples.copy()

for n, line in enumerate(open(sys.argv[1])):
        if n ==  0:  #Header filae for biomark is 17 lines
            continue
        else:
            dat = line.rstrip().split(',')
            sample_name = dat[0] #modify as neeeded
            ct = dat[18] #modify as needed
            gene = dat[5] #modify as needed
            d_genes[gene][sample_name] = ct

sample_list = sorted(list(d_all_samples.keys()))
gene_list = sorted(list(d_genes.keys()))
# Writing File
print('\t' + '\t'.join(sample_list))
for n, x in enumerate(gene_list):
    print(gene_list[n], end = ''),
    for m, y in enumerate(sample_list):
        print('\t' + str(d_genes[x][y]), end = '')
    print('')
