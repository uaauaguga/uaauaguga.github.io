from Bio import AlignIO
from io import StringIO
from tqdm import tqdm
with open('Rfam.seed',encoding="ISO-8859-1") as f:
    content = None
    families = {}
    accession = 'NA'
    for lineno, line in tqdm(enumerate(f)):
        if line.startswith('# STOCKHOLM 1.0'):
            if content:
                content = ''.join(content)
                rfam = AlignIO.read(StringIO(content), 'stockholm')
                families[accession] = rfam
            content = []
        elif line.startswith('#=GF AC'):
            accession = line.split()[-1]
            
        content.append(line)

import os
if not os.path.exists('fasta'):
    os.makedirs('fasta')
for accession in tqdm(families.keys()):
    with open('fasta/{}.fa'.format(accession), 'w') as f:
        for record in families[accession]:
            f.write('>{}\n'.format(record.id))
            f.write(str(record.seq).replace('-', '').replace("U","T") + '\n')
