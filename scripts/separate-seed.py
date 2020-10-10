from Bio import AlignIO
from io import StringIO
from tqdm import tqdm
inSeed = "Rfam/seeds/concatenated/Rfam.seed" 
outDir = "Rfam/seeds/separated" 
with open(inSeed,encoding="ISO-8859-1") as f:
    content = None
    accession = 'NA'
    for lineno, line in tqdm(enumerate(f)):
        if line.startswith('# STOCKHOLM 1.0'):
            if content:
                content = ''.join(content)
                with open(outDir+"/{}.seed".format(accession),"w") as f:
                    f.write(content)
            content = []
        elif line.startswith('#=GF AC'):
            accession = line.split()[-1]
            
        content.append(line)

