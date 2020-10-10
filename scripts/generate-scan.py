from string import Template
template = Template(open("templates/scan-command.txt").read())
threads=16
output="alignment/infernal/arabidopsis-thaliana/rfam-hits.txt"
model="Rfam/covariance-models/concatenated/Rfam.cm"
fasta="genome/arabidopsis-thaliana/fasta/genome.fasta"
log="log/scanning/arabidopsis.log"
script="scripts/generated/arabidopsis-scanning.sh"
command = template.substitute(threads=threads,output=output,model=model,fasta=fasta,log=log)
with open(script,"w") as f:
    print(command,file=f)
print("See {}".format(script))
