from string import Template
import sys
import os
if len(sys.argv)==1:
    print("The organism should be specified.")
    sys.exit(0)
organism=sys.argv[1]
query="bpRNA"
queryPath="known-structure/{}/{}.fa".format(query,query)
template = Template(open("templates/blat-command.txt").read())
db="genome/{}/fasta/genome.fasta".format(organism)
ooc="genome/{}/ooc/11.ooc".format(organism)
output="alignment/blat/{}/{}.psl".format(organism,query)
if not os.path.exists("alignment/blat/{}".format(organism)):
   os.mkdir("alignment/blat/{}".format(organism)) 
script="scripts/generated/blat/{}.sh".format(organism)
command = template.substitute(query=queryPath,db=db,ooc=ooc,output=output)
with open(script,"w") as f:
    print(command,file=f)
print("See {}".format(script))
