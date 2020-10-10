from string import Template
import sys
if len(sys.argv)==1:
    print("The organism should be specified.")
    sys.exit(0)
organism=sys.argv[1]
template = Template(open("templates/make-ooc.txt").read())
db="genome/{}/fasta/genome.fasta".format(organism)
ooc="genome/{}/ooc/11.ooc".format(organism)
script="scripts/generated/blat-make-ooc/{}.sh".format(organism)
command = template.substitute(db=db,ooc=ooc)
with open(script,"w") as f:
    print(command,file=f)
print("See {}".format(script))
