#!/bin/bash
for rfam in $(cat data/Rfam/info/rfam-id-st50-lt10.txt);do
make rfamid=${rfam}
done
