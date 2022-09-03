#!/bin/bash
if [ -f $1 ];then
cat $1 | awk '$0~/^>/{print;next;}{gsub("[^ACGTU]","N",$0);print;}' 
else
echo "input path $1 does not exists"
fi
