#! /bin/bash

cp $1 temp.loci

csplit --digit=2 --quiet --prefix=split_$1 temp.loci "////+1" "{*}"

for x in split*; do
	
	sed -i -e 's/\/\/[\s\-\*\|0-9]*//g' $x;
	sed -i -e 's/L/>L/g' $x;
	sed -i -e 's/     /\n/g' $x;
	
done