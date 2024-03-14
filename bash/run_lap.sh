#!/bin/bash
nameI="LA_S"	# File head name
nameF="_R0"	# File back name
counter="u-x"
final=6	# Total data number

for (( i = 1; i <= $final; i++ ))
do
	fileName=$nameI$(($i))$nameF
	chmod $counter $fileName
	./$fileName
done