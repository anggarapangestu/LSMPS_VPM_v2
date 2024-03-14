#!/bin/bash
nameH="TF"	# File head name
# nameI="_DIR"	# File head name
# nameI="_NEU"	# File head name
# nameI="_COM1"	# File head name
nameI="_COM2"	# File head name

final=2	# Total data number

for (( i = 0; i <= $final; i++ ))
do
	fileName=$nameH$(($i))$nameI
	./$fileName
done