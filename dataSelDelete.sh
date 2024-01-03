#!/bin/bash
name="particle_state_"	# File head name
ext=".csv"		# File extension

cur=20		# Current data set iteration interval
tar=40		# Target data set iteration interval
dig=5		# Number of maximum digit

final=294	# Total data number

for (( i = 0; i <= $final; i++ ))
do
	data=$(($i * $cur))
	# Delete the file 
	if (( $data % $tar != 0 ))
	then
		zero=""
		lim=$(($data))
		for (( j = 0; j < $dig; j++ ))
		do
			if (($lim != 0))
			then
				lim=$(($lim / 10))
			else
				zero=$zero$(("0"))
			fi
		done
		#echo $zero$data	# Check the value of data
		fileName=$name$zero$data$ext
		rm $fileName
	fi
done
