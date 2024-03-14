#!/bin/bash
name="particle_state_"	# File head name
ext=".csv"		# File extension

beg=0		# Starting iteration number
fin=8000	# Final iteration number
inv=10		# Current data set iteration interval

dig=4		# Number of maximum digit initial
trgDig=5	# Number of maximum digit final
addDig=$(($trgDig-$dig))	# Digit different

# Start iterating the data
for (( i = $beg; i <= $fin; i+=inv ))
do
	# Generate the iteration number
	# num=$(($i * $inv))
	num=$i

	# Add the leading zeros
	zero=""
	lim=$(($num))
	for (( j = 0; j < $dig; j++ ))
	do
		if (($lim != 0))
		then
			lim=$(($lim / 10))
		else
			zero=$zero$(("0"))
		fi
	done

	# Get the char of number
	if (($num == 0))
	then
		data=$zero	# Check the value of data
	else
		data=$zero$num	# Check the value of data
	fi
	fileName=$name$data$ext
	trgFileName=$name
	for (( j = 0; j < $addDig; j++ ))
	do
		trgFileName=$trgFileName$(("0"))
	done
	trgFileName=$trgFileName$data$ext

	# echo $fileName
	# echo $trgFileName
	# rm $fileName
	mv $fileName $trgFileName
done
