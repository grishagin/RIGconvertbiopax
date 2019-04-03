#!/bin/bash
fun_names=(*.R)

# add the first opening statement
# overwrite the file 
printf 'export(' > '../NAMESPACE'
# append the first filename on the same line
ff=${fun_names[@]:0:1}
printf '"%s"\n' ${ff::-2} >> '../NAMESPACE'

# loop through all functions' file names
for f in ${fun_names[@]:1}
	do 
	# store output of a grep on each of those filenames
	# to determine which of them do not belong (start witn "inner")
	temp=$(printf '%s' "$f" | grep ^internal)
	# if said output is longer than 0, ignore
	if [ ${#temp} \< 1 ]
	then
		# else append to file
		printf '\t\t,"%s"\n' ${f::-2} >> '../NAMESPACE'
	fi
done

# append a closing bracket
printf ')' >> '../NAMESPACE'
