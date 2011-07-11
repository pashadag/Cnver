#!/bin/bash
if [ $# -ne 3 ] 
then
	echo $0 forgiveness cutsize filename
	exit
fi


#want to keep going untill its not diff
filename=$3.supsmoo
count=0
cat $3 | grep -v "#" | awk '{print $1,$2,$3,int($5)}' | sort -k 2n > ${filename}0
d=1
while [ $d -ne 0 ]
do
	old_count=$count
	count=`expr $count + 1`
	#merge_close_intervals ${filename}$old_count $1 > ${filename}$count 2> /dev/null
	cat ${filename}$old_count | ${CNVER_FOLDER}/src/merge_intervals.py $1 | grep -v "#" > ${filename}$count
	diff ${filename}$old_count ${filename}$count > /dev/null 2> /dev/null
	d=$?
	if [ `expr $old_count % 10` -ne 0 ]
	then
		rm ${filename}$old_count
	fi
done
cat ${filename}$count | awk -v cutsize=$2 '{if (($3-$2)>=cutsize) {print $0}}' 
