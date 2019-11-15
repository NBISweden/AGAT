#!/bin/bash

#This script
cleanIntermediateFile="yes"

if [[ $1 == "no" ]];then
	cleanIntermediateFile="no"
fi

for i in {0..50}_*;do

	if [[ -f $i ]];then
		#go through only the test files (not the correct output)
		if [[ ! $i =~ ^[[:digit:]]+_correct ]];then
			echo -e "\nTest of $i";
			testperfect="no"

			#case with comon locus needed and prokaryote!!
			nb=${i%"_test.gff"}
			if (( $nb == 28 ));then
				gff3_sp_gxf_to_gff3.pl --gff $i -o test.gff3 -c Name  &> /dev/null
			#case with prokaryote mode needed
			elif (( $nb == 8  ));then
				gff3_sp_gxf_to_gff3.pl --gff $i -o test.gff3   &> /dev/null
			#others
			else
				gff3_sp_gxf_to_gff3.pl --gff $i -o test.gff3 --merge_loci &> /dev/null
			fi

			#get the expected name of the correct output file we will have to check against
			pref=$(echo $i | cut -d'_' -f1)
			fileok=${pref}_correct_output.gff

			# Check against the correct output
			if [ ! -f $fileok ];then
				echo "We didnt find any correct output to check against for $i ( $fileok ) "
			else
				resu=$(diff test.gff3 $fileok)
				if [[ $resu != "" ]];then
					echo -e "There is differences between the correct reference output and the current output:\n$resu"
				else
					echo "test1 ok !"
					testperfect="yes"
				fi
			fi

			#echo "check against itself"
			#case with prokaryote!!
			nb=${i%"_test.gff"}
			if (( $nb == 8 || $nb == 28 ));then
				gff3_sp_gxf_to_gff3.pl --gff test.gff3 -o test2.gff3  &> /dev/null
			else
				gff3_sp_gxf_to_gff3.pl --gff test.gff3 -o test2.gff3  --merge_loci &> /dev/null
			fi
			resu=$(diff test2.gff3 $fileok)
			if [[ $resu != "" ]];then
					echo -e "There is differences between the original current output and the output of this file processed again:\n$resu"
			else
					echo "test2 ok !"
					if [[ $testperfect == "yes" ]];then
						echo "All test perfect !"
					fi
			fi

		fi
	fi
done

#Intermediate file cleaned
if [[ $cleanIntermediateFile == "yes" ]];then
	rm test.gff3
	rm test2.gff3
fi
