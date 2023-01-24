#!/bin/bash

# The script annotates peaks using HOMER command line tool and saves
# annotated peak files into a specified directory:

INPUTS=/home/daniele/Desktop/III_course/II_semester/Kursinis_darbas/Tbx5_analysis/Inputs/BED/*;
RESULTS=/home/daniele/Desktop/IV_course/I_semester/Kursinis_projektas/Tbx5_analysis_II/Peak_annotations;
HOMER=/home/daniele/HOMER/bin

# Retrieving sequences based on BED file peak ranges:
for	bed_file in ${INPUTS};
	do
		name=$(basename ${bed_file});
		echo "Retrieving sequences from $(basename ${bed_file}) file...";
		annotation_file=$(basename ${name} bed)"annot";
        echo "Changing file name into ${annotation_file}...";
        
        sed 's/ /_/g' result.txt
        perl ${HOMER}/annotatePeaks.pl ${bed_file} mm10 > ${RESULTS}/${annotation_file}
        #sed 's/\t/,/g' ${RESULTS}/${annotation_file} > ${RESULTS}/$(basename ${annotation_file} annot)"csv"
        #rm -rf ${RESULTS}/${annotation_file}
	done
