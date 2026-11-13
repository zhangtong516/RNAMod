#!/bin/bash
# create_config_template.sh
# Create configuration file template for m6A repeat analysis

cat > m6a_repeat_config.txt << 'EOF'
species	peak_file	signal_file	repeat_file
Human_H9_G	data/lifted_peaks/H9_G__m6A_hg38ToHg38.lifted.bed	data/bigwig/H9_G_m6A.bw	data/annotations/hg38_repeats.txt
Human_H9_I	data/lifted_peaks/H9_I__m6A_hg38ToHg38.lifted.bed	data/bigwig/H9_I_m6A.bw	data/annotations/hg38_repeats.txt
Monkey_LG_I	data/lifted_peaks/LG_I__m6A_nomLeu3ToHg38.lifted.bed	data/bigwig/LG_I_m6A.bw	data/annotations/nomLeu3_repeats.txt
Mouse_OU_G	data/lifted_peaks/OU_G__m6A_nomLeu3ToHg38.lifted.bed	data/bigwig/OU_G_m6A.bw	data/annotations/nomLeu3_repeats.txt
Mouse_OU_I	data/lifted_peaks/OU_I__m6A_nomLeu3ToHg38.lifted.bed	data/bigwig/OU_I_m6A.bw	data/annotations/nomLeu3_repeats.txt
Mouse_PM_G	data/lifted_peaks/PM_G__m6A_nasLar1ToHg38.lifted.bed	data/bigwig/PM_G_m6A.bw	data/annotations/nasLar1_repeats.txt
Mouse_PM_I	data/lifted_peaks/PM_I__m6A_nasLar1ToHg38.lifted.bed	data/bigwig/PM_I_m6A.bw	data/annotations/nasLar1_repeats.txt
EOF

echo "Configuration template created: m6a_repeat_config.txt"
echo "Please update file paths according to your data structure"