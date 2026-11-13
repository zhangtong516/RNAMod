#!/bin/bash
###
 # @Author: Zhang Tong
 # @Email: zhangtong516@gmail.com
 # @Company: GIS
 # @Date: 2025-08-06 14:24:10
 # @LastEditors: Zhang Tong
 # @LastEditTime: 2025-08-06 14:24:42
### 
# setup_repeat_annotations.sh

mkdir -p data/repeats

# Download for your species (modify as needed)
echo "Downloading repeat annotations..."

# Human
wget -O data/repeats/hg38_rmsk.txt.gz http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz
gunzip data/repeats/hg38_rmsk.txt.gz

# Mouse  
wget -O data/repeats/mm10_rmsk.txt.gz http://hgdownload.cse.ucsc.edu/goldenPath/mm10/database/rmsk.txt.gz
gunzip data/repeats/mm10_rmsk.txt.gz

# Macaque
wget -O data/repeats/rheMac10_rmsk.txt.gz http://hgdownload.cse.ucsc.edu/goldenPath/rheMac10/database/rmsk.txt.gz
gunzip data/repeats/rheMac10_rmsk.txt.gz

echo "Converting to RepeatMasker format..."

# Convert UCSC format to standard RepeatMasker .out format
python3 << 'EOF'
import pandas as pd

def convert_ucsc_to_repeatmasker(ucsc_file, output_file):
    """Convert UCSC rmsk format to RepeatMasker .out format"""
    
    # Read UCSC file (skip first line which is header)
    df = pd.read_csv(ucsc_file, sep='\t', skiprows=1, header=None)
    
    with open(output_file, 'w') as f:
        # Write RepeatMasker header
        f.write("   SW  perc perc perc  query      position in query           matching       repeat            position in  repeat\n")
        f.write("score  div. del. ins.  sequence    begin     end    (left)    repeat         class/family         begin  end (left)   ID\n") 
        f.write("\n")
        
        # Process each row
        for _, row in df.iterrows():
            sw_score = row[1] if pd.notna(row[1]) else 0
            div = row[2]/10 if pd.notna(row[2]) else 0
            del_rate = row[3]/10 if pd.notna(row[3]) else 0  
            ins = row[4]/10 if pd.notna(row[4]) else 0
            chrom = row[5]
            start = row[6] + 1  # Convert to 1-based
            end = row[7]
            left = row[8]
            strand = 'C' if row[9] == '-' else '+'
            rep_name = row[10]
            rep_class = row[11] 
            rep_family = row[12]
            rep_start = row[13] if pd.notna(row[13]) else 1
            rep_end = row[14] if pd.notna(row[14]) else 1
            rep_left = row[15] if pd.notna(row[15]) else 0
            rep_id = row[16] if pd.notna(row[16]) else 0
            
            # Format line
            line = f"{sw_score:>5.0f} {div:>4.1f} {del_rate:>4.1f} {ins:>4.1f}  {chrom:<15} {start:>8} {end:>8} ({left:>8}) {strand} {rep_name:<15} {rep_family}/{rep_class:<18} {rep_start:>6} {rep_end:>6} ({rep_left:>6}) {rep_id}\n"
            f.write(line)

# Convert files
convert_ucsc_to_repeatmasker('data/repeats/hg38_rmsk.txt', 'data/repeats/human_repeats.out')
convert_ucsc_to_repeatmasker('data/repeats/mm10_rmsk.txt', 'data/repeats/mouse_repeats.out') 
convert_ucsc_to_repeatmasker('data/repeats/rheMac10_rmsk.txt', 'data/repeats/monkey_repeats.out')

print("Conversion completed!")
print("Files created:")
print("- data/repeats/human_repeats.out")
print("- data/repeats/mouse_repeats.out")
print("- data/repeats/monkey_repeats.out")
EOF

# Clean up intermediate files
rm data/repeats/*_rmsk.txt

echo "Setup complete! RepeatMasker annotations ready for analysis."