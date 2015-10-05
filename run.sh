#!bin/bash

ALNS="/Users/julianzaugg/Documents/University/Phd/Projects/NES/RESULTS/JTT/NES_JTT_5_to_30_results/Alignments/"
QUERY="/Users/julianzaugg/Documents/University/Phd/Projects/NES/Data/NES_potential_fasta.txt"
OUTPUT="/Users/julianzaugg/Documents/University/Phd/Projects/NES/Results/"
BGS="/Users/julianzaugg/Documents/University/Phd/Projects/NES/Data/human_proteome_reviewed.txt"
# BGF="/Users/julianzaugg/Documents/University/Phd/Projects/NES/Results/hp_group_scores.txt"

python /Users/julianzaugg/Documents/University/Phd/Code\ Projects/pwm_sr/pwm_sr.py \
-alns $ALNS \
-bgs $BGS \
-q $QUERY \
-o $OUTPUT \
-tt .5

# -bgf $BGF \