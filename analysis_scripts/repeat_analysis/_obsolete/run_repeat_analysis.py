'''
Author: Zhang Tong
Email: zhangtong516@gmail.com
Company: GIS
Date: 2025-08-07 14:37:04
LastEditors: Zhang Tong
LastEditTime: 2025-08-07 15:01:06
'''
# Update your species_data configuration
species_data = {
    'H9_G': {
        'peaks': 'data/peaks/H9_G__m6A_peaks.narrowPeak',
        'bigwig': 'data/bigwig/H9_G__m6A_treat_pileup.bw'
    },
    'H9_I': {
        'peaks': 'data/peaks/H9_I__m6A_peaks.narrowPeak', 
        'bigwig': 'data/bigwig/H9_I__m6A_treat_pileup.bw'
    },
    'LG_I': {
        'peaks': 'data/peaks/LG_I__m6A_peaks.narrowPeak',
        'bigwig': 'data/bigwig/LG_I__m6A_treat_pileup.bw'
    },
    'OU_I': {
        'peaks': 'data/peaks/OU_I__m6A_peaks.narrowPeak', 
        'bigwig': 'data/bigwig/OU_I__m6A_treat_pileup.bw'
    },
    'OU_G': {
        'peaks': 'data/peaks/OU_G__m6A_peaks.narrowPeak', 
        'bigwig': 'data/bigwig/OU_G__m6A_treat_pileup.bw'
    },
    'PM_I': {
        'peaks': 'data/peaks/PM_I__m6A_peaks.narrowPeak', 
        'bigwig': 'data/bigwig/PM_I__m6A_treat_pileup.bw'
    },
    'PM_G': {
        'peaks': 'data/peaks/PM_G__m6A_peaks.narrowPeak', 
        'bigwig': 'data/bigwig/PM_G__m6A_treat_pileup.bw'
    }
}

repeat_annotations = {
    'H9_G': 'data/repeats/hg38_repeats.out',
    'H9_I': 'data/repeats/hg38_repeats.out',
    'LG_I': 'data/repeats/nomLeu3_repeats.out',
    'OU_G': 'data/repeats/ponAbe3_repeats.out',
    'OU_I': 'data/repeats/ponAbe3_repeats.out',
    'PM_I': 'data/repeats/nasLar1_repeats.out',
    'PM_G': 'data/repeats/nasLar1_repeats.out'
}

# Run your analysis
analyzer = M6ARepeatAnalyzer(species_data, repeat_annotations)
results = analyzer.analyze_all_species()