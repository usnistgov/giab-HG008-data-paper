#!/usr/bin/bash

# PacBio Revio
samtools stats data/HG008-N-P_PacBio-HiFi-Revio_20240125_35x_GRCh38-GIABv3.bam > data/samtools-PB-Revio/HG008-N-P_PacBio-HiFi-Revio_20240125_35x_GRCh38-GIABv3_samtools_stats.txt
samtools stats data/HG008-T_PacBio-HiFi-Revio_20240125_116x_GRCh38-GIABv3.bam > data/samtools-PB-Revio/HG008-T_PacBio-HiFi-Revio_20240125_116x_GRCh38-GIABv3_samtools_stats.txt

# BCM Revio
samtools stats data/HG008-N-D_PacBio-HiFi-Revio_20240313_68x_GRCh38-GIABv3.bam > data/samtools-BCM-Revio/HG008-N-D_PacBio-HiFi-Revio_20240313_68x_GRCh38-GIABv3_samtools_stats.txt
samtools stats data/HG008-T_PacBio-HiFi-Revio_20240313_106x_GRCh38-GIABv3.bam > data/samtools-BCM-Revio/HG008-T_PacBio-HiFi-Revio_20240313_106x_GRCh38-GIABv3_samtools_stats.txt

# UCSC ONT UL
samtools stats data/HG008-T_GRCh38_GIABv3_ONT-UL-R10.4.1-dorado0.4.3_sup4.2.0_5mCG_5hmCG_54x_UCSC_20231031.bam > data/samtools-UCSC-ONT-ULandStd/ul/HG008-T_GRCh38_GIABv3_ONT-UL-R10.4.1-dorado0.4.3_sup4.2.0_5mCG_5hmCG_54x_UCSC_20231031_samtools_stats.txt

# UCSC ONT Std
samtools stats data/HG008-T_GRCh38_GIABv3_ONT-R10.4.1-doradov0.3.4-5mCG-5hmC-63x_UCSC_20230905.bam > data/samtools-UCSC-ONT-ULandStd/std/HG008-T_GRCh38_GIABv3_ONT-R10.4.1-doradov0.3.4-5mCG-5hmC-63x_UCSC_20230905_samtools_stats.txt

# NE ONT Std
samtools stats data/HG008-N-P_GRCh38-GIABv3_ONT-R1041-dorado_0.5.3_5mC_5hmC_41x.bam > data/samtools-NE-ONT-std/HG008-N-P_GRCh38-GIABv3_ONT-R1041-dorado_0.5.3_5mC_5hmC_41x_samtools_stats.txt
samtools stats data/HG008-N-D_GRCh38-GIABv3_ONT-R1041-dorado_0.5.3_5mC_5hmC_94x.bam > data/samtools-NE-ONT-std/HG008-N-D_GRCh38-GIABv3_ONT-R1041-dorado_0.5.3_5mC_5hmC_94x_samtools_stats.txt
