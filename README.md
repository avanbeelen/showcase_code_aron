# showcase_code_aron
This repository has two directories and one .doxc document:
1. 'external_project_pipeline/';
2. 'variant_calling/';
3. 'BPEXA eindverslag_groep2.docx'.

The first directory is a third-year group project where my project group and I have improved an existing pipeline in terms of: modularization, parallelization and documentation.The pipeline finds genes, in plant pathogenic fungi, that translates to 'effector' proteins; small proteins that negatively influences the immune system of economically important crops, and generates a heatmap with genome x effector to give more insight in host specificity of these fungi.The directory already has output data, so the pipeline does not have to be executed. We had a 8,0/10 for this project

The second directory is a third-year individual project where I had to build a pipeline that called variants in the Leptopilina clavipes from Illumina paired-end reads data. I had written code that: did quality control on the paired-end reads, trimmed the paired-end reads using the sliding-window technique, and made a report with the called variants. I used and implemented bowtie2, samtools and bcftools to call the variants and paralleled it in Snakemake. I had a 8,9/10 for this project

The .doxc document is the paper that my project group and I made for the external project. Please, feel free to take a look; it has a clear contrast analysis of the old code and the new, improved code. The paper gives a better understanding what improvements have been made.
