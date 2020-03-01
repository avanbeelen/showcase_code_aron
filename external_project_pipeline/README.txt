This program is made by:                                Date: Januari 10 2020
Aron van Beelen
Timo de Graaf
Aldo Vree
Lieke Vree
Mentor: Like Fokkens (University of Amsterdam)

This is the README of the pipeline that finds effectors in genomes of fungi. The pipeline only works in a Unix based System e.g. MacOS or Ubuntu (preferably Ubuntu).

------------------------
# INSTALLATION GUIDE
------------------------

1. Manually install SignalP-4.1, SignalP-5.0 and Augustus3.3
	a. Open the 'Signalp_Tools' directory in the navigator.
	b. Open the terminal in the 'SignalP_Tools' directory and type:
		'tar -xvzf signalp-5.0b.Linux.tar.gz signalp-5.0b'
		'tar -xvzf signalp-4.1g.Linux.tar.gz signalp-4.1'
		'tar -xzf augustus.current.tar.gz'
					or (for MacOS)
		'tar -xvzf signalp-5.0b.Darwin.tar.gz signalp-5.0b'
		'tar -xvzf signalp-4.1g.Darwin.tar.gz signalp-4.1'
		'tar -xzf augustus.current.tar.gz'
	   This results in two directories: 'signalp-4.1', 'signalp-5.0b' and augustus-3.3.3.
	c. Move the three directories to a desired path.
	   For example: '/home/user/tools'.
	   Later on, this is the path that you have to refer to in the config file.
	d. An extra step for signalp-4.1 has to be taken. Move to the signalp-4.1 directory (in the specified path) and open the 'signalp' file.
		
	   In the 'GENERAL SETTINGS', change the full path to the signalp-4.1 directory (in the specified path) on your system:
	   For example: '/home/user/tools/signalp-4.1'
	e. An extra step for augustus has to be taken. Move to the augustus-3.3.3 directory (in the specified path) and open the terminal. Type 'make' and hit Enter.
	   Wait. Type sudo make install and hit Enter.

2. Run INSTALL.sh
INSTALL.sh is a bash file which installs the requirements to use the pipeline. It downloads and/or installs:
BioPython3.6, bcbio-gff, Augustus, Snakemake, R and blastn.
	a. Go back to the pipeline directory and open the terminal in that directory.	
	b. Type 'bash INSTALL.sh' in the terminal and hit Enter.
		* Type password and type 'y' if necessary.
	
3. Manually install R packages
	a. Type 'R' in the terminal and hit Enter. This opens the R terminal.
	b. Type and hit Enter:
		*install.packages("dendextend")
			Question: select CRAN mirror --> choose a desired location. (We picked UK London).
		*install.packages("gplots")
		*install.packages("extrafont")
		*install.packages("ade4")
		*if (!requireNamespace("BiocManager", quietly = TRUE))
    			install.packages("BiocManager")
		 BiocManager::install(version = "3.10")
		*BiocManager::install("ctc")
		
------------------------
# PIPELINE USAGE GUIDE
------------------------

The pipeline runs with Snakemake; see installation guide above.

First, change the parameters in the config file. Read the commentary in the config file carefully.

1. Open the terminal on your computer (Unix based system).
2. Move to the directory where the pipeline is stored.
3. Check if the input files in the genome directory are without a dot (.) in the name. 
For example: sample.fasta is allowed but test.sample.fasta gives an error because there is a dot in the name.
4. Change, if necessary, the parameters in the config file.
5. Type 'snakemake' in the terminal and hit enter.

# The config.ini file
In the config file are all the parameters of the pipeline. If you want to change a parameter, change it in the config
file and run the pipeline again. All the possible options for each parameter are given in the comments.
