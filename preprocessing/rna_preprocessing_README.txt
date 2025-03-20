Build Conda environment:
After creating and activating your environment, add to your environment:

conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda update -n base -c defaults conda
conda install -c conda-forge xopen
conda install -c bioconda fastqc samtools bowtie2 cutadapt trim-galore subread
pip install cufflinks

preprocessing.zsh (pair-end reads)
1. Check my test1.sh script to see how you should set up all your runs with wget
- You can get these as a bash script with the European Nucleotide Archive, or just aggregate them manually (sorry, GEO isn't working)
2. Give execution permissions to preprocessing.zsh:
chmod +x path/to/file/rna_preprocessing.zsh (or atac_preprocessing.zsh)
3. Define $DIRECTORY and $GET_DATA (Don’t put any spaces around the = sign; make sure the characters touch)
OPTION A: Open rna_preprocessing.zsh (or atac_preprocessing.zsh) in a text editor, or with nano in the terminal
- Uncomment & update GET_DATA as the path/to/data.zsh (technically can be a sh file, too. Terminal will read the wget commands the same)
- Uncoment & update DIRECTORY as the path/to/where_you_want_your_stuff_put
OPTION B:
- define $DIRECTORY and $GET_DATA in the command line with the 'export' command prior to running rna_ preprocessing.zsh (or atac_preprocessing.zsh)
4. Execute rna_preprocessing.zsh/atac_preprocessing.zsh in the terminal (within your conda environment):
zsh path/to/preprocessing.zsh (or atac_preprocessing.zsh)