# PSEN1-FAD-TF-Targets

## Build Conda Environment:
After creating and activating your environment, add to your environment:

conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda update -n base -c defaults conda
conda install -c conda-forge xopen
conda install -c bioconda fastqc samtools bowtie2 cutadapt trim-galore subread
pip install cufflinks

## preprocessing.zsh:
1. Check my test.1sh script to see how you should set up all your runs with wget
- Make sure to add -P raw_reads to it
2. Give execution permissions to preprocessing.zsh
chmod +x path/to/file/preprocessing.zsh
3. In preprocessing.zsh, open in a text editor (or with nano in the terminal)
- Update GET_DATA as the path/to/data.zsh (technically can be a sh file, too. Terminal will read the wget commands the same)
- Update DIRECTORY as the path/to/where_you_want_your_stuff_put
- Don’t put any spaces around the = sign; make sure the characters touch
4. Execute preprocessing.zsh in the terminal (within your conda environment)
zsh path/to/preprocessing.zsh
