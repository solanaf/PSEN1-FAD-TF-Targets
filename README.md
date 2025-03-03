# PSEN1-FAD-TF-Targets

## Build Conda Environment:
After creating and activating your environment, add to your environment:

conda config --add channels defaults <br/>
conda config --add channels bioconda <br/>
conda config --add channels conda-forge <br/>
conda update -n base -c defaults conda <br/>
conda install -c conda-forge xopen <br/>
conda install -c bioconda fastqc samtools bowtie2 cutadapt trim-galore subread <br/>
pip install cufflinks <br/>

## preprocessing.zsh:
1. Check my test.1sh script to see how you should set up all your runs with wget <br/>
- Make sure to add -P raw_reads to it <br/>
2. Give execution permissions to preprocessing.zsh <br/>
chmod +x path/to/file/preprocessing.zsh <br/>
3. In preprocessing.zsh, open in a text editor (or with nano in the terminal) <br/>
- Update GET_DATA as the path/to/data.zsh (technically can be a sh file, too. Terminal will read the wget commands the same) <br/>
- Update DIRECTORY as the path/to/where_you_want_your_stuff_put <br/>
- Don’t put any spaces around the = sign; make sure the characters touch <br/>
4. Execute preprocessing.zsh in the terminal (within your conda environment) <br/>
zsh path/to/preprocessing.zsh
