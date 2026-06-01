## How to use BUSCO (for transcriptome assembly QC)
Here is the BUSCO [website](https://busco.ezlab.org/). And the latest BUSCO [paper](https://academic.oup.com/nar/article/53/D1/D516/7899526). BUSCO is a good tool for transcriptome quality control because it surveys whether common evolutionarily-conserved genes were assembled correctly. First I need to open my virtual environment and then install BUSCO. More instructions [here](https://docs.alliancecan.ca/wiki/BUSCO#Installation).
```
cd ~
module load python
source my_env/bin/activate
pip install --no-index busco==6.0.0
```
Once installed, if you want more info: 
```
busco --help
```
The Digital Alliance website also suggests freezing the python environment to create a dependencies file. This is helpful because it generates a "snapshot" of what the environment looked like when you first installed. My snapshot is probably not super relevant, since I'm using my_env, which is an environment for all my bioinformatics python tools. (So this .txt file contains stuff not relevant to BUSCO.)
```
pip freeze > ~/busco-requirements.txt
```
In order to run BUSCO, you need to download relevant datasets. Browse them [here](https://busco-data.ezlab.org/v5/data/). The best dataset for my purposes is the tetrapoda_odb12 dataset. To download: 
```
cd projects/rrg-ben/for_Sam
busco --download tetrapoda_odb12
```
BUSCO automatically makes a file called busco_downloads. In busco_downloads/lineages/ you will find your dataset.
```
To actually run a BUSCO analysis:
```

