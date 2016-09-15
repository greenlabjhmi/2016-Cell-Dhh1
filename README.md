# 2016-Cell-Dhh1

<p align="center">
  <img src = "https://github.com/greenlabjhmi/2016-Cell-Dhh1/blob/master/Abstract.png?raw=true" />
</p>

Analysis scripts related to the processing and visualization of the data in "The DEAD-box helicase Dhh1p couples mRNA decay and translation by monitoring codon optimality." Raw data from this project is available at [the following page on the GEO Archive](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE81269)

In this project are four main files:
* Pipeline.py: A script used to download raw FASTQ files and process and convert them into WIG files using the settings outlined in the paper
* Param.in: A list of parameters parsed by Pipeline.py to determine what the settings for data downloading and processing shoudl be.
* DataGen.py: A script used to extract information from the generated WIG files into processed data tables
* Plot.R: A script used to generate the figures in the paper

There is also an additional helper file (Gene.py) which is helps define a Gene class which is capable of holding all the information required for the analysis performed.
