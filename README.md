# README #

### What is this repository for? ###

This is a fork of proteomics-scripts (https://github.com/NomuraRG/proteomics-scripts). Whis is described 
as "General scripts for proteomics data processing"

We are working to harmonize processing of IP2 data (http://manual.integratedproteomics.com/#/) and 
Proteome Discover (https://www.thermofisher.com/order/catalog/product/OPTON-30795) data. 

### How do I get set up? ###

#### Dependencies ####

Package and enviornment managment is done with conda.
Miniconda is sufficient. If you want to use the notebooks use Conda. You can also install conda
after Miniconda is installed. 

`$ conda install conda`

https://conda.io/docs/user-guide/install/index.html#regular-installation

#### Install and test ####

To install create the environment do the following.

```
$ git clone https://github.com/mdjones/proteomics-scripts.git
$ cd proteomics-scripts
$ conda create --name  nb-cpact nose sqlalchemy pandas matplotlib seaborn pip
$ source activate nb-cpact
$ nosetests -s
```
Note: Look into pyproteome


### Contribution guidelines ###

* Writing tests
* Code review
* Other guidelines

### Who do I talk to? ###


# proteomics-scripts
General scripts for proteomics data processing
