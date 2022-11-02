# scBionfoCourse - Week 2

**import_data_R.R** contains ways to create *Seurat* object.


## Useful Sources 

1. [Seurat Vignette](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html)

2. [Conda Cheatsheet](https://docs.conda.io/projects/conda/en/4.6.0/_downloads/52a95608c49671267e40c689e0bc00ca/conda-cheatsheet.pdf)


### Setup Ubuntu Command Line (optional for windows users)

[Ubuntu Command Line](https://www.youtube.com/watch?v=LLlfLpvQg04)

## Setup Project Environment

1. Download *.yaml* file from the project email.

2. Open command line and, go to the directory contains *.yaml* file using:

```
cd /path/to/.yaml
```

3. Initiate conda

```
conda init
```

4. Create conda environment using the *.yaml*

```
conda env create -f environment.yml
```

5. Start the conda environment with:

```
conda activate single-cell
```

6. Install *devtools* by starting R:

```
R
install.packages("devtools")
```

7. *CellChat*  and install *CellChat* using *devtools*:

```
devtools::install_github("sqjin/CellChat")
```

8. *IMPORTANT!* test if you have all the libraries!

9. If you have major difficulties setting up your system, you can as well set up your environment using
the *environment_minimal.yml* file and install all packages by hand.

### Disclaimer

Some of the code snippets were taken or created based on the sources above.