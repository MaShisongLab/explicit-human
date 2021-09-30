# explicit-human

**Identify TF regulators for human gene modules.**

The EXPLICIT-Human package is developed based on the [EXPLICIT](https://github.com/MaShisongLab/explicit) package we published before ([Geng *et al.* 2021](https://github.com/MaShisongLab/explicit-human#Reference)). It has been used to construct a gene expression predictor for human, which use the expression of 1,613 TFs to predict the expression of 23,346 non-TF genes. It also infers TF regulator for gene modules participating in diverse human pathways. Here is a brief tutorial for the tool. For more details, please refer to the original [EXPLICIT](https://github.com/MaShisongLab/explicit) package. 

## Table of Contents
- [Install](https://github.com/MaShisongLab/explicit-human#Install)
- [Usage](https://github.com/MaShisongLab/explicit-human#Usage)
   - Identify TF regulators for gene modules
   - Draw chord diagrams to visualize the TF-target genes interaction
   - Create gene expression predictor using custom datasets
- [Reference](https://github.com/MaShisongLab/explicit-human#Reference)

## Install
This package requires [Perl](https://www.activestate.com/products/perl/downloads/), [R](https://www.r-project.org/), and the [circlize](https://www.rdocumentation.org/packages/circlize/) package in R. [MATLAB](https://www.mathworks.com/products/matlab.html) is optional, only required for creating custom gene expression predictors. After installing these software, just download the package to a local computer and start using it.

## Usage

### 1. Infer TF regulators for gene modules

#### (1). Prepare the module file
The file used to store gene modules information is `modules_to_analyze.txt`. The file contains two tab-separated columns for gene ids and module names, respectively. For gene ids, only standard human ensembl gene ids are supported. Multiple modules can be analyzed at the same time. The file `modules_to_analyze.txt` came with the package is preloaded with 1,416 human co-expression modules as described in Wang *et al.*, which can be used as test cases.  
```
Gene_Name   ModuleID
ENSG00000187033	Module0037
ENSG00000114349	Module0037
ENSG00000113262	Module0037
ENSG00000139053	Module0037
ENSG00000116703	Module0037
ENSG00000173976	Module0037
ENSG00000167791	Module0037
ENSG00000105392	Module0037
.........   ........
```
#### (2). Conduct enrichment assay to identify TF regulators for the modules
The Perl script `getHumanRegulatorTFs.pl` will take the modules from the file `modules_to_analyze.txt`, conduct enrichment assays to identify potential TF regulators, and save the results to a file named `results.regulator.tfs.txt`.

Open a command line window or shell terminal, navigate to the home directory of the EXPLICIT-Human package, and type in the following command:
```shell
perl getHumanRegulatorTFs.pl
```

The resulted file `results.regulator.tfs.txt` can be opened and viewed in EXCEL. It lists the potential TF regulators for every input modules. 

### 2. Draw chord diagrams to visualize the TF-target genes interaction

The `getChordDiagram` function can be used to draw chord diagrams for the modules in R. It extracts the TF-target gene pairs from the file `results.regulator.tfs.txt` for the input module, and then draws a chord Diagram accordingly. It has four input variables:<br>
`module` - the name of the module <br>`ratio` - the relative size of the target gene area occupies <br>`tfnum` - the maximum number of TF genes to be included in the diagram <br>`targetnum` - the maximum number of target genes to be included in the diagram.<br><br>

Open an R console and change the working directory to the home directory of the EXPLICIT package. Within the R console, type in the following commands:
```R
# R code
# Load the scripts that define the getChordDiagram function.
source("Rscripts.R")  

# The function requires the 'circlize' package
library("circlize")

# To draw a chord diagram for Module0037
getChordDiagram( module="Module0037", ratio = 1, tfnum = 50, targetnum = 15)

# Change the relative size of the target gene area
getChordDiagram( module="Module0037", ratio = 0.5, tfnum = 50, targetnum = 15)

# To draw chord diagrams for other modules
getChordDiagram( module="Module0018", ratio = 1, tfnum = 50, targetnum = 15)
getChordDiagram( module="Module0116", ratio = 1, tfnum = 50, targetnum = 15)
getChordDiagram( module="Module0038", ratio = 1, tfnum = 50, targetnum = 15)
```

### 3. Create custom gene expression predictor 

Please refer to the [EXPLICIT](https://github.com/MaShisongLab/explicit#3-create-custom-gene-expression-predictor) package for detailed procedures on how to create your own gene expression predictor. 

## Reference

To be added.

