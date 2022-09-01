# explicit-human

**Identify TF regulators for human gene modules.**

The EXPLICIT-Human package is developed based on the [EXPLICIT](https://github.com/MaShisongLab/explicit) package we published previously ([Geng *et al.* 2021](https://github.com/MaShisongLab/explicit-human#References)). It first constructs a gene expression predictor for human, which uses the expression of 1,613 TF genes to predict the expression of 23,346 non-TF genes. For every TF gene, it also predicts the gene's expression value using the expression values of all other 1,612 TF genes. The model then infers potential TF regulators for both non-TF and TF genes, as well as for gene modules participating in diverse human pathways. The EXPLICIT-Human package is a tool for users to infer potential TF regulators for their own gene modules of interest. For more information, please refer to [Wang et al, 2022](https://github.com/MaShisongLab/explicit-human#References). 

## Table of Contents
- [Install](https://github.com/MaShisongLab/explicit-human#Install)
- [Usage](https://github.com/MaShisongLab/explicit-human#Usage)
   - Identify TF regulators for human gene modules
   - Visualize TF-target gene interactions using Chord diagrams
   - Build a custom gene expression predictor
- [References](https://github.com/MaShisongLab/explicit-human#References)

## Install
[Perl](https://www.activestate.com/products/perl/downloads/), [R](https://www.r-project.org/), and the [circlize](https://www.rdocumentation.org/packages/circlize/) package in R are required for the package. [MATLAB](https://www.mathworks.com/products/matlab.html) is required only for building custom gene expression predictors. Once these softwares are installed, just download the package and start using it.

## Usage

### 1. Infer TF regulators for human gene modules

#### (1). Prepare the module file
The file `modules_to_analyze.txt` is used to store gene module information. It has two tab-separated columns for gene ids and module names, respectively. Only human Ensembl gene ids are supported. The `modules_to_analyze.txt` file came with the package can be used as a test case, which is preloaded with 1,416 human gene co-expression modules as described in [Wang *et al*, 2022](https://github.com/MaShisongLab/explicit-human#Reference).  
```shell
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
#### (2). Identify potential TF regulators for the modules
The Perl script `getHumanRegulatorTFs.pl` takes the modules from `modules_to_analyze.txt`, performs enrichment assays to identify potential TF regulators, and saves the results to a file `results.regulator.tfs.txt`. The file can be opened and previewed in EXCEL, which lists the potential TF regulators for the input modules. 
```shell
perl getHumanRegulatorTFs.pl
``` 

### 2. Visualize TF-target gene interactions using Chord diagrams

The `getChordDiagram` function in R draws Chord diagrams to visualize the interactions between TF regulators and their target genes. The function extracts TF-target gene pairs from `results.regulator.tfs.txt` for a input module and draws a Chord Diagram accordingly. <br><br>
`getChordDiagram(module, ratio, tfnum, targetnum)` <br/>
`module` - module name <br>`ratio` - relative size of the target gene area <br>`tfnum` - maximum number of TF genes to be included in the diagram <br>`targetnum` - maximum number of target genes to be included in the diagram<br>

Open an R console, change the working directory to the home directory of the EXPLICIT-Human package, and type in the following commands:
```R
# R code
# Load the scripts that define the getChordDiagram function.
source("Rscripts.R")  

# The 'circlize' package is required.
library("circlize")

# Draw a Chord diagram for Module0037.
getChordDiagram( module="Module0037", ratio = 1, tfnum = 50, targetnum = 15)

# Change the relative size of the target gene area.
getChordDiagram( module="Module0037", ratio = 0.5, tfnum = 50, targetnum = 15)

# Draw Chord diagrams for other modules.
getChordDiagram( module="Module0018", ratio = 1, tfnum = 50, targetnum = 15)
getChordDiagram( module="Module0116", ratio = 1, tfnum = 25, targetnum = 20)
getChordDiagram( module="Module0038", ratio = 1, tfnum = 20, targetnum = 15)
```

### 3. Build a custom gene expression predictor 

We used the procedure below to build the EXPLICIT-Human predictor model. One can also construct custom gene expression predictors using custom gene expression matrices accordingly. 

We used a human gene expression matrix extracted from the ARCHS4 database ([Lachmann *et al*, 2018](https://github.com/MaShisongLab/explicit-human#Reference)) to build the EXPLICIT-Human model. A human gene expression file `human_transcript_v7.h5` was downloaded from [the ARCHS4 database](https://maayanlab.cloud/archs4/download.html) and processed into CPM gene expression values. The expression values are then log-transformed via log<sub>2</sub>(CPM + 1). After quality control, 24,959 human genes in 59,097 bulk RNA-seq samples are selected to generate a human gene expression matrix, contained within a file `human_expression_extracted_from_archs4_v7.h5`, which is available via [Figshare](https://figshare.com/s/ec58e5b149c3060e1a6f). 

```shell
human_expression_extracted_from_archs4_v7.h5
  ├─ expression  					(59,097 samples [row] X 24,959 genes [column])
  ├─ gene_name						(24,959 genes)
  ├─ sample_geo_accession			(59,097 samples)
  ├─ idx_tf_gene					(specify TF genes used for model construction)
  └─ idx_non_tf_gene 				(specify non-TF genes used for model construction)
```

Download the file `human_expression_extracted_from_archs4_v7.h5` from [Figshare](https://figshare.com/s/ec58e5b149c3060e1a6f) and place it in the home directory of the EXPLICIT-Human package. Build the EXPLICIT-Human model using the following commands within MATLAB.

#### (1). Build a predictor model for non-TF genes
The following steps build a gene expression predictor for non-TF genes. <i> The analysis requires ~ 150G computer memory.</i>
```matlab
% MATLAB code
% Navigate to and start within the home directory of the EXPLICIT-Human package.

%%%%%%%%%%%%%
% Step A. Create the predictor

% Obtain the expression matrix and gene names.
mtx = h5read('human_expression_extracted_from_archs4_v7.h5','/expression');
gene_name = h5read('human_expression_extracted_from_archs4_v7.h5','/gene_name');

% idx_tf specifies which of the 24,959 genes are TFs to be used.
% idx_non_tf specifies which of the 24,959 genes are non-TFs to be used.
idx_tf = h5read('human_expression_extracted_from_archs4_v7.h5','/idx_tf_gene') == 1;
idx_non_tf = h5read('human_expression_extracted_from_archs4_v7.h5','/idx_non_tf_gene') == 1;

% Obtain the TF expression matrix and non_TF expression matrix
tf_mtx = mtx(:,idx_tf);
non_tf_mtx = mtx(:,idx_non_tf);

% Obtain the TF gene names and non_TF gene names
tf_name = gene_name(idx_tf);
non_tf_name = gene_name(idx_non_tf);

% Build the predictor model
mdl = explicit(tf_mtx, non_tf_mtx, tf_name, non_tf_name);

% Inspect the model
mdl

% The first 5 significant TF-target gene pairs
mdl.SigEdges(1:5,:)

% The SigEdges (significant interaction TF-target gene pairs) with pValue <= 1e-12 
% are extracted and saved to a file named "Hs.SigEdges.tf.to.nontf.1e-12.txt".
i = mdl.SigEdges{:,4} <= 1e-12;
writetable(mdl.SigEdges(i,:),'Hs.SigEdges.tf.to.nontf.1e-12.txt','Delimiter','tab')
```

The outputs are:
```shell
mdl

mdl = 

  explicit with properties:

                          beta: [1614x23346 double]
                   beta_pvalue: [1614x23346 double]
                       TF_name: [1x1614 string]
                   Target_name: [1x23346 string]
                         NRMSE: 0.1165
         Correlation_by_sample: [59097x1 double]
    Correlation_by_target_gene: [23346x1 double]
                           SST: [23346x1 double]
                           SSR: [23346x1 double]
                           SSE: [23346x1 double]
                         Fstat: [23346x1 double]
                       Fpvalue: [23346x1 double]
                      SigEdges: [7823519x4 table]

mdl.SigEdges(1:5,:)

ans =

  5x4 table

          Gene                  TF             beta      beta_pvalue
    _________________    _________________    _______    ___________

    "ENSG00000211752"    "ENSG00000221923"     0.0131      4.25e-07 
    "ENSG00000211752"    "ENSG00000214189"    -0.0135     1.986e-06 
    "ENSG00000211752"    "ENSG00000164916"    -0.0208     1.071e-07 
    "ENSG00000211752"    "ENSG00000152784"     0.0109     1.206e-09 
    "ENSG00000211752"    "ENSG00000171425"     0.0352     1.287e-18 
```
Next, we test how the number of training samples used for model construction affects the model's predicting performance.
```matlab
%%%%%%%%%%%%
% Step B. Investigate how the number of training samples affects the predictor model
mdl_eosn = explicit_eosn( tf_mtx, non_tf_mtx, 3000)
mdl_eosn
mdl_eosn.stat
```
The outputs are:
```shell
mdl_eosn

mdl_eosn = 

  explicit_eosn with properties:

    TestSampleNumber: 3000
                stat: [40x6 table]

mdl_eosn.stat

ans =

  40x6 table

    Run    TrainingSampleNum    R_training    NRMSE_training    R_test     NRMSE_test
    ___    _________________    __________    ______________    _______    __________

     1            1700           0.99965         0.019595       0.79484     0.62377  
     2            1725           0.99956         0.021752       0.82348      0.5603  
     3            1750           0.99943          0.02472       0.84881     0.49454  
     4            1775           0.99933         0.026958       0.86489     0.45907  
     5            1800           0.99921         0.029116        0.8769     0.42681  
     6            1850           0.99908         0.031462       0.89195     0.39348  
     7            1900           0.99883         0.035441       0.90475     0.36014  
     8            1950           0.99865         0.038212       0.91447     0.33681  
     9            2000           0.99844          0.04104       0.92365     0.31191  
    10            2100           0.99812         0.044893       0.93231     0.28931  
    11            2200           0.99776         0.049115       0.94091     0.26649  
    12            2300           0.99746         0.052181       0.94591     0.25271  
    13            2400           0.99706         0.056297       0.95201     0.23531  
    14            2500           0.99682         0.058445        0.9539     0.23092  
    15            3000           0.99558         0.068997       0.96518     0.19602  
    16            3500           0.99463         0.076022       0.97017     0.18018  
    17            4000           0.99374         0.082015       0.97408     0.16704  
    18            4500           0.99312         0.086087       0.97573     0.16109  
    19            5000           0.99265         0.088947       0.97722     0.15553  
    20            6000            0.9918         0.093665       0.97931     0.14812  
    21            7000           0.99103          0.09804       0.98097     0.14174  
    22            8000           0.99076         0.099361       0.98164     0.13914  
    23            9000           0.99024          0.10209       0.98238     0.13618  
    24           10000           0.99011           0.1027        0.9829     0.13414  
    25           12000           0.98952          0.10568       0.98354     0.13153  
    26           14000           0.98909          0.10775       0.98408     0.12932  
    27           16000           0.98886          0.10877       0.98443      0.1279  
    28           18000           0.98864          0.10988       0.98471     0.12668  
    29           20000           0.98832          0.11126       0.98497     0.12574  
    30           22000           0.98824          0.11169       0.98513     0.12501  
    31           24000           0.98804          0.11255       0.98521     0.12463  
    32           26000           0.98795          0.11302       0.98543     0.12373  
    33           28000           0.98786          0.11336       0.98547     0.12351  
    34           30000           0.98775          0.11393       0.98557     0.12303  
    35           35000           0.98759          0.11455       0.98575     0.12236  
    36           40000            0.9875          0.11497       0.98588     0.12178  
    37           45000           0.98736          0.11562         0.986     0.12129  
    38           50000           0.98725          0.11605       0.98603     0.12113  
    39           55000           0.98718          0.11635       0.98612     0.12075  
    40           56097           0.98716          0.11645       0.98613      0.1207  
```
Next, we preform a Leave One Out Cross-Validation (LOOCV) test on the model. The samples from a GEO study [GSE115736](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE115736) ([Choi *et al*, 2018](https://github.com/MaShisongLab/explicit-human#References)) are held out during model construction, the model is re-trained with all other samples, and the held-out samples are then used to test the prediction performance of the re-trained model.
```matlab
%%%%%%%%%%%%
% Step C. LOOCV test 
% Perform a LOOCV test using samples from GSE115736.
% idx_test specifies which samples to be held out during model construction.
% idx_training specifies which samples to be used during model construction.
sample_ids_gse115736 =  string (['GSM3188494' ; 'GSM3188495' ; ...
'GSM3188496' ; 'GSM3188497' ; 'GSM3188498' ; 'GSM3188499' ; ... 
'GSM3188500' ; 'GSM3188501' ; 'GSM3188502' ; 'GSM3188503' ; ...
'GSM3188504' ; 'GSM3188505' ; 'GSM3188506' ; 'GSM3188507' ; ...
'GSM3188508' ; 'GSM3188509' ; 'GSM3188510' ; 'GSM3188511' ; ...
'GSM3188512' ; 'GSM3188513' ; 'GSM3188514' ; 'GSM3188515' ; ...
'GSM3188516' ; 'GSM3188517' ; 'GSM3188518' ; 'GSM3188519' ; ...
'GSM3188520' ; 'GSM3188521' ; 'GSM3188522' ; 'GSM3188523' ; ...
'GSM3188524' ; 'GSM3188525' ; 'GSM3188526' ; 'GSM3188527' ; ...
'GSM3188528' ; 'GSM3188529' ; 'GSM3188530' ; 'GSM3188531' ; ...
'GSM3188532' ; 'GSM3188533' ; 'GSM3188534' ; 'GSM3188535']) ;
sample_geo_accession = h5read('human_expression_extracted_from_archs4_v7.h5','/sample_geo_accession');
idx_test = contains(sample_geo_accession, sample_ids_gse115736);
idx_training = ~idx_test;
test_sample_ids = sample_geo_accession(idx_test);
loocv_gse115736 = explicit_cv(tf_mtx(idx_training,:),non_tf_mtx(idx_training,:),tf_mtx(idx_test,:),non_tf_mtx(idx_test,:),test_sample_ids);

% Inspect the results
loocv_gse115736
loocv_gse115736.Test_Sample_Stat
```
The outputs are:
```shell
loocv_gse115736

loocv_gse115736 = 

  explicit_cv with properties:

             Training_Sample_Num: 59059
               Actual_Target_Exp: [38x23346 double]
            Predicted_Target_Exp: [38x23346 double]
                 Test_Sample_Num: 38
                Test_Sample_Stat: [38x3 table]
           Test_Sample_NRMSE_all: 0.1540
    Test_Sample_Mean_Correlation: 0.9795

loocv_gse115736.Test_Sample_Stat

ans =

  38x3 table

       Sample       R_test     NRMSE_test
    ____________    _______    __________

    'GSM3188494'    0.98334     0.13836  
    'GSM3188495'    0.98164     0.14615  
    'GSM3188496'    0.98326     0.14061  
    'GSM3188497'    0.98839     0.11478  
    'GSM3188498'    0.98634     0.12389  
    'GSM3188499'    0.95612     0.22899  
    'GSM3188500'     0.9809     0.14743  
    'GSM3188501'    0.98285     0.14364  
    'GSM3188503'    0.96585      0.2124  
    'GSM3188504'    0.98357      0.1346  
    'GSM3188506'    0.97226     0.17746  
    'GSM3188509'    0.98733     0.11995  
    'GSM3188510'    0.93942     0.27076  
    'GSM3188511'    0.98123     0.14674  
    'GSM3188512'    0.98444     0.13571  
    'GSM3188513'    0.97995     0.15309  
    'GSM3188514'     0.9843     0.13217  
    'GSM3188515'    0.98051     0.14901  
    'GSM3188516'    0.99026     0.10565  
    'GSM3188517'    0.98781     0.11852  
    'GSM3188518'     0.9784     0.15778  
    'GSM3188519'    0.98591     0.13074  
    'GSM3188520'    0.98097     0.14517  
    'GSM3188521'    0.98153     0.14603  
    'GSM3188522'    0.97253       0.189  
    'GSM3188523'    0.99031     0.10485  
    'GSM3188524'    0.98477     0.13123  
    'GSM3188525'    0.98323     0.14212  
    'GSM3188526'    0.98058     0.15026  
    'GSM3188527'    0.98353     0.13449  
    'GSM3188528'    0.98976     0.10847  
    'GSM3188529'      0.987     0.12229  
    'GSM3188530'    0.96982     0.19323  
    'GSM3188531'    0.97955     0.15435  
    'GSM3188532'    0.97006     0.19989  
    'GSM3188533'    0.98193     0.14184  
    'GSM3188534'     0.9829     0.14039  
    'GSM3188535'    0.95963     0.21534  
```

#### (2). Build a predictor model for TF genes
This predictor model for TF genes predicts the expresson of every TF gene based on the expression of the remaining 1,612 TF genes. It further infers TF regulators for TF genes.

```matlab
% Note that this step might take a long time to run.
mdl_tfxtf = explicit_tfxtf(tf_mtx, tf_name);

% Inspect the model
mdl_tfxtf

% The SigEdges with pValue <= 1e-12 are extracted and saved to a file.
i2 = mdl_tfxtf.SigEdges{:,4} <= 1e-12;
writetable(mdl_tfxtf.SigEdges(i2,:),'Hs.SigEdges.tf.to.tf.1e-12.txt','Delimiter','tab')
```
The outputs are:
```shell
mdl_tfxtf

mdl_tfxtf = 

  explicit_tfxtf with properties:

                   beta: [1614x1613 double]
            beta_pvalue: [1614x1613 double]
                TF_name: [1613x1 string]
                  NRMSE: 0.1139
    Correlation_by_gene: [1613x1 double]
               SigEdges: [597224x4 table]
```
#### (3). Combine SigEdges from both models
The SigEdges from both models are combined to obtain the SigEdges used for TF regulator inference in Step 1.  
```matlab
% Combine the SigEdges from the two models to obtain all SigEdges for the EXPLICIT-Human model. 
% They are the SigEdges used in Step 1 to infer TF regulators for gene modules. 
% Note that the file Hs.SigEdges.1e-12.txt are split into two files: Hs.SigEdges.part1.txt and 
% Hs.SigEdges.part2.txt so that their file sizes are under 100M. These two files are 
% contained within the data folder.
writetable([mdl.SigEdges(i,:) ; mdl_tfxtf.SigEdges(i2,:)],'Hs.SigEdges.1e-12.txt','Delimiter','tab')
```

## References

Choi J, Baldwin TM, Wong M, Bolden JE, Fairfax KA, Lucas EC, Cole R, Biben C, Morgan C, Ramsay KA, et al. 2018. Haemopedia RNA-seq: a database of gene expression during haematopoiesis in mice and humans. Nucleic Acids Research  47:D780-D785.

Geng H, Wang M, Gong J, Xu Y, and Ma S. 2021. An Arabidopsis expression predictor enables inference of transcriptional regulators for gene modules. Plant J  107:597-612.

Lachmann A, Torre D, Keenan AB, Jagodnik KM, Lee HJ, Wang L, Silverstein MC, and Ma'ayan A. 2018. Massive mining of publicly available RNA-seq data from human and mouse. Nature Communications  9:1366.

Wang Y, Zhang Y, Yu N, Li B, Gong J, Mei Y, Bao J, Ma S. 2022. Decoding transcriptional regulation via a human gene expression predictor. *submitted* 

