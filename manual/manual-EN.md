# PhyloSelect User Manual

## **1. About PhyloSelect**

### **1.1 Software overview**

PhyloSelect is an integrated command-line software package for phylogenetic analysis and molecular evolutionary research. Within a unified workflow, it integrates target CDS recovery, gene-level and site-level selection analysis, environmental association analysis, and structure-level comparison. It is designed to provide a consistent, reusable, and extensible workflow for multi-level evolutionary analysis of protein-coding genes.

In addition to classical molecular evolutionary methods, PhyloSelect supports sequence-context features extracted by the large-scale genomic foundation model Evo2, including sequence scores and site entropy metrics. These features can be combined with traditional codon-model results from tools such as codeml to assist in identifying candidate genes, candidate sites, and localized evolutionary heterogeneity.

### **1.2 Main functions**

PhyloSelect consists of five functional modules, covering the main steps from target protein-coding sequence acquisition to selection signal detection, environmental association analysis, and structure-level comparison. These modules can be used together as a recommended workflow or run independently according to specific research aims.

| **Module**    | **Description**                                              |
| ------------- | ------------------------------------------------------------ |
| **GeneMiner** | Recovers target genes from raw sequencing data, including shallow sequencing data. |
| **SiteView**  | Performs site-level evolutionary analysis for a single protein-coding gene by combining Evo2 entropy scores and codeml site models to identify local variation patterns and candidate positively selected sites. |
| **Selection** | Performs gene-level selection analysis for multiple protein-coding genes and compares the M0 and free-ratio models to evaluate branch-level heterogeneity in selective pressure. |
| **EnvAssoc**  | Tests statistical associations between branch-level selection signals and environmental variables using aBSREL and PGLS. |
| **Docking**   | Compares relative docking score patterns between candidate proteins and substrates or products based on predicted protein structures and molecular docking results, providing structure-level evidence for evaluating potential differences. |

### **1.3 Typical workflow**

In common analysis scenarios, PhyloSelect can be used as follows:

**GeneMiner → SiteView / Selection → EnvAssoc → Docking**

First, GeneMiner is used to prepare target CDS sequences. Then, SiteView or Selection is used to analyze site-level or branch-level selection signals. If environmental factors need to be further evaluated, EnvAssoc can be used for association analysis. If structure-level evidence is required, Docking can be used to perform molecular docking analysis for candidate genes.

This workflow is only a general reference. In practice, users can choose the starting module according to their available input data. For example, users with prepared CDS sequences can start directly from Selection or SiteView, while users with prepared structure files can run the Docking module directly.



---


## 2 **Installation and requirements**

This chapter describes the runtime requirements and installation methods for PhyloSelect. To simplify deployment, PhyloSelect manages major dependencies through Conda environment configuration or software packaging, so users usually do not need to manually configure most third-party tools. PhyloSelect can currently be installed either from a Conda channel or from the GitHub source code.

### **2.1 Recommended runtime environment**

PhyloSelect currently supports Linux and macOS systems. The following environment is recommended:

- Operating system: Linux or macOS
- Python: ≥ 3.10
- Memory: ≥ 8 GB recommended
- Disk space: 3–5 GB of available space recommended for the software, bundled components, and runtime cache, excluding user-generated analysis data

Some PhyloSelect modules, especially Evo2-related sequence scoring and feature extraction functions, require API access to remote model services. Therefore, internet access is required when running these modules.

On high-performance computing (HPC) platforms, compute nodes may have restricted internet access. If compute nodes do not have external network access, Evo2-related modules may not run properly. Users should check the network policy before submitting jobs or run these analyses in an environment with internet access.

### **2.2 Installation methods**

Users can choose either of the following installation methods.

#### **Option 1: Install from the Conda channel**

This option is recommended for users who want a quick deployment in an isolated environment.

```bash
conda create -n phyloselect python=3.10
conda activate phyloselect
conda install evanstone::phyloselect
```

#### **Option 2: Install from GitHub source code**

This option is suitable for users who need the latest version, want to debug the program, or plan to contribute to development. The runtime environment is built using the `environment.yaml` file provided in the repository.

```bash
git clone https://github.com/evanstone/phyloselect.git
cd phyloselect
conda env create -f environment.yaml
conda activate phyloselect
pip install .
```

### **2.3 Installation verification**

After installation, use the following commands to check whether the program can be called and whether the runtime environment is properly configured:

```bash
# 1. Check whether PhyloSelect is installed correctly
phyloselect --version

# 2. Check major dependencies and related service status
phyloselect check
```

**Note:** If the `check` command returns `Ready`, the current runtime environment is configured and ready for analysis. Otherwise, follow the prompts to check whether the Conda environment is activated, whether dependencies are installed, and whether the required network services are accessible.



---


## 3 **Input file format**

This chapter describes the input data requirements and basic usage rules for different PhyloSelect modules. Because each module is designed for a different analytical task, the required input data types and file organization differ across modules. Users should carefully check whether the input files meet the required format before running each module.

### **3.1 Overview of module input types**

Different PhyloSelect modules correspond to different analysis stages and therefore require different input files. The table below summarizes the main input data and common file formats for each module. More detailed file naming rules, field requirements, and parameter descriptions are provided in the corresponding module sections.

| **Module**    | **Main input data**                                          | **Common file formats** |
|---|---|---|
| **GeneMiner** | Sequencing data; target gene references; sample list         | FASTQ, FASTA, TSV       |
| **Selection** | Multiple protein-coding sequences; TreeMap; phylogenetic trees | FASTA, CSV, Newick      |
| **SiteView**  | Single protein-coding sequence; phylogenetic tree            | FASTA, Newick           |
| **EnvAssoc**  | Codon alignment; phylogenetic tree; environmental table      | FASTA, Newick, CSV      |
| **Docking**   | Protein structures; docking configuration; phylogenetic tree | PDB/CIF, CSV, Newick    |

### **3.2 Input requirements for different modules**

#### **3.2.1 GeneMiner module**

The GeneMiner module recovers target gene sequences from sequencing data under the guidance of homologous reference sequences. In PhyloSelect, this module is mainly used to obtain target protein-coding sequences (CDS) and provide input data for downstream modules such as Selection and SiteView. By default, GeneMiner sequentially performs read filtering, target gene assembly, flank trimming, and result merging to generate a target CDS dataset.

The GeneMiner module requires a sample list file in TSV format and target gene reference sequences in FASTA format. The sample list format is `<species_name><Tab><data_file_1>` for single-end data or `<species_name><Tab><data_file_1><Tab><data_file_2>` for paired-end data. Each row represents one sample. Absolute paths are recommended for data files.

For paired-end sequencing data, the sample list can be configured as follows:

```
Bupleurum_chinense	/home/user/GeneMiner2/DEMO/DEMO3/DATA/PLANT/Bupleurum_chinense_1.fq.gz	/home/user/GeneMiner2/DEMO/DEMO3/DATA/PLANT/Bupleurum_chinense_2.fq.gz
Bupleurum_fruticosum	/home/user/GeneMiner2/DEMO/DEMO3/DATA/PLANT/Bupleurum_fruticosum_1.fq.gz	/home/user/GeneMiner2/DEMO/DEMO3/DATA/PLANT/Bupleurum_fruticosum_2.fq.gz
Bupleurum_krylovianum	/home/user/GeneMiner2/DEMO/DEMO3/DATA/PLANT/Bupleurum_krylovianum_1.fq.gz	/home/user/GeneMiner2/DEMO/DEMO3/DATA/PLANT/Bupleurum_krylovianum_2.fq.gz
Bupleurum_malconense	/home/user/GeneMiner2/DEMO/DEMO3/DATA/PLANT/Bupleurum_malconense_1.fq.gz	/home/user/GeneMiner2/DEMO/DEMO3/DATA/PLANT/Bupleurum_malconense_2.fq.gz
Bupleurum_wenchuanense	/home/user/GeneMiner2/DEMO/DEMO3/DATA/PLANT/Bupleurum_wenchuanense_1.fq.gz	/home/user/GeneMiner2/DEMO/DEMO3/DATA/PLANT/Bupleurum_wenchuanense_2.fq.gz
Bupleurum_yunnanense	/home/user/GeneMiner2/DEMO/DEMO3/DATA/PLANT/Bupleurum_yunnanense_1.fq.gz	/home/user/GeneMiner2/DEMO/DEMO3/DATA/PLANT/Bupleurum_yunnanense_2.fq.gz
```

Reference sequences should be placed in an independent folder. Each gene can have one or more reference sequences. For each gene, all corresponding reference sequences should be stored in a file named `<gene_name>.fasta`. For example, if the target genes are `matK` and `psbA`, create `matK.fasta` and `psbA.fasta` in an empty reference folder and store the corresponding reference sequences in these files.

Assume that the sample list is saved as `/home/user/samples.tsv`, the Angiosperms353 reference gene set is saved in `/home/user/Angiosperms353`, and the expected output folder is `/home/user/output`. The command for running GeneMiner with default settings is:

```bash
phyloselect geneminer -f /home/user/samples.tsv -r /home/user/Angiosperms353 -o /home/user/output
```

#### **3.2.2 SiteView module**

The SiteView module is used for site-level analysis of a single gene. The input is a CDS sequence file. **Input sequences should be standard coding sequences, their lengths should be multiples of 3, and they should not contain premature stop codons**. Users can also provide a phylogenetic tree file in Newick format. If no tree file is provided, a tree can be inferred automatically using `--infer-tree`.

A typical example command is:

```bash
phyloselect siteview -s DEMO/DEMO1/DATA/rbcL.fasta -t DEMO/DEMO1/DATA/rbcL.treefile -o DEMO/DEMO1/results
```

#### **3.2.3 Selection module**

The Selection module is used to perform branch-level selection pressure analysis across multiple genes. The basic input is a directory containing multiple CDS sequence files, where each sequence file corresponds to one independent gene. The program treats the sequence files in the directory as analysis targets and performs batch branch-level evolutionary analysis. File names in the input directory are used as task identifiers and result directory names.

For example, the input directory `DEMO/DEMO2/DATA/FASTA` can be organized as follows:

```text
.
├── 4HPAAS.fasta
├── 4HPAR1.fasta
└── 4HPAR2.fasta
```
By default, the Selection module performs multiple sequence alignment, phylogenetic tree construction, and branch-level selection pressure analysis for each input gene. This mode is suitable when users only provide target CDS sequences and have not prepared phylogenetic trees in advance.

If users already have preconstructed phylogenetic trees, these tree files can also be provided directly. In this case, a `tree_file_map` file must also be provided to specify the corresponding phylogenetic tree for each sequence file. The program will use the specified tree topology for each gene according to this mapping relationship. If a gene in `tree_file_map` does not have a specified tree file, or if the specified tree file does not exist, the program will automatically reconstruct a phylogenetic tree for that gene.

The `tree_file_map` file uses CSV format and should contain at least the following two columns:

- `sequence`: path to the input sequence file
- `tree`: path to the corresponding phylogenetic tree file

To ensure correct analysis, species names in the sequence files should match the labels in the corresponding tree files. The `sequence` and `tree` paths recorded in the `tree_file_map` file should also strictly match the actual file locations.

For example, `DEMO/DEMO2/DATA/tree_file_map.csv` can contain:

```text
sequence,tree
DEMO/DEMO2/data/4HPAAS.fasta,DEMO/DEMO2/tree_file/4HPAAS.treefile
DEMO/DEMO2/data/4HPAR1.fasta,DEMO/DEMO2/tree_file/4HPAR1.treefile
DEMO/DEMO2/data/4HPAR2.fasta,DEMO/DEMO2/tree_file/4HPAR2.treefile
```

Assume that the target CDS sequence directory is `home/user/DEMO/DEMO2/DATA/FASTA`, the tree mapping file is `/home/user/DEMO/DEMO2/DATA/tree_file_map.csv`, and the output directory is `/home/user/output`. The command for running the Selection module with user-provided trees is:

```
phyloselect selection -i home/user/DEMO/DEMO2/DATA/FASTA -o /home/user/output
```

#### **3.2.4 EnvAssoc module**

The EnvAssoc module is used to analyze the relationship between branch-level selection signals and environmental factors within a phylogenetic framework. This module takes a codon alignment file, a phylogenetic tree file, and an environmental data table as input. It first calls aBSREL to identify branch-level selection signals and then uses PGLS to test statistical associations between branch-level selection metrics and environmental variables.

The basic input includes one codon alignment file, one phylogenetic tree file corresponding to that alignment, and one environmental data table. The codon alignment file supports FASTA or PHYLIP format. The phylogenetic tree file should be in Newick format, and terminal labels in the tree should match the sequence names in the alignment file. The environmental data table uses CSV format and must contain a column named `taxon`, which is used to match the tree terminal labels and sequence names.

For example, part of the environmental data table `DEMO/DEMO3/DATA/env.csv` is shown below:

```text
taxon,elevation,bio_15,GHFD,WRB2_CODE,AWC
Rhodiola_quadrifida,1626,60.41747,21,18,0
Rhodiola_crenulata,990,78.01089,9,19,159
Rhodiola_tibetica,3909,92.29301,56,8,150
```

The `taxon` column is required, and its content should exactly match the terminal labels in the phylogenetic tree and the sequence names in the alignment file. The remaining columns are environmental variables to be tested, such as elevation, climate variables, soil factors, or other ecological factors. Before running the analysis, users should check whether the environmental variables contain missing values, nonnumeric characters, or obvious outliers, as these may affect the PGLS results.

A typical example command is:

```bash
phyloselect envassoc -s DEMO/DEMO3/DATA/4HPAAS3.fasta -t DEMO/DEMO3/DATA/4HPAAS3.treefile -e DEMO/DEMO3/DATA/env.csv -o DEMO/DEMO3/results
```

If users need to adjust the behavior of aBSREL or PGLS, additional parameters can be set, including the genetic code, branch set to be tested, multiple-hit model, whether to enable synonymous rate variation (SRV), significance threshold, and minimum number of taxa. Detailed descriptions of these parameters are provided in the parameter section below.

#### **3.2.5 Docking module**

The Docking module is used to compare relative ligand-binding patterns among different genes or among proteins from different species at the structural level. It provides reference information for interpreting potential functional differences in candidate genes or candidate sites. This module takes protein structure files and a docking configuration file as input, performs molecular docking for specified substrates and products, and organizes and displays docking results across different lineages together with a phylogenetic tree.

The basic input includes one docking configuration file and one phylogenetic tree file. The docking configuration file uses CSV format and specifies the receptor structure directory, substrate and product information, optional cofactors, and reference structure information for each gene. The phylogenetic tree is used to display docking results in a phylogenetic context.

Each row in the docking configuration file corresponds to one gene and should contain at least the following columns:

- `Gene`: gene name or identifier
- `ReceptorDir`: protein structure directory corresponding to the gene
- `Substrate`: substrate name
- `Product`: product name

The configuration file can also include the following optional columns:

- `Cofactor`: cofactor name
- `Reference`: reference structure ID or reference ligand information

The `Reference` column can be omitted. If `Reference` is not provided, the program will use an automatic pocket detection strategy to search for potential binding pockets during downstream analysis.

For example, part of the configuration file `DEMO/DEMO4/DATA/docking_config.csv` is shown below:

```text
Gene,ReceptorDir,Substrate,Product,Cofactor,Reference
UGT1,DEMO/DEMO4/DATA/receptors/UGT1,Tyrosol,Salidroside,UDP-Glucose,8ITA
UGT2,DEMO/DEMO4/DATA/receptors/UGT2,Tyrosol,Salidroside,UDP-Glucose,8ITA
UGT3,DEMO/DEMO4/DATA/receptors/UGT3,Tyrosol,Salidroside,UDP-Glucose,8ITA
...
```

Here, `ReceptorDir` points to the protein structure directory for the corresponding gene. Each directory usually contains structure files for multiple species. The file format is `.cif`, and each file corresponds to one species. For example, the `UGT1` directory can be organized as follows:

```text
UGT1/
├── rhodiola_amabilis.cif
├── rhodiola_bupleuroides.cif
├── rhodiola_chrysanthemifolia.cif
...
```
To ensure comparability of results, the species sets in different gene directories should be as consistent as possible. In addition, the species names in the structure file names should match the terminal labels in the phylogenetic tree, so that docking results can be integrated and interpreted in a phylogenetic context.

A typical example command is:

```bash
phyloselect docking -c DEMO/DEMO4/DATA/docking_config.csv -t DEMO/DEMO4/DATA/UGT.tree -o DEMO/DEMO4/results
```



---

## 4. Parameters

PhyloSelect uses a modular command-line structure. Different modules have different parameter combinations. This section describes the parameters for each module and distinguishes required and optional parameters.

### **4.1 Common options**

The following options are shared by multiple modules and are used to control basic program behavior:

- `--threads`: Specify the number of threads used for analysis.
- `--force`: Overwrite existing results if the output directory already exists.
- `--keep-intermediate`: Keep temporary and intermediate files generated during the run for debugging or inspection.
- `--skip-plot`: Skip figure generation and output only tables or model-based results.

### 4.2 GeneMiner

Our team has developed and released the latest version of GeneMiner2. In PhyloSelect, the GeneMiner module only calls its target CDS extraction function to prepare input sequences for downstream Selection and SiteView analyses. Therefore, this section only introduces the three required parameters for this function. Users who need the graphical interface or the full functions of GeneMiner2 can refer to  [official GeneMiner2 documentation](https://github.com/sculab/GeneMiner2).

- `-f FILE`: Sample list in TSV format. Each row corresponds to one sample. The format is `<species_name><Tab><data_file_1>` for single-end data and `<species_name><Tab><data_file_1><Tab><data_file_2>` for paired-end data.
- `-r DIR`: Target gene reference directory. Each target gene should correspond to one FASTA file, preferably named `<gene_name>.fasta`.
- `-o DIR`: Output directory for intermediate files and recovered target CDS results.

If no specific action is provided, the program runs the default target CDS extraction workflow. Users can also specify actions at the end of the command.

Advanced parameters such as k-mer size, number of parallel processes, trimming mode, merging strategy, multiple sequence alignment program, and tree-building program usually do not need to be modified for routine CDS extraction. For detailed parameter descriptions, see the latest [GeneMiner2 command-line manual](https://github.com/sculab/GeneMiner2/blob/master/manual/ZH_CN/command_line.md)

### 4.3 SiteView

```text
usage: phyloselect siteview [-h] -s FILE -o OUTPUT_DIR [-t FILE] [-g SPECIES [SPECIES ...]] [-n INT] [--conservation-threshold FLOAT] [--no-run-paml] [--entropy-file FILE] [-p INT] [--force] [--keep-intermediate] [--infer-tree] [--codon-aligned] [--skip-plot]

options:
  -h, --help            show this help message and exit
  -s FILE, --seq FILE   Input CDS FASTA file (use --codon-aligned if already aligned).
  -o OUTPUT_DIR, --output-dir OUTPUT_DIR
                        Output directory.
  -t FILE, --tree FILE  Phylogenetic tree in Newick format (use --infer-tree if not provided).
  -g SPECIES [SPECIES ...], --outgroups SPECIES [SPECIES ...]
                        One or more outgroup taxa used for rooting.
  -n INT, --label-max-chars INT
                        Maximum number of characters shown for each leaf label in the tree.
  --conservation-threshold FLOAT
                        Entropy threshold used to distinguish conserved and poorly conserved sites.
  --no-run-paml         Skip codeml site-model analyses and only perform Evo2-based scoring.
  --entropy-file FILE   Precomputed entropy matrix file, mainly for debugging or reproducing previous results.
  -p INT, --threads INT
                        Total number of threads used for codeml model runs and related analysis steps.
  --force               Overwrite the output directory if it already exists.
  --keep-intermediate   Keep temporary and intermediate files generated during the run.
  --infer-tree          Infer a phylogenetic tree automatically when no tree file is provided.
  --codon-aligned       Treat the input FASTA file as an already codon-aligned sequence file and skip alignment-related preprocessing
  --skip-plot           Skip figure generation and output only tabular and model-based results.
```
**Required parameters**

- `seq`: Input FASTA file for a single CDS. If the input file is already codon-aligned, use `--codon-aligned` to skip alignment-related preprocessing.
- `output-dir`: Output directory for analysis results, intermediate files, and figures.

**Optional parameters**

- `tree`: Input phylogenetic tree in Newick format. Tree labels should match the sequence names in the input file.
- `infer-tree`: Infer a phylogenetic tree automatically when no tree file is provided.
- `outgroups`: Specify one or more outgroup taxa for rooting.
- `label-max-chars`: Set the maximum number of characters shown for each leaf label in the tree.
- `codon-aligned`: Treat the input file as an already codon-aligned sequence file and skip alignment preprocessing.
- `conservation-threshold`: Site entropy threshold used to distinguish relatively conserved sites from relatively variable sites. The default value is 1.4.
- `no-run-paml`: Skip codeml site-model analysis and only perform Evo2-related scoring and visualization.
- `entropy-file`: Provide a precomputed entropy file, mainly for debugging, reproduction, or testing.


### 4.4 Selection 

```text
usage: phyloselect selection [-h] -i PATH -o DIR [-m FILE] [-p INT] [--force] [--keep-intermediate] [--skip-plot]

options:
  -h, --help            show this help message and exit
  -i PATH, --input PATH
                        Input FASTA directory or a semicolon-separated list of FASTA files.
  -o DIR, --output-dir DIR
                        Output directory.
  -m FILE, --tree-file-map FILE
                        Optional file mapping each FASTA file to its corresponding phylogenetic tree.
  -p INT, --threads INT
                        Total number of threads used for codeml model runs and related analysis steps.
  --force               Overwrite the output directory if it already exists.
  --keep-intermediate   Keep temporary and intermediate files generated during the run.
  --skip-plot           Skip figure generation and output only tabular and model-based results.
```

**Required parameters**

- `input`: Input sequence files. This parameter accepts either a directory containing multiple sequence files or multiple sequence file paths separated by semicolons. Each sequence file is treated as an independent gene.
- `output-dir`: Output directory for analysis results and intermediate files.

**Optional parameters**

- `tree-file-map`: Tree mapping file used to specify the tree file corresponding to each input sequence file. If this parameter is not provided, the program will automatically construct a phylogenetic tree for each gene. See [Section 3.2.3](#323-selection-module) for the mapping file format.

### 4.5 EnvAssoc module

```text
usage: phyloselect envassoc [-h] -s FILE -t FILE -e FILE -o DIR [--code CODE] [--branches BRANCH_SET] [--multiple-hits MODEL]
                            [--srv YES/NO] [--timeout-sec SECONDS] [--no-log-transform] [--alpha FLOAT] [--min-n INT] [-p INT]
                            [--force]

options:
  -h, --help            show this help message and exit
  -s FILE, --seq FILE   Codon-aligned sequence file used for aBSREL analysis.
  -t FILE, --tree FILE  Phylogenetic tree in Newick format matching the input alignment.
  -e FILE, --env FILE   Environmental metadata table in CSV format.
  -o DIR, --output-dir DIR
                        Output directory.
  --code CODE           Genetic code used by HyPhy aBSREL (default: Universal).
  --branches BRANCH_SET
                        Branch set to be tested in aBSREL (default: All).
  --multiple-hits MODEL
                        Multiple-hit substitution model used by aBSREL (default: None).
  --srv YES/NO          Enable synonymous rate variation in aBSREL (default: No).
  --timeout-sec SECONDS
                        Maximum runtime allowed for aBSREL.
  --no-log-transform    Disable log1p transformation of omega_weighted before PGLS analysis.
  --alpha FLOAT         Significance threshold for FDR-adjusted Q-values (default: 0.05).
  --min-n INT           Minimum number of taxa required for each environmental variable.
  -p INT, --threads INT
                        Number of threads used for analysis.
  --force               Overwrite the output directory if it already exists.
```

**Required parameters**

- `seq`: Input codon alignment for aBSREL analysis. The file can be in FASTA or PHYLIP format and should be a codon alignment suitable for branch-level analysis.
- `tree`: Input phylogenetic tree in Newick format. Tree labels should match the taxon names in the alignment file.
- `env`: Environmental data table in CSV format. This table must contain a `taxon` column, whose values should match the taxon names in the tree and alignment file. See [Section 3.2.4](#324-envassoc-module) for an input format example.
- `output-dir`: Output directory for analysis results and intermediate files.

**Optional parameters**

- `code`: Specify the genetic code used by HyPhy aBSREL. The default is `Universal`.
- `branches`: Specify the branch set to be tested in aBSREL. The default is `All`.
- `multiple-hits`: Specify the multiple-hit substitution model used by aBSREL.
- `srv`: Whether to enable synonymous rate variation (SRV).
- `timeout-sec`: Set the maximum runtime for aBSREL in seconds.
- `no-log-transform`: Do not apply log1p transformation to `omega_weighted` before PGLS analysis.
- `alpha`: Significance threshold for FDR-adjusted Q-values. The default is 0.05.
- `min-n`: Minimum number of taxa required for each environmental variable to be included in PGLS analysis. The default is 5.


### 4.6 Docking 模块

```text
usage: phyloselect docking [-h] -c FILE -t FILE -o DIR [--force] [--keep-intermediate]

options:
  -h, --help            show this help message and exit
  -c FILE, --dock-config FILE
                        CSV file mapping scripts to substrate/product CIDs.
  -t FILE, --tree FILE  Phylogenetic tree in Newick format.
  -o DIR, --output-dir DIR
                        Output directory.
  --force               Overwrite the output directory if it already exists.
  --keep-intermediate   Keep temporary and intermediate files.
```

**Required parameters**

- `docking-config`: Docking configuration file in CSV format. This file specifies gene names, receptor structure directories, substrates, products, and optional cofactors or reference structure information. See [Section 3.2.5](#325-docking-module) for the configuration file format.
- `tree`: Input phylogenetic tree in Newick format. Species names in the tree should match the species names corresponding to the receptor structure files.
- `output-dir`: Output directory for docking results and intermediate files.




---
## 5. Outputs

### **5.1 SiteView module outputs**

After the analysis is completed, SiteView results are saved in the corresponding task directory under the output directory. A typical directory structure is shown below:

```text
siteview/
├── file_input/
├── evo_output/
├── paml_output/
├── plots/
```

The meanings of the main subdirectories are as follows:

- `file_input/`: Stores input and preprocessing files used by SiteView, including the codon alignment and phylogenetic tree. If no tree is provided by the user, the automatically inferred tree is also saved in this directory.
- `evo_output/`: Stores site-level features calculated by Evo2, such as site entropy or related scoring information. These results can be used to evaluate site conservation and variation.
- `paml_output/`: Stores PAML/codeml site-model results, including parameter estimates, log-likelihood values, LRT results, and candidate positively selected sites.
- `plots/`: Stores figures generated by SiteView, including site-level conservation plots, entropy curves, or integrated evolutionary profiles.

The following files are usually the core results to inspect first:

| **File**                           | **Description**                                              |
| ---------------------------------- | ------------------------------------------------------------ |
| file_input/example.codon.aln.fasta | Codon-aligned sequence file                                  |
| file_input/example.paml.tree       | Input or automatically inferred phylogenetic tree            |
| evo_output/example_entropy.csv     | Site entropy table recording entropy and related variation metrics for each site in the alignment |
| paml_output/site_test_summary.csv  | Summary of PAML site-model results, including parameter estimates, model comparison statistics, and positively selected sites |
| plots/example_site_evolution.png   | Integrated visualization of site-level evolutionary patterns, including the phylogenetic tree, amino acid conservation heatmap, and site-wise variation profile |

Users are advised to first inspect `plots/example_site_evolution.png` to obtain an overview of the gene’s evolutionary pattern across species, and then examine `evo_output/example_entropy.csv` and `paml_output/site_test_summary.csv` for site-level scores and model-based selection results.

### **5.2 Selection module outputs**

After the analysis is completed, Selection results are saved in the corresponding task directory under the output directory. A typical directory structure is shown below:

```text
selection/
├── Evo_dNdS.png
├── gene1/
│   ├── file_input/
│   ├── paml_output/
│   └── result/
│       ├── *_NLL_score.csv
│       └── *_omega.csv
├── gene2/
└── ...
```

The meanings of the main files and subdirectories are as follows:

- `Evo_dNdS.png`: Integrated summary figure showing Evo2 scores and dN/dS patterns across all genes.
- `gene1/`, `gene2/`, etc.: Independent result directories for each gene.
  - `file_input/`: Input and preprocessing files for the corresponding gene.
  - `paml_output/`: Raw PAML/codeml output files for the corresponding gene.
  - `result/`: Final statistical results for the corresponding gene.
    - `*_NLL_score.csv`: Model score or likelihood results for the gene.
    - `*_omega.csv`: ω (dN/dS) estimates for the gene.

Users are advised to first inspect `Evo_dNdS.png` to obtain an overview of selection patterns across genes. For genes of interest, users can then enter the corresponding gene directory and focus on `*_omega.csv` and `*_NLL_score.csv` under `result/`, together with model outputs in `paml_output/`.

### **5.3 EnvAssoc module outputs**

After the analysis is completed, EnvAssoc results are saved in the corresponding task directory under the output directory. A typical directory structure is shown below:

```text
envassoc/
├── *.ABSREL.json
├── TableS1_aBSREL_full_branch_results.csv
├── Table1_aBSREL_branch_summary.csv
└── Table2_PGLS_environment_association.csv
```

The meanings of the output files are as follows:

- `*.ABSREL.json`: Raw HyPhy aBSREL output file containing detailed branch-level model results, including ω distributions and significance test statistics for each branch. This file is suitable for detailed parsing or secondary analysis.
- `TableS1_aBSREL_full_branch_results.csv`: Complete aBSREL branch result table summarizing selection signals and statistical metrics for all branches.
- `Table1_aBSREL_branch_summary.csv`: Summary table of aBSREL branch-level results, including `ω_weighted`, LRT, P-value, Q-value, and filtering status for each branch. This table provides an overview of branch-level selection signals and indicates which branches are retained for downstream environmental association analysis.
- `Table2_PGLS_environment_association.csv`: PGLS result table recording statistical associations between environmental variables and branch-level selection metrics, including regression coefficients and significance levels.

Users are advised to first inspect `Table1_aBSREL_branch_summary.csv` to identify branch-level selection signals and filtering status. Then, `Table2_PGLS_environment_association.csv` can be used to evaluate whether these signals are statistically associated with environmental variables. For detailed inspection, users can further refer to `TableS1_aBSREL_full_branch_results.csv` and the raw `*.ABSREL.json` output.

### **5.4 Docking module outputs**

After the analysis is completed, Docking results are saved in the corresponding task directory under the output directory. A typical directory structure is shown below:

```text
```text
docking/
├── ligands/
├── receptors/
├── pockets/
├── docking_results/
├── docking_summary.csv
├── LigandBindingProfile.png
├── SubstrateProductPreference.png
└── CofactorBindingProfile.png
```

The meanings of the main files and subdirectories are as follows:

- `ligands/`: Stores ligand structure files used in docking analysis, including substrates, products, and optional cofactors.
- `receptors/`: Stores receptor structure files for candidate proteins. Each gene is usually organized as an independent subdirectory containing structures from different species.
- `pockets/`: Stores binding pocket detection results. If no reference structure is provided in the configuration file, the program uses automatic pocket detection to identify potential docking regions.
- `docking_results/`: Stores raw docking outputs and intermediate files for receptor–ligand combinations.
- `docking_summary.csv`: Summarizes docking scores for candidate proteins and their corresponding ligands.
- `LigandBindingProfile.png`: Shows the overall docking score profiles between candidate proteins and specified ligands, including substrates and products.
- `SubstrateProductPreference.png`: Shows relative substrate–product docking preference patterns across candidate proteins or species.

- `CofactorBindingProfile.png`: Normalized heatmap of cofactor docking score patterns. This file is generated when cofactor information is provided in the input configuration.

**Core results**

| **File**                       | **Description**                                              |
| ------------------------------ | ------------------------------------------------------------ |
| LigandBindingProfile.png       | Overall docking score profile for candidate proteins and specified ligands |
| SubstrateProductPreference.png | Relative docking preference between substrate and product    |
| docking_summary.csv            | Summary table of docking scores for downstream comparison    |

Users are advised to first inspect `LigandBindingProfile.png` to obtain an overview of ligand-related docking score patterns, and then inspect `SubstrateProductPreference.png` to compare substrate–product preference across candidate proteins or species. After reviewing the visualization results, users can examine `docking_summary.csv` for the corresponding numerical docking scores.

If users need to trace individual docking runs or verify result sources, the raw output files can be found in `docking_results/`.



---

## **6 Running examples**

This chapter provides two types of examples for PhyloSelect: quick-start examples and manuscript case-study examples.

For the `quickstart/` data included in the repository, please enter the PhyloSelect repository root directory before running the commands to ensure that relative paths in the example data can be correctly resolved. Data required for the manuscript case-study analysis should be downloaded from Figshare in advance, and input paths should be configured according to the instructions below.

### **6.1 Quick start**

This section provides minimal running examples for the four main PhyloSelect analysis modules. These examples help users quickly check whether the software can run properly and provide an initial understanding of the basic inputs and outputs of each module. For more complete data descriptions, result interpretation, and example figures, please continue to [6.2 Manuscript case-study analysis](#62-manuscript-case-study-analysis).

Before running the quickstart examples, please prepare the example data in your working directory：

- If PhyloSelect was installed via Conda, copy the built-in `quickstart/` dataset to the current directory:

```bash
cp -r $(python -c "import phyloselect, pathlib; print(pathlib.Path(phyloselect.__file__).parent / 'quickstart')") .
```

- Alternatively, you can clone the GitHub repository to obtain the complete example data:

```bash
git clone https://github.com/scu-shiyi/PhyloSelect.git
cd PhyloSelect
```

If PhyloSelect was installed from GitHub source code, there is no need to clone the repository again. Simply `cd PhyloSelect`.

#### 6.1.1 Selection 模块

**Example Data:**

- CDS sequence of a single gene: `quickstart/sequences/gene1.fasta`

- Corresponding phylogenetic tree file: `quickstart/trees/test1.nwk`

**Run Command:**

```text
cd phyloselect
phyloselect siteview \
-s quickstart/sequences/gene1.fasta \
-t quickstart/trees/test1.nwk \
-o outputdir 
```

**Main Outputs:**

- `gene2EvolutionarySites.png` : Provides a quick overview of the gene’s evolutionary pattern across species, including phylogenetic relationships, site conservation, and variation trends.
- `SiteTestSummary.csv` : Summarizes analysis results across different site models, useful for identifying potential positively selected sites.

#### 6.1.2 SiteView 模块

**Example Data:**

- Directory of CDS sequences for multiple genes : `quickstart/sequences`

> **Note:** We do not provide the phylogenetic trees in this example. If you have your own phylogenetic trees, prepare a configuration file mapping sequences to their corresponding tree files(refer to the [user manual](manual.md)), and specify it using the `--tree-file-map` option.

**Run Command:**

```bash
phyloselect selection -i quickstart/sequences -o outputdir
```

**Main Outputs:**

- `Evo_dNdS.png`：Shows the Evo scores and dN/dS patterns of different genes across the phylogenetic tree.
- Individual folders for each gene containing:
  - `*_omega.csv` – ω (dN/dS) estimates for different branches.
  - `*_omega.csv` :  ω (dN/dS) estimates for different branches.

#### 6.1.3 EnvAssoc 模块

**Example Data:**

- CDS sequence file: `quickstart/sequences/gene3.fasta`
- Environmental trait matrix: `quickstart/config/env_traits.csv`
- Corresponding phylogenetic tree file: `quickstart/trees/test1.nwk`

**Run command:**

```bash
phyloselect envassoc \
-s quickstart/sequences/gene3.fasta \
-e quickstart/config/env_traits.csv \
-t quickstart/trees/test1.nwk \
-o outputdir
```

**Main Outputs:**

- `Table1_aBSREL_branch_summary.csv` : Summary table of aBSREL branch-level results. This table provides an overview of branch-level selection signals and indicates which branches are retained for downstream environmental association analysis.
- `Table2_PGLS_environment_association.csv` : association results between branch-level selection signals and environmental variables based on PGLS analysis.

#### 6.1.4 Docking 模块

**Example Data:**

- Docking configuration file: `quickstart/config/docking_config.csv`
- Phylogenetic tree file for result visualization: `quickstart/trees/test2.nwk`

> **Note:** The docking configuration file should include the target genes, paths to their modeled protein structures, substrates, products, key cofactors, and reference proteins for active-pocket definition. Please refer to the [user manual](manual.md) for the required file format.

**Run command:**

```bash
phyloselect docking \
-c quickstart/config/docking_config.csv \
-t quickstart/trees/test2.nwk \
-o outputdir
```

**Main Outputs：**

- `TotalBindingEnergy.png` : visualization of binding energy patterns across the tested receptors or lineages.
- `DockingResults.csv` : detailed docking results, including receptor–ligand combinations and binding energy values.
- Additional visualization files, such as `SubstrateProductPreference.png` and `CofactorBindingEnergy.png`, are also generated.

### 6.2 Manuscript case-study analysis

This section provides a complete reproduction workflow for the example analyses presented in the associated publication. These datasets are not bundled directly with the GitHub repository, but are available through Figshare. Users should first download the example dataset from Figshare at https://doi.org/10.6084/m9.figshare.32333805.v2, and then follow the instructions in this section to enter the corresponding directories and run the commands.The published example uses UGT genes from *Rhodiola* as the primary dataset to demonstrate the typical outputs of PhyloSelect in gene-level selection analysis, site-level evolutionary profiling, environmental association analysis, and structural docking comparison. 

#### 6.2.1 Gene-level selection analysis

The Selection module is used to perform gene-level selection analysis across multiple UGT genes by comparing the M0 and free-ratio models to evaluate whether different genes exhibit branch-level heterogeneity in selective pressure. In addition, this module incorporates Evo2 log-likelihood scores for sequence evaluation. These scores are generated by a large-scale genomic language model based on learned sequence-context patterns and can serve as auxiliary indicators for comparing sequence constraint and overall evolutionary characteristics among genes.

Assuming that the CDS sequences of multiple UGT genes are stored in `home/user/DEMO/DEMO1/DATA/FASTA`, run the following command:

```bash
phyloselect selection -i home/user/DEMO/DEMO1/DATA/FASTA  -o /home/user/demo/DEMO/DEMO1/results
```
The analysis will generate the following result:

`Evo_dNdS.png`: A gene-level evolutionary analysis summary generated by the Selection module, showing Evo2 sequence-context scores alongside ω values estimated under the codeml M0 model for candidate genes.

![Selection output](../DEMO/DEMO1/results/Evo_dNdS.png)

Asterisks in the figure indicate that the likelihood ratio test comparing the free-ratio model with the M0 model reached statistical significance, suggesting branch-level heterogeneity in selective pressure for that gene (`*`, P < 0.05; `**`, P < 0.01).

#### **6.2.2 Site-level evolutionary analysis**

After identifying candidate genes with significant branch-level heterogeneity in selective pressure from the gene-level analysis, the SiteView module can be used for site-level evolutionary analysis of individual genes. This module integrates Evo2 site-level sequence scores with codeml site-model results to identify candidate positively selected sites and examine their positional distribution along the sequence.

Here, UGT2 is used as an example. Assuming that the corresponding CDS sequence file is located at `/home/user/DEMO/DEMO2/DATA/UGT2.fasta`, and the phylogenetic tree file is located at `home/user/DEMO/DEMO2/DATA/Rhodiola.tree`, run the following command:

```
phyloselect siteview \
-s /home/user/DEMO/DEMO2/DATA/UGT2.fasta \
-t home/user/DEMO/DEMO2/DATA/Rhodiola.tree \
-o /home/user/demo/DEMO/DEMO2/results
```

After completion, the following results will be generated:

- `site_test_summary.csv`: A summary of codeml site-model analyses, including parameter estimates, likelihood values, and BEB-supported candidate positively selected sites under models such as M0, M3, M7, and M8.

| Model         | np   | lnL      | Parameters                                                   | Positively selected sites                                    |
| ------------- | ---- | -------- | ------------------------------------------------------------ | ------------------------------------------------------------ |
| M0: one ratio | 28   | -5359.51 | ω = 0.244                                                    | Not applicable                                               |
| M3: discrete  | 32   | -5206.25 | p0 = 0.455, ω0 = 0.000;  p1 = 0.444, ω1 = 0.196; p2 = 0.101, ω2 = 1.792 | 66A\*, 77V\*, 86L\*\*, 98A\*\*, 106L\*\*, 151T\*, 166V\*\*, 182M\*\*, 183A\*\*, 193I\*\*, 210V\*\*, 216Y\*\*, 242L\*\*, 245Q\*, 272Q\*\*, 286R\*, 364S\*\*, 421T\*\*, 427V\*, 431A\*, 436D\*\*, 442G\*, 456Y\*, 464D\*, 472T\*\* |
| M7: beta      | 29   | -5223.15 | p = 0.110, q = 0.398                                         | Not allowed                                                  |
| M8: beta & ω  | 31   | -5206.19 | p0 = 0.911;  p = 0.383, q = 2.948;  p1 = 0.089, ω = 1.893    | 98A\*\*, 182M\*\*, 193I\*, 242L\*, 364S\*                    |

​	Note: `*` and `**` indicate posterior probability support levels from BEB analysis. `*` indicates posterior 	probability > 0.95, and `**` indicates posterior probability > 0.99.

- `UGT2EvolutionarySites.png`: A visualization showing site-level variation patterns in UGT2 within a phylogenetic framework, including Evo2 site-level score distributions and candidate positively selected sites.

  ![SiteView output](../DEMO/DEMO2/results/UGT2/plots/UGT2EvolutionarySites.png)

  The lower track in the figure shows posterior probabilities of positive selection across codon sites under the M8 model. Red markers indicate candidate positively selected sites identified by Bayes Empirical Bayes (BEB) analysis, where `*` indicates posterior probability > 0.95 and `**` indicates posterior probability > 0.99.

#### **6.2.3 Environmental association analysis**

After identifying candidate genes with branch-level heterogeneity in selective pressure, the EnvAssoc module can be used to explore whether these evolutionary signals are associated with environmental factors. This module first uses aBSREL to detect episodic diversifying selection signals along branches in the phylogenetic tree, and then applies phylogenetic generalized least squares (PGLS) to test statistical associations between branch-level selection metrics and environmental variables, thereby providing ecological context for adaptive evolution.

Here, UGT1 is used as an example. Assuming that the CDS sequence file is located at `/home/user/DEMO/DEMO3/DATA/UGT1.fasta`, the phylogenetic tree file at `/home/user/DEMO/DEMO3/DATA/UGT1.tree`, and the environmental data table at `/home/user/DEMO/DEMO3/DATA/R_env.csv`, run:

```bash
phyloselect envassoc \
-s /home/user/DEMO/DEMO3/DATA/UGT1.fasta \
-t /home/user/DEMO/DEMO3/DATA/UGT1.tree \
-e /home/user/DEMO/DEMO3/DATA/R_env.csv \
-o /home/user/DEMO/DEMO3/results
```
After completion, the following outputs will be generated:

1. `Table1_aBSREL_branch_summary.csv`: Summary table of branch-level aBSREL results.

| Taxon        | ω_weighted | LRT    | P-value | Q-value | Selection status |
|--------------------|-----|----------|--------------------------|------------------|--------------------|
| R. kirilowii | 20.26      | 37.574 | <0.001  | <0.001  | Excluded         |
| R. tibetica  | 1e+10      | 0.503  | 0.329   | 1       | Excluded         |
| R. hobsonii  | 0.584      | 4.126  | 0.047   | 1       | Retained         |
| R. amabilis  | 0.3849     | 0      | 1       | 1       | Retained |

Note:  `ω_weighted` represents the branch-level weighted ω value estimated by aBSREL, and LRT denotes the likelihood ratio test statistic. P-value and Q-value represent the raw significance level and the FDR-adjusted significance level, respectively. `Selection status` summarizes both statistical significance and quality-control outcomes, with three possible categories: `Significant` indicates that the branch has a Q-value < 0.05 and passes quality filtering; `Retained` indicates that the branch passes quality filtering but does not reach statistical significance; `Excluded` indicates that the branch is excluded from downstream interpretation due to extreme ω estimates, unstable dN/dS values, or weak substitution signals. Although `R. kirilowii` showed a significant aBSREL signal, it was excluded from subsequent interpretation because its branch-level ω estimate was unstable.

2. `Table2_PGLS_environment_association.csv`: Results of PGLS environmental association analysis.

| Environmental factor | N    | β        | SE      | t      | P-value | Q-value | Significant after FDR |
| -------------------- | ---- | -------- | ------- | ------ | ------- | ------- | --------------------- |
| GHIL_0.05            | 16   | -3.25    | 0.8891  | -3.656 | 0.003   | 0.327   | FALSE                 |
| bio_5_0.05           | 16   | -0.03814 | 0.01414 | -2.697 | 0.017   | 0.785   | FALSE                 |
| bio_10_0.05          | 16   | -0.0446  | 0.01708 | -2.611 | 0.021   | 0.785   | FALSE                 |
| GPD_0.05             | 16   | -0.05588 | 0.0235  | -2.378 | 0.032   | 0.785   | FALSE                 |

Note: N represents the number of branches or taxa included in the analysis for a given environmental variable; β represents the regression coefficient; SE represents the standard error; t represents the t statistic; P-value represents the raw significance level; and Q-value represents the FDR-adjusted significance level. `Significant after FDR` indicates whether the environmental factor remains statistically significant after multiple-testing correction.

In this example, all environmental variables have `Significant after FDR = FALSE`, indicating that no environmental associations remained significant after FDR correction in the UGT1 dataset.

**Supplementary example: saikosaponin biosynthesis-related gene in *Bupleurum chinense***

Because the UGT1 example did not yield significant environmental associations after FDR correction, an additional example using a saikosaponin biosynthesis-related gene from *Bupleurum chinense* is provided to illustrate EnvAssoc outputs when significant signals are present.

Run:

```
phyloselect envassoc \
-s /home/user/DEMO/DEMO3/DATA/Bupleurum_chinense.fasta \
-t /home/user/DEMO/DEMO3/DATA/Bupleurum_chinense.tree \
-e /home/user/DEMO/DEMO3/DATA/B_env.csv \
-o /home/user/DEMO/DEMO3/results
```

After the analysis is completed, the following outputs will be generated:

1. `Table1_aBSREL_branch_summary.csv`: Summary table of branch-level aBSREL results. Representative entries are shown below.

| Taxon | ω_weighted | LRT    | P-value | Q-value | Selection status |
| ----- | ---------- | ------ | ------- | ------- | ---------------- |
| S81   | 1.607      | 19.988 | <0.001  | <0.001  | Significant      |
| B51   | 65.44      | 9.303  | 0.003   | 0.133   | Excluded         |
| F61   | 0.4905     | 4.747  | 0.034   | 1       | Retained         |
| F41   | 1e+10      | 3.593  | 0.061   | 1       | Excluded         |

2. `Table2_PGLS_environment_association.csv`: Results of the PGLS environmental association analysis, showing statistical associations between branch-level selection metrics and environmental variables. Representative entries are shown below.

| Environmental factor | N    | β         | SE        | t      | P-value | Q-value | Significant after FDR |
| -------------------- | ---- | --------- | --------- | ------ | ------- | ------- | --------------------- |
| bio_16               | 17   | 0.001425  | 0.0002451 | 5.813  | <0.001  | <0.001  | TRUE                  |
| bio_18               | 17   | -0.001214 | 0.0005024 | -2.417 | 0.029   | 0.099   | FALSE                 |
| bio_8                | 17   | -0.04023  | 0.01672   | -2.406 | 0.029   | 0.099   | FALSE                 |

In this supplementary case, aBSREL detected positive selection on branch S81, and PGLS analysis showed that precipitation of the wettest quarter (`bio_16`) remained significant after FDR correction. This demonstrates that EnvAssoc can identify candidate environmental factors associated with gene evolutionary changes when stronger branch-level selection signals are present.

#### **6.2.4 Molecular docking analysis**

The Docking module is used to compare relative docking score patterns between candidate proteins and their substrates, products, or cofactors at the structural level. This module takes protein structure files and a docking configuration file as input, and organizes docking results across different species and genes within a phylogenetic framework to facilitate structural interpretation of candidate genes.

Here, UGT genes are used as an example. Assuming that the docking configuration file is located at `/home/user/DEMO/DEMO4/DATA/docking_config.csv`, the phylogenetic tree file at `/home/user/DEMO/DEMO4/DATA/UGT.tree`, and receptor structure files are stored under `/home/user/DEMO/DEMO4/DATA/receptors`, run:

```bash
phyloselect docking \
-c /home/user/DEMO/DEMO4/DATA/docking_config.csv \
-t /home/user/DEMO/DEMO4/DATA/UGT.tree \
-o /home/user/DEMO/DEMO4/results
```

The following main outputs will be generated:

1. `LigandBindingProfile.png`: A visualization of overall docking score patterns between candidate proteins and specified substrates/products. This figure can be used to compare relative ligand-binding structural patterns across genes or species.

​	![docking output](../DEMO/DEMO4/results/LigandBindingProfile.png)

2. `docking_summary.csv`: A summary table of docking scores for candidate proteins against substrates, products, and optional cofactors, which can be used for detailed comparison and downstream analysis.




---

