# PhyloSelect

An Integrative Framework for Multi-Dimensional Phylogenetic and Functional Gene Analysis

## Introduction

PhyloSelect is an integrated command-line software package for phylogenetic analysis and molecular evolutionary research. Within a unified workflow, it integrates target CDS recovery, gene-level and site-level selection analysis, environmental association analysis, and structure-level comparison. It is designed to provide a consistent, reusable, and extensible workflow for multi-level evolutionary analysis of protein-coding genes.

In addition to classical molecular evolutionary methods, PhyloSelect supports sequence-context features extracted by the large-scale genomic foundation model Evo2, including sequence scores and site entropy metrics. These features can be combined with traditional codon-model results from tools such as codeml to assist in identifying candidate genes, candidate sites, and localized evolutionary heterogeneity.

PhyloSelect consists of five functional modules, covering the main steps from target protein-coding sequence acquisition to selection signal detection, environmental association analysis, and structure-level comparison. These modules can be used together as a recommended workflow or run independently according to specific research aims.

| **Module**    | **Description**                                              |
| ------------- | ------------------------------------------------------------ |
| **GeneMiner** | Recovers target genes from raw sequencing data, including shallow sequencing data. |
| **SiteView**  | Performs site-level evolutionary analysis for a single protein-coding gene by combining Evo2 entropy scores and codeml site models to identify local variation patterns and candidate positively selected sites. |
| **Selection** | Performs gene-level selection analysis for multiple protein-coding genes and compares the M0 and free-ratio models to evaluate branch-level heterogeneity in selective pressure. |
| **EnvAssoc**  | Tests statistical associations between branch-level selection signals and environmental variables using aBSREL and PGLS. |
| **Docking**   | Compares relative docking score patterns between candidate proteins and substrates or products based on predicted protein structures and molecular docking results, providing structure-level evidence for evaluating potential differences. |

![PhyloSelect workflow](pipeline.png)

For the **complete user manual**, please refer to [User Manuall](./manual/manual-EN.md).

---

## Module Overview

The table below provides a brief overview of the input file types required for each module. For detailed file formats, parameter descriptions, and complete usage instructions, please refer to the [User Manuall](./manual/manual-EN.md).

| **Module**    | **Main input data**                                          | **Common file formats** |
| ------------- | ------------------------------------------------------------ | ----------------------- |
| **GeneMiner** | Sequencing data; target gene references; sample list         | FASTQ, FASTA, TSV       |
| **Selection** | Multiple protein-coding sequences; TreeMap; phylogenetic trees | FASTA, CSV, Newick      |
| **SiteView**  | Single protein-coding sequence; phylogenetic tree            | FASTA, Newick           |
| **EnvAssoc**  | Codon alignment; phylogenetic tree; environmental table      | FASTA, Newick, CSV      |
| **Docking**   | Protein structures; docking configuration; phylogenetic tree | PDB/CIF, CSV, Newick    |

---

## Installation

### System Requirements

PhyloSelect currently supports Linux and macOS systems. The following environment is recommended:

- Operating system: Linux or macOS
- Python: ≥ 3.10
- Memory: ≥ 8 GB recommended
- Disk space: 3–5 GB of available space recommended for the software, bundled components, and runtime cache, excluding user-generated analysis data

### Option 1. Conda installation

```bash
conda create -n phyloselect
conda activate phyloselect
conda install evanstone::phyloselect
```

### Option 2. Source installation

```bash
git clone https://github.com/evanstone/phyloselect.git
cd phyloselect
conda env create -f environment.yaml
conda activate phyloselect
pip install .
```

### Installation verification:

``````bash
# Check program installation
phyloselect --version
# Check dependencies and service availability
phyloselect check
``````

If `phyloselect check` reports that all required dependencies are available, the current environment is ready for analysis. Otherwise, please check whether the Conda environment has been correctly activated, whether all dependencies have been installed, and whether internet access is available for modules that require remote services.

---

## Quick Start

Before running the quickstart examples, please prepare the example data in your working directory.

- If PhyloSelect was installed via Conda, copy the built-in `quickstart/` dataset to the current directory:

```
cp -r $(python -c "import phyloselect, pathlib; print(pathlib.Path(phyloselect.__file__).parent / 'quickstart')") .
```

- If PhyloSelect was installed from the GitHub repository, enter the repository directory directly:

```
cd phyloselect
```

### Example 1: SiteView analysis

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
-  `SiteTestSummary.csv` : Summarizes analysis results across different site models, useful for identifying potential positively selected sites.

### Example 2: Selection analysis

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

### Example 3: EnvAssoc analysis

**Example Data:**

- CDS sequence file: `quickstart/sequences/gene3.fasta`
- Environmental trait matrix: `quickstart/config/env_traits.csv`
-  Corresponding phylogenetic tree file: `quickstart/trees/test1.nwk`

**Run command:**

```bash
phyloselect envassoc \
-s quickstart/sequences/gene3.fasta \
-e quickstart/config/env_traits.csv \
-t quickstart/trees/test1.nwk \
-o outputdir
```

**Main Outputs:**

- `Table1_aBSREL_significant_branches.csv` : significant branches identified by aBSREL as showing evidence of episodic positive selection.
- `Table2_PGLS_environment_association.csv` : association results between branch-level selection signals and environmental variables based on PGLS analysis.

### Example 4: Docking analysis

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

---

## Full Demo Data

The `quickstart/` directory is intended for installation testing and small-scale demonstration.Complete demo datasets are provided in the `DEMO/` directory. For detailed instructions on how to run each demo, interpret the outputs, and understand the analysis results, please refer to the [PhyloSelect User Manual](manual.md).

---

## Citation

If you use PhyloSelect or any of its modules, please cite the corresponding tools and methods used in your analysis:

1. Yi Shi et al. PhyloSelect: An Integrated Workflow for Multi-Level Phylogenetic and Functional Gene Analysis.
2. Yu X, Tang Z, Zhang Z, et al. GeneMiner2: Accurate and Automated Recovery of Genes From Genome Skimming Data. *Molecular Ecology Resources*. 2026;26(2):e70111. doi:10.1111/1755-0998.70111.
3. Brixi G, Durrant MG, Ku J, et al. Genome modelling and design across all domains of life with Evo 2. *Nature*. 2026. doi:10.1038/s41586-026-10176-5.
4. Yang Z. PAML 4: Phylogenetic Analysis by Maximum Likelihood. *Molecular Biology and Evolution*. 2007;24(8):1586–1591. doi:10.1093/molbev/msm088.
5. Kosakovsky Pond SL, Poon AFY, Velazquez R, et al. HyPhy 2.5: A customizable platform for evolutionary hypothesis testing using phylogenies. *Molecular Biology and Evolution*. 2020;37(1):295–299. doi:10.1093/molbev/msz197.
6. Smith MD, Wertheim JO, Weaver S, et al. Less Is More: An Adaptive Branch-Site Random Effects Model for Efficient Detection of Episodic Diversifying Selection. *Molecular Biology and Evolution*. 2015;32(5):1342–1353. doi:10.1093/molbev/msv022.
7. Trott O, Olson AJ. AutoDock Vina: improving the speed and accuracy of docking with a new scoring function, efficient optimization, and multithreading. *Journal of Computational Chemistry*. 2010;31(2):455–461. doi:10.1002/jcc.21334.
8. Minh BQ, Schmidt HA, Chernomor O, et al. IQ-TREE 2: New Models and Efficient Methods for Phylogenetic Inference in the Genomic Era. *Molecular Biology and Evolution*. 2020;37(5):1530–1534. doi:10.1093/molbev/msaa015.
9. Capella-Gutiérrez S, Silla-Martínez JM, Gabaldón T. trimAl: a tool for automated alignment trimming in large-scale phylogenetic analyses. *Bioinformatics*. 2009;25(15):1972–1973. doi:10.1093/bioinformatics/btp348.
10. Edgar RC. Muscle5: High-accuracy alignment ensembles enable unbiased assessments of sequence homology and phylogeny. *Nature Communications*. 2022;13:6968. doi:10.1038/s41467-022-34630-w.
11. Suyama M, Torrents D, Bork P. PAL2NAL: robust conversion of protein sequence alignments into the corresponding codon alignments. *Nucleic Acids Research*. 2006;34(Web Server issue):W609–W612. doi:10.1093/nar/gkl315.
12. O’Boyle NM, Banck M, James CA, et al. Open Babel: An open chemical toolbox. *Journal of Cheminformatics*. 2011;3:33. doi:10.1186/1758-2946-3-33.
13. Le Guilloux V, Schmidtke P, Tuffery P. Fpocket: An open source platform for ligand pocket detection. *BMC Bioinformatics*. 2009;10:168. doi:10.1186/1471-2105-10-168.



## License

PhyloSelect is released under the MIT License. See [LICENSE](LICENSE) for details.
