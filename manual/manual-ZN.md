# PhyloSelect用户手册

## 1. 关于 PhyloSelect
### 1.1 软件简介
PhyloSelect 是一款面向系统发育分析与分子进化研究的集成化命令行软件。它在统一工作框架下整合目标 CDS 序列恢复、基因层与位点层选择分析、环境关联分析以及结构层面比较等功能，旨在为蛋白编码基因的多层次进化分析提供一致、可复用且便于扩展的流程。

除经典分子进化分析方法外，PhyloSelect 还支持引入大规模基因组基础模型 Evo2 提取的序列上下文特征，包括序列评分和位点熵值等指标，并将其与 codeml 等传统密码子模型结果结合，用于辅助识别候选基因、候选位点和局部进化异质性。

### 1.2 主要功能

PhyloSelect 由五个功能模块组成，覆盖目标蛋白编码序列获取、选择信号检测、环境关联分析和结构层面比较等主要分析步骤。各模块既可以按照推荐流程组合运行，也可以根据具体研究目的独立使用。

| 模块          | 用途说明                                                     |
| ------------- | ------------------------------------------------------------ |
| **GeneMiner** | 用于从原始测序数据（支持浅层测序）中提取目标基因             |
| **SiteView**  | 对单个蛋白编码基因开展位点层进化分析，结合 Evo2 熵值分数和 codeml 位点模型识别局部变异模式及候选正选择位点 |
| **Selection** | 对多个蛋白编码基因开展基因层选择分析，比较 M0 与 free-ratio 模型以评估分支层选择压力异质性 |
| **EnvAssoc**  | 基于 aBSREL 和 PGLS 检验分支层选择信号与环境变量之间的统计关联 |
| **Docking**   | 基于预测蛋白结构和分子对接结果，比较候选蛋白与底物或产物的相对对接评分模式，为结构层面的差异评估提供参考 |

### 1.3 典型分析流程

在常见分析场景中，PhyloSelect 可按照以下流程使用：

**GeneMiner → SiteView / Selection → EnvAssoc → Docking**

首先，使用 **GeneMiner** 准备目标 CDS 序列；随后使用 **SiteView** 或 **Selection** 模块分析位点或分支层面的选择信号；若需要进一步讨论环境因素，可继续使用 **EnvAssoc** 模块进行关联分析；若需要从结构层面补充证据，则可对候选基因使用 **Docking** 模块开展分子对接分析。

该流程仅作为一般参考。实际使用时，用户可根据已有输入数据选择起始模块；例如，已有 CDS 序列时可直接从 Selection 或 SiteView 开始分析，已有结构文件时也可直接运行 Docking 模块。



---


## 2 安装与运行环境

本章介绍 PhyloSelect 的运行环境要求及安装方法。为简化部署过程，PhyloSelect 已通过 Conda 环境配置或软件包方式管理主要依赖组件，用户通常无需手动配置多数第三方软件。当前软件支持通过 Conda 仓库安装或从 GitHub 源代码安装两种方式。


### 2.1 推荐运行环境

PhyloSelect 当前支持 Linux 和 macOS 系统，建议在满足以下条件的环境中运行：

- 操作系统：Linux 或 macOS
- Python：≥ 3.10
- 内存：建议 ≥ 8 GB
- 磁盘空间：建议预留 3–5 GB 可用空间（用于软件本体、封装组件及运行时缓存，不包含用户分析产生的数据）

PhyloSelect 的部分模块，尤其是 Evo2 相关序列评分与特征提取功能，需要通过 API 调用远程模型服务。因此，在运行相关模块时必须具备互联网访问权限。

在高性能计算平台（HPC）中，计算节点常存在外网访问受限的情况。若计算节点默认禁网，则 Evo2 相关模块可能无法正常运行。建议用户在任务提交前确认节点网络策略，或优先在具备外网访问能力的环境中完成相关分析。


### 2.2 安装方法

用户可根据需求选择以下任意一种方法完成安装：

#### 方式一：通过 Conda 仓库安装

适用于希望快速完成部署的用户，建议在独立环境中安装。

```bash
conda create -n phyloselect python=3.10
conda activate phyloselect
conda install evanstone::phyloselect
```

#### 方式二：从 GitHub 源码安装

该方式适用于需要获取最新版本、调试程序或参与开发的场景。安装时依赖仓库中提供的 `environment.yaml`文件构建运行环境。

```bash
git clone https://github.com/evanstone/phyloselect.git
cd phyloselect
conda env create -f environment.yaml
conda activate phyloselect
pip install .
```


### 2.3 安装验证

安装完成后，可使用以下命令检查程序是否可正常调用，以及运行环境是否配置完整：
```bash
# 1. 确认程序是否正确安装
phyloselect --version

# 2. 检查主要依赖和相关服务状态
phyloselect check
```

**提示：** 若 `check` 命令返回状态为 `Ready`，则说明当前运行环境已完成配置，可用于正式分析，否则请根据提示检查 Conda 环境是否正确激活、依赖是否完整安装，以及网络是否可访问相关服务。



---


## 3 输入文件格式

本章介绍 PhyloSelect 在不同分析模块中的输入数据要求及基本使用规则。由于各模块面向的分析任务不同，其输入数据类型与组织方式存在差异。建议用户在运行具体模块前，仔细确认输入数据格式是否符合要求。

### 3.1 模块输入类型概览

PhyloSelect 的不同模块面向不同分析阶段，因此所需输入文件也有所不同。下表概括了各模块的主要输入数据和常见文件格式。更详细的文件命名规则、字段要求和参数说明见各模块对应章节。

| **模块**      | **主要输入数据**                                             | **常见文件格式**     |
|---|---|---|
| **GeneMiner** | 测序数据文件；目标基因参考序列文件；样本列表文件             | FASTQ、FASTA、TSV    |
| **Selection** | 多个蛋白编码基因 CDS 序列文件；TreeMap 映射表；系统发育树文件 | FASTA、CSV、Newick   |
| **SiteView**  | 单个蛋白编码基因 CDS 序列文件；系统发育树文件                | FASTA、Newick        |
| **EnvAssoc**  | 密码子比对文件；系统发育树文件；环境变量表                   | FASTA、Newick、CSV   |
| **Docking**   | 蛋白结构文件；Docking 配置文件（包含参考蛋白 ID、底物与产物信息、配体路径及对接参数等）；系统发育树文件 | PDB/CIF、CSV、Newick |

### 3.2 不同模块的输入要求

#### 3.2.1 GeneMiner 模块

GeneMiner 模块用于在同源参考序列指导下，从测序数据中恢复目标基因序列。在 PhyloSelect 中，该模块主要用于获取目标蛋白编码序列（CDS），并为后续 Selection 和 SiteView 等模块提供输入数据。默认设置下，GeneMiner 将依次执行 reads 过滤、目标基因组装、侧翼序列修剪和结果合并等步骤，以生成目标 CDS 数据集。

GeneMiner 模块需要提供 TSV 格式的样本列表文件和FASTA 格式的目标基因参考序列。样本列表的具体格式是`<物种名><Tab><数据文件1>`（单端）或者`<物种名><Tab><数据文件1><Tab><数据文件2>`（双端），每一行代表一个样本。其中的数据文件建议采用绝对路径。

以双端测序数据为例，样本列表配置如下：

```
Bupleurum_chinense	/home/user/GeneMiner2/DEMO/DEMO3/DATA/PLANT/Bupleurum_chinense_1.fq.gz	/home/user/GeneMiner2/DEMO/DEMO3/DATA/PLANT/Bupleurum_chinense_2.fq.gz
Bupleurum_fruticosum	/home/user/GeneMiner2/DEMO/DEMO3/DATA/PLANT/Bupleurum_fruticosum_1.fq.gz	/home/user/GeneMiner2/DEMO/DEMO3/DATA/PLANT/Bupleurum_fruticosum_2.fq.gz
Bupleurum_krylovianum	/home/user/GeneMiner2/DEMO/DEMO3/DATA/PLANT/Bupleurum_krylovianum_1.fq.gz	/home/user/GeneMiner2/DEMO/DEMO3/DATA/PLANT/Bupleurum_krylovianum_2.fq.gz
Bupleurum_malconense	/home/user/GeneMiner2/DEMO/DEMO3/DATA/PLANT/Bupleurum_malconense_1.fq.gz	/home/user/GeneMiner2/DEMO/DEMO3/DATA/PLANT/Bupleurum_malconense_2.fq.gz
Bupleurum_wenchuanense	/home/user/GeneMiner2/DEMO/DEMO3/DATA/PLANT/Bupleurum_wenchuanense_1.fq.gz	/home/user/GeneMiner2/DEMO/DEMO3/DATA/PLANT/Bupleurum_wenchuanense_2.fq.gz
Bupleurum_yunnanense	/home/user/GeneMiner2/DEMO/DEMO3/DATA/PLANT/Bupleurum_yunnanense_1.fq.gz	/home/user/GeneMiner2/DEMO/DEMO3/DATA/PLANT/Bupleurum_yunnanense_2.fq.gz
```

参考序列需要放在一个独立的文件夹下，每个基因可以有一条或多条参考序列。对于每个基因，将对应的所有参考序列存在`<基因名>.fasta`下。假设需要提取matK和psbA两个基因，则需要在一个空白文件夹下创建`matK.fasta`和`psaA.fasta`两个文件，分别保存对应的参考序列。

假设样本列表保存在`/home/user/samples.tsv`，被子植物353参考基因保存在`/home/user/Angiosperms353`，期望的输出文件夹为`/home/user/output`，用默认设置运行GeneMiner2的命令如下：：

```bash
phyloselect geneminer -f /home/user/samples.tsv -r /home/user/Angiosperms353 -o /home/user/output
```

#### 3.2.2 SiteView 模块

SiteView 模块用于单基因的位点层分析，输入为一个 CDS 序列文件，**输入序列应为标准编码序列，长度为 3 的倍数，且序列内部不应包含提前终止密码子**。用户还可提供 Newick 格式的系统发育树文件；若未提供，则可通过`--infer-tree`自动完成建树。

下面给出一个典型示例命令：

```bash
phyloselect siteview -s DEMO/DEMO1/DATA/rbcL.fasta -t DEMO/DEMO1/DATA/rbcL.treefile -o DEMO/DEMO1/results
```

#### 3.2.3 Selection 模块

Selection 模块用于在多个基因上开展分支层面的选择压力分析。该模块的基本输入为一个包含多个 CDS 序列文件的目录，其中每个序列文件对应一个独立基因。程序将以该目录中的序列文件作为分析对象，进行批量的分支层进化分析。目录中的文件名将作为后续结果目录和任务标识使用。

例如，输入目录`DEMO/DEMO2/DATA/FASTA`的结构如下：
```text
.
├── 4HPAAS.fasta
├── 4HPAR1.fasta
└── 4HPAR2.fasta
```
默认情况下，Selection 模块将自动对每个输入基因进行多序列比对、系统发育树构建以及分支层选择压力分析，适用于用户仅提供目标 CDS 序列而未预先准备系统发育树的场景。

如果用户已经具备预先构建好的系统发育树，也可以直接提供树文件进行分析。在这种情况下，除树文件外，还需要同时提供一个 `tree_file_map` 文件，用于指定每个序列文件对应的系统发育树。程序将根据该映射关系，对对应基因使用用户指定的树拓扑开展分析。若 `tree_file_map` 中某个基因未指定树文件，或指定的树文件不存在，则程序将自动为该基因重新构建系统发育树。

`tree_file_map` 文件采用CSV格式，至少应包含以下两列：

- sequence：输入序列文件路径
- tree：对应的系统发育树文件路径

为保证分析结果的正确性，序列文件中的物种名称应与对应树文件中的标签保持一致；`tree_file_map` 文件中的 sequence 和 tree 路径也应与实际文件位置严格匹配。

例如，`DEMO/DEMO2/DATA/tree_file_map.csv`内容如下：

```text
sequence,tree
DEMO/DEMO2/data/4HPAAS.fasta,DEMO/DEMO2/tree_file/4HPAAS.treefile
DEMO/DEMO2/data/4HPAR1.fasta,DEMO/DEMO2/tree_file/4HPAR1.treefile
DEMO/DEMO2/data/4HPAR2.fasta,DEMO/DEMO2/tree_file/4HPAR2.treefile
```

假设目标 CDS 序列目录为`home/user/DEMO/DEMO2/DATA/FASTA`，系统发育树映射文件为`/home/user/DEMO/DEMO2/DATA/tree_file_map.csv`，输出目录为：`/home/user/output`：使用用户提供的系统发育树运行 Selection 模块的命令如下：

```bash
phyloselect selection -i home/user/DEMO/DEMO2/DATA/FASTA -m /home/user/DEMO/DEMO2/DATA/tree_file_map.csv -o /home/user/output
```

如果不提供 `tree_file_map` 文件，程序将自动为所有输入基因构建系统发育树并完成后续选择压力分析：

```
phyloselect selection -i home/user/DEMO/DEMO2/DATA/FASTA -o /home/user/output
```

#### 3.2.4 EnvAssoc 模块

EnvAssoc 模块用于在系统发育框架下分析分支层选择信号与环境因子之间的关系。该模块以密码子比对文件、系统发育树文件和环境数据表为输入，首先调用 aBSREL 识别分支层选择信号，随后使用 PGLS 检验分支层选择指标与环境变量之间的统计关联。

该模块的基本输入包括一个密码子比对文件、一个与该比对文件对应的系统发育树文件，以及一个环境数据表。密码子比对文件支持 FASTA 或 PHYLIP 格式；系统发育树文件应为 Newick 格式，且树中的末端标签需与比对文件中的序列名称保持一致。环境数据表采用 CSV 格式，且必须包含名为 `taxon` 的列，用于与树末端标签和序列名称进行匹配。

例如，环境数据表 `DEMO/DEMO3/DATA/env.csv` 的部分内容如下：

```text
taxon,elevation,bio_15,GHFD,WRB2_CODE,AWC
Rhodiola_quadrifida,1626,60.41747,21,18,0
Rhodiola_crenulata,990,78.01089,9,19,159
Rhodiola_tibetica,3909,92.29301,56,8,150
```

其中，`taxon` 列为必需列，其内容应与系统发育树末端标签及比对文件中的序列名称完全一致；其余列为待检验的环境变量，可包括海拔、气候变量、土壤因子或其他生态因子。运行分析前，建议用户检查环境变量是否存在缺失值、非数值字符或明显异常值，以避免影响 PGLS 分析结果。

下面给出一个典型示例命令：

```bash
phyloselect envassoc -s DEMO/DEMO3/DATA/4HPAAS3.fasta -t DEMO/DEMO3/DATA/4HPAAS3.treefile -e DEMO/DEMO3/DATA/env.csv -o DEMO/DEMO3/results
```

若需要调整 aBSREL 或 PGLS 的分析行为，用户还可进一步设置遗传密码表、待检验分支范围、多重替换模型、是否启用同义位点速率变异（SRV）、显著性阈值以及最小物种数等参数。有关这些参数的详细说明，见后续参数章节。


#### 3.2.5 Docking 模块

Docking 模块用于从结构层面比较不同基因或不同物种对应蛋白的相对配体结合模式，为候选基因或候选位点的功能差异解释提供参考。该模块以蛋白结构文件和对接配置文件为输入，分别针对指定底物和产物进行分子对接，并结合系统发育树对不同谱系中的对接结果进行组织和展示。

该模块的基本输入包括一个对接配置文件和一个系统发育树文件。对接配置文件采用 CSV 格式，用于指定每个基因对应的受体结构目录、底物与产物信息，以及可选的辅助因子和参考结构信息。系统发育树用于在系统发育背景下展示不同物种的对接结果。

对接配置文件中每一行对应一个基因，至少应包含以下列：

- `Gene`：基因名称或标识
- `ReceptorDir`：该基因对应的蛋白结构文件目录
- `Substrate`：底物名称
- `Product`：产物名称

此外，配置文件还可包含以下可选列：

- `Cofactor`：辅助因子名称
- `Reference`：参考结构编号或参考配体信息

其中，`Reference` 列可省略。若未提供 `Reference`，程序将在后续分析中使用自动口袋识别策略搜索潜在结合口袋。

例如，配置文件 `DEMO/DEMO4/DATA/docking_config.csv` 的部分内容如下：

```text
Gene,ReceptorDir,Substrate,Product,Cofactor,Reference
UGT1,DEMO/DEMO4/DATA/receptors/UGT1,Tyrosol,Salidroside,UDP-Glucose,8ITA
UGT2,DEMO/DEMO4/DATA/receptors/UGT2,Tyrosol,Salidroside,UDP-Glucose,8ITA
UGT3,DEMO/DEMO4/DATA/receptors/UGT3,Tyrosol,Salidroside,UDP-Glucose,8ITA
...
```

其中，ReceptorDir 指向该基因对应的受体结构文件目录。每个目录下通常包含多个物种的结构文件，文件格式为 .cif，每个文件对应一个物种。例如，UGT1 目录可组织为：

```text
UGT1/
├── rhodiola_amabilis.cif
├── rhodiola_bupleuroides.cif
├── rhodiola_chrysanthemifolia.cif
...
```
为保证结果的可比性，建议不同基因目录中的物种集合尽可能一致；同时，结构文件名中的物种名称应与系统发育树中的末端标签保持一致，以便后续在系统发育背景下整合和解释对接结果。

下面给出一个典型示例命令：
```bash
phyloselect docking -c DEMO/DEMO4/DATA/docking_config.csv -t DEMO/DEMO4/DATA/UGT.tree -o DEMO/DEMO4/results
```



---

## 4. 参数说明

PhyloSelect 采用模块化命令行结构，不同模块具有不同的参数组合。本节按照模块分别介绍各参数的含义，并区分必需参数与可选参数。

### 4.1 Common options

以下参数在多个模块中通用，用于控制程序的基本运行行为：

- `--threads`: 指定分析使用的线程数。
- `--force`: 若输出目录已存在，则覆盖原有结果。
- `--keep-intermediate`: 保留运行过程中生成的临时文件和中间文件，便于排错或复查。
- `--skip-plot`: 跳过图形绘制，仅输出表格或模型结果。

### 4.2 GeneMiner 模块

本团队已开发并发布 GeneMiner2 的最新版本。PhyloSelect 中的 GeneMiner 模块仅调用其中的目标 CDS 提取功能，用于为后续 Selection 和 SiteView 分析准备输入序列。因此，本节仅介绍运行该功能所需的三个必需参数。若用户需要使用 GeneMiner2 的图形界面版本或完整功能，可[参考 GeneMiner2 官方说明](https://github.com/sculab/GeneMiner2)。

- `-f FILE`: 样本列表文件，TSV 格式。每行对应一个样本，单端数据格式为 <物种名><Tab><数据文件1>，双端数据格式为 <物种名><Tab><数据文件1><Tab><数据文件2>。
- `-r DIR `: 目标基因参考序列目录。每个目标基因对应一个 FASTA 文件，文件名建议为 <基因名>.fasta。
- `-o DIR`: 输出目录，用于保存 GeneMiner 模块生成的中间文件和目标 CDS 序列结果。

若不指定具体 action，程序将使用默认流程完成目标 CDS 提取。用户也可在命令末尾指定需要执行的步骤，例如：

其他高级参数，如 k-mer 大小、并行线程数、trim 模式、合并策略、多序列比对程序和建树程序等，通常无需在常规 CDS 提取中修改。详细参数说明请参见 最新版[GeneMiner2 的命令行说明](https://github.com/sculab/GeneMiner2/blob/master/manual/ZH_CN/command_line.md)。

### 4.3 SiteView 模块

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
**必需参数**:

- `seq`: 输入单个 CDS FASTA 文件。若输入文件已经是密码子比对结果，可同时指定 --codon-aligned，以跳过比对相关预处理。
- `output_dir`: 输出目录，用于保存分析结果、中间文件和图形结果。

**可选参数**:

- `tree`: 输入 Newick 格式的系统发育树文件。树中的标签应与输入序列名称一致。
- `infer-tree`: 当未提供树文件时，自动构建系统发育树。
- `outgroups`: 指定一个或多个外群物种名称，用于树定根。
- `--conservation-threshold`: 设置树图中叶节点标签显示的最大字符数。
- `codon-aligned`: 声明输入文件已经是密码子比对结果，程序将跳过比对预处理步骤。
- `conservation-threshold`: 位点熵阈值，用于区分相对保守位点与相对变异位点，默认值为 1.4。
- `no-run-paml`: 跳过 codeml 位点模型分析，仅执行 Evo2 相关评分和可视化。
- `entropy-csv`: 提供预先计算好的熵值文件，主要用于调试、复现实验或测试。


### 4.4 Selection 模块

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

**必需参数**：

- `input`: 输入序列文件。该参数可接受两种形式：一个包含多个序列文件的目录，或以分号分隔的多个序列文件路径。程序将把每个序列文件视为一个独立基因进行分析。
- `output-dir`: 输出目录，用于保存分析结果和中间文件。


**可选参数**:

- `tree-file-map`: 系统发育树映射文件，用于指定每个输入序列文件对应的树文件。若未提供该参数，程序将为各基因自动构建系统发育树。映射文件格式见 [3.2.3 节](#323-selection-模块)。

### 4.4 EnvAssoc module

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

**必须参数**:

- `seq`: 输入密码子比对文件，用于 aBSREL 分析。该文件可为 FASTA 或 PHYLIP 格式，并应为适用于分支层分析的密码子比对结果。
- `tree`: 输入 Newick 格式的系统发育树文件。树中的标签应与比对文件中的物种名称一致。
- `env`: 输入环境数据表，格式为 CSV。该表必须包含 taxon 列，其内容应与系统发育树和比对文件中的物种名称一致。输入格式示例见 [3.2.4 节](#324-envassoc-模块) 。
- `output_dir`: 输出目录，用于保存分析结果和中间文件。


**可选参数**:

- `code`: 指定 HyPhy aBSREL 使用的遗传密码表，默认值为 Universal。
- `branches`: 指定 aBSREL 中待检验的分支集合，默认值为 All。
- `multiple-hits`: 指定 aBSREL 使用的多重替换模型。
- `srv`: 是否启用同义位点速率变异（SRV）。
- `timeout-sec`: 设置 aBSREL 的最长运行时间，单位为秒。
- `no-log-transform`: 不对 omega_weighted 进行 log1p 转换后再执行 PGLS 分析。
- `alpha`: 设置 FDR 校正后 Q 值的显著性阈值，默认值为 0.05。
- `min-n`: 设置环境变量纳入 PGLS 分析所需的最小物种数，默认值为 5。


### 4.5 Docking 模块

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

**Required parameters**:
- `docking-config`: 输入对接配置文件，格式为 CSV。该文件用于指定基因名称、受体结构目录、底物、产物以及可选的辅助因子或参考结构信息。配置文件格式示例见 [3.2.5 节](#325-docking-模块)。
- `tree`: 输入 Newick 格式的系统发育树文件。树中的物种名称应与受体结构文件对应的物种名称一致。
- `output-dir`: 输出目录，用于保存对接结果和中间文件。




---
## 5. 输出结果

### 5.1 SiteView 模块输出

分析完成后，SiteView 模块的结果将保存在输出目录对应的任务子目录下。典型目录结构如下：
```text
siteview/
├── file_input/
├── evo_output/
├── paml_output/
├── plots/
```

各子目录的含义如下：
- `file_input/`：保存 SiteView 分析所使用的输入与预处理文件，包括密码子比对结果和系统发育树文件。若用户未提供树文件，程序自动推断得到的树也将保存在该目录中。用于保存分析输入和预处理结果，便于结果复现。
- `evo_output/`：保存 Evo2 计算得到的位点层特征结果，例如每个位点的熵值或相关评分信息。这些结果可用于辅助评估位点的保守性与变异性。
- `paml_output/`：保存基于 PAML/codeml 的位点模型分析结果，包括不同模型的参数估计、对数似然值、LRT 检验结果以及候选正选择位点信息。
- `plots/`：保存 SiteView 自动生成的可视化结果，包括位点层保守性分布图、熵值曲线图或综合进化图等。该目录中的图形文件适合用于结果浏览、比较和展示。

其中，以下文件通常是用户最需要优先关注的核心结果：

| 文件                          | 说明 |
|-----------------------------------|-------------|
| `file_input/example.codon.aln.fasta` | 密码子比对后的序列文件 |
| `file_input/example.paml.tree`      | 输入或自动推断的系统发育树 |
| `evo_output/example_entropy.csv`  | 位点熵值统计结果文件，记录序列比对中各个位点的信息熵及相关变异性指标，用于识别保守位点与高变异位点 |
| `paml_output/site_test_summary.csv` | PAML 位点模型分析结果汇总表，包含不同位点模型的参数估计结果、模型比较统计量，以及正选择位点检测结果等关键信息。 |
| `plots/example_site_evolution.png` | 位点演化模式综合可视化结果，包括系统发育树、氨基酸位点保守性热图以及逐位点变异性分布。 |

建议优先查看`plots/example_site_evolution.png`，快速了解该基因在不同物种中的整体演化模式，包括系统发育关系、位点保守性分布以及变异趋势，再结合 `evo_output/example_entropy.csv` 以了解位点模型检验结果。

### 5.2 Selection 模块输出

分析完成后，Selection 模块的结果将保存在输出目录对应的任务子目录下。典型目录结构如下：

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
各文件和子目录的含义如下：

- `Evo_dNdS.png`：所有基因的 Evo 评分与 dN/dS 模式综合图，用于整体结果概览。
- `gene1/`, `gene2`/等：每个基因对应一个独立结果目录。
  - `file_input/`：该基因的输入与预处理文件。
  - `paml_output/`：该基因的 PAML/codeml 原始输出结果。
  - `result/`：该基因的最终统计结果。
    - `*_NLL_score.csv`：该基因的模型评分或似然结果。
		- `*_omega.csv`：该基因的 ω（dN/dS）估计结果。


建议优先查看 `Evo_dNdS.png` 以获得整体选择模式概览；随后针对感兴趣的基因，进入对应基因目录，重点查看 `result/` 中的 `*_omega.csv` 和 `*_NLL_score.csv` 文件，并结合 `paml_output/` 中的模型输出结果进行进一步分析。


### 5.3 EnvAssoc 模块输出

分析完成后，EnvAssoc 模块的结果将保存在输出目录对应的任务子目录下。典型目录结构如下：
```text
envassoc/
├── *.ABSREL.json
├── TableS1_aBSREL_full_branch_results.csv
├── Table1_aBSREL_branch_summary.csv
└── Table2_PGLS_environment_association.csv
```
各文件的含义如下：
- `*.ABSREL.json`：HyPhy aBSREL 的原始输出文件，包含分支水平的详细模型结果，包括每条分支的 ω 分布、显著性检验统计量等信息。该文件适用于深入解析或二次分析。
- `TableS1_aBSREL_full_branch_results.csv`：aBSREL 分析的完整分支结果表，汇总所有分支的选择信号及统计指标。
- `Table1_aBSREL_branch_summary.csv`： aBSREL 分支层结果汇总表，包含各分支的 `ω_weighted`、LRT、P 值、Q 值及筛选状态。该表用于概览各分支的选择信号，并标记哪些分支被保留用于后续环境关联分析。
- `Table2_PGLS_environment_association.csv`：PGLS 分析结果，记录环境变量与选择信号之间的统计关联关系，包括回归系数、显著性水平等信息。

建议优先查看 `Table1_aBSREL_significant_branches.csv` 以识别显著的选择分支；随后结合 `Table2_PGLS_environment_association.csv` 分析这些选择信号是否与环境变量显著相关；如需进一步验证或深入分析，可参考 `TableS1_aBSREL_full_branch_results.csv` 及 `*.ABSREL.json` 原始结果。


### 5.4 Docking 模块输出

分析完成后，Docking 模块的结果将保存在输出目录对应的任务子目录下。典型目录结构如下：

```text
docking/
├── ligands/
├── receptors/
├── pockets/
├── docking_results/
├── docking_summary.csv
├── CofactorBindingEnergy.png
└── InhibitionDiff.png
```
各文件和子目录的含义如下：
- `ligands/`：配体结构文件目录。
保存输入配置文件中涉及的所有配体文件，包括底物、产物及辅助因子等对应的 `sdf`、`pdb` 和 `pdbqt` 文件，用于后续分子对接分析。
- `receptors/`：受体结构文件目录。
保存所有基因对应的受体结构文件夹，其中包括目标基因以及参考蛋白（若配置文件中提供）。每个基因对应一个独立子目录，目录中包含所有物种对应的 `pdb` 和 `pdbqt` 文件。
- `pockets/`：分子口袋识别结果目录。
保存用于对接分析的结合口袋信息。若配置文件中未提供参考蛋白，程序将默认使用 `fpocket` 自动寻找最合适的分子口袋，并将结果保存在该目录中。
- `docking_results/`：原始分子对接结果目录。
保存所有物种、所有蛋白与所有相关配体组合的分子对接结果，包括不同基因和物种层面的原始打分与中间输出文件。
- `docking_summary.csv`：对接结果汇总表。
汇总所有物种中各蛋白与相应配体结合的结合能结果。
- `CofactorBindingEnergy.png`：辅助因子结合能标准化热图。
当输入配置文件中包含辅助因子信息时，程序将提取各物种、各基因对应的辅助因子结合能，并在归一化后生成热图，用于比较不同基因和不同物种之间辅助因子结合模式的相对差异。
- `InhibitionDiff.png`：底物与产物结合差异的归一化热图。
用于展示不同基因在不同物种中的底物与产物结合差异模式，是比较结构功能分化的重要图形结果。

**核心结果**

| 文件                       | 说明            |
|--------------------------------|------------------------|
| `docking_summary.csv` | 所有物种中各基因对应蛋白与相应配体结合能的结果  |
| `CofactorBindingEnergy.png`    | 辅助因子结合能的归一化热图，仅在提供辅助因子时输出 |
| `InhibitionDiff.png` | 底物与产物结合差异的归一化热图       |

建议首先查看可视化结果以获得整体模式：
- 若配置文件中提供了辅助因子，可优先查看 `CofactorBindingEnergy.png`，评估不同基因和物种之间的辅助因子结合模式；
- 查看 `InhibitionDiff.png`，分析底物与产物结合差异，从而识别潜在的结构功能分化。


在获得整体趋势后，可进一步查看 `docking_summary.csv `以获取具体的结合能数值，并进行精细比较或下游统计分析。 如需追溯具体对接细节或验证结果来源，可进入 `docking_results/ `目录查看对应的原始分子对接输出文件。



---

## 6. 运行示例

本章提供 PhyloSelect 的两类示例：快速上手示例和文章完整案例。

对于仓库内置的 `quickstart/` 示例数据，请先进入 PhyloSelect 仓库根目录后再运行命令，以确保示例数据中的相对路径能够被正确解析。案例复现分析所需数据需提前从 Figshare 下载，并按照对应说明配置输入路径。

### 6.1 快速开始

本节提供 PhyloSelect 四个主要分析模块的最简运行示例，用于帮助用户快速检查软件是否能够正常运行，并初步了解各模块的基本输入和输出。若需要查看更完整的数据说明、结果解释和示例图表，请继续阅读 [6.2 论文完整案例](#6.2-论文完整案例)。

在运行快速开始示例前，请先准备 `quickstart/` 示例数据，并进入 PhyloSelect 仓库根目录。

- 如果通过 Conda 安装 PhyloSelect，可使用以下命令将软件包内置的 `quickstart/` 目录复制到当前目录：

```bash
cp -r $(python -c "import phyloselect, pathlib; print(pathlib.Path(phyloselect.__file__).parent / 'quickstart')") .
```

​	也可以克隆 GitHub 仓库以获取完整示例数据：

```bash
git clone https://github.com/scu-shiyi/PhyloSelect.git
cd PhyloSelect
```

- 如果通过 GitHub 源代码安装 PhyloSelect，则无需再次克隆仓库，直接进入已有的 PhyloSelect 仓库根目录即可。

#### 6.1.1 Selection 模块

**示例数据：**

- 多个基因的 CDS 序列目录：`quickstart/sequences`

该快速示例默认不提供系统发育树文件。若用户已有每个基因对应的系统发育树，可按照前文说明准备序列文件与树文件的映射表，并通过 `--tree-file-map` 参数指定。

**运行命令：**

```bash
phyloselect selection -i quickstart/sequences -o outputdir2
```

**主要输出：**

- `Evo_dNdS.png`：展示不同基因的 Evo 评分与 dN/dS 模式，可用于快速比较多个基因的选择压力差异。
- 每个基因对应的独立结果目录：保存该基因的模型结果和分支选择压力结果。
  - `*_omega.csv`：记录不同分支上的 ω（dN/dS）估计值。
  - `*_NLL_score.csv`：记录不同模型的似然值或模型比较结果。

#### 6.1.2 SiteView 模块

**示例数据：**

- 单个基因的 CDS 序列文件：`quickstart/sequences/gene2.fasta`
- 对应的系统发育树文件：`quickstart/trees/test1.nwk`

**运行命令：**

```bash
phyloselect siteview -s quickstart/sequences/gene2.fasta -t quickstart/trees/test1.nwk -o outputdir1
```

**主要输出：**

- `gene2EvolutionarySites.png`：用于快速查看该基因在不同物种中的整体位点演化模式，包括系统发育关系、位点保守性和变异趋势。
- `SiteTestSummary.csv`：汇总不同位点模型的分析结果，可用于识别潜在正选择位点。

#### 6.1.3 EnvAssoc 模块

**示例数据：**

- CDS 序列文件：`quickstart/sequences/gene2.fasta`
- 环境变量矩阵：`quickstart/config/env_traits.csv`
- 对应的系统发育树文件：`quickstart/trees/test1.nwk`

**运行命令：**

```bash
phyloselect envassoc -s quickstart/sequences/gene2.fasta -e quickstart/config/env_traits.csv -t quickstart/trees/test1.nwk -o outputdir3
```

**主要输出：**

- `Table1_aBSREL_significant_branches.csv`：记录 aBSREL 检测到的显著分支，可用于识别具有 episodic diversifying selection 信号的分支。
- `Table2_PGLS_environment_association.csv`：记录分支层选择信号与环境变量之间的 PGLS 关联分析结果。

#### 6.1.4 Docking 模块

**示例数据：**

- 对接配置文件：`quickstart/config/docking_config.csv`
- 用于结果展示的系统发育树文件：`quickstart/trees/test2.nwk`

对接配置文件应包含目标基因、蛋白结构文件路径、底物、产物、辅助因子以及参考结构等信息。具体字段要求见 [3.2.5 Docking 模块](#325-docking-模块)。

**运行命令：**

```bash
phyloselect docking -c quickstart/config/docking_config.csv -t quickstart/trees/test2.nwk -o outputdir4
```

**主要输出：**

- `docking_summary.csv`：汇总所有受体与相关配体组合的结合能结果。
- `InhibitionDiff.png`：展示底物与产物结合差异的归一化热图，可用于比较不同基因或不同物种之间的结构功能差异。
- `CofactorBindingEnergy.png`：当输入配置文件中包含辅助因子信息时输出，用于展示辅助因子结合能的相对差异。

### 6.2 论文完整案例

本节提供论文中示例分析的完整复现流程。该部分数据不随 GitHub 仓库直接内置，而是通过 Figshare 提供。用户需要先下载 Figshare 中的[示例数据](https://doi.org/10.6084/m9.figshare.32333805.v2)，并根据本节说明进入相应目录运行命令。论文示例以红景天属（*Rhodiola*）UGT 基因为主要数据集，展示 PhyloSelect 在基因层选择分析、位点层进化特征分析、环境关联分析和结构对接比较中的典型输出。

#### 6.2.1 **基因层选择分析**

Selection 模块用于对多个 UGT 基因进行基因层选择分析，比较 M0 与 free-ratio 模型，以评估不同基因是否存在分支层选择压力异质性；同时，该模块引入 Evo2 对序列的对数似然评分，该评分由大规模基因组语言模型基于其学习到的序列上下文规律计算，可作为辅助指标用于比较不同基因的序列约束程度和整体进化特征。


假设多个 UGT 基因的 CDS 序列存放在`home/user/DEMO/DEMO1/DATA/FASTA`，运行以下命令：

```bash
phyloselect selection -i home/user/DEMO/DEMO1/DATA/FASTA  -o /home/user/demo/DEMO/DEMO1/results
```
会得到以下结果：

`Evo_dNdS.png`：Selection 模块生成的基因层进化分析结果，展示候选基因的 Evo2 序列上下文评分以及 codeml M0 模型估计的 ω 值。

[`Selection output`](../images/Evo_dNdS.png)

图中的星号表示 free-ratio 模型与 M0 模型比较的似然比检验达到显著水平，提示该基因在不同分支间存在选择压力异质性（`*`，P < 0.05；`**`，P < 0.01）。

#### 6.2.2 位点层进化分析

在基因层分析识别出具有显著选择压力异质性的候选基因后，可进一步使用 SiteView 模块对单个基因开展位点层进化分析。该模块整合 Evo2 位点级序列评分与 codeml 位点模型分析结果，用于识别候选正选择位点，并观察其在序列中的空间分布模式。

这里以 UGT2 为例，假设对应的 CDS 序列文件存放在 `/home/user/DEMO/DEMO2/DATA/UGT2.fasta`，进化树文件存放在`home/user/DEMO/DEMO2/DATA/Rhodiola.tree`，运行以下命令：

```
phyloselect siteview -s /home/user/DEMO/DEMO2/DATA/UGT2.fasta -t home/user/DEMO/DEMO2/DATA/Rhodiola.tree -o /home/user/demo/DEMO/DEMO2/results
```

运行完成后，会得到以下结果：

- `site_test_summary.csv` ：汇总 codeml 位点模型的分析结果，包括 M0、M3、M7 和 M8 等模型的参数估计、似然值以及 BEB 支持的候选正选择位点。

| Model         | np   | lnL      | Parameters                                                   | Positively selected sites                                    |
| ------------- | ---- | -------- | ------------------------------------------------------------ | ------------------------------------------------------------ |
| M0: one ratio | 28   | -5359.51 | ω = 0.244                                                    | Not applicable                                               |
| M3: discrete  | 32   | -5206.25 | p0 = 0.455, ω0 = 0.000;  p1 = 0.444, ω1 = 0.196; p2 = 0.101, ω2 = 1.792 | 66A\*, 77V\*, 86L\*\*, 98A\*\*, 106L\*\*, 151T\*, 166V\*\*, 182M\*\*, 183A\*\*, 193I\*\*, 210V\*\*, 216Y\*\*, 242L\*\*, 245Q\*, 272Q\*\*, 286R\*, 364S\*\*, 421T\*\*, 427V\*, 431A\*, 436D\*\*, 442G\*, 456Y\*, 464D\*, 472T\*\* |
| M7: beta      | 29   | -5223.15 | p = 0.110, q = 0.398                                         | Not allowed                                                  |
| M8: beta & ω  | 31   | -5206.19 | p0 = 0.911;  p = 0.383, q = 2.948;  p1 = 0.089, ω = 1.893    | 98A\*\*, 182M\*\*, 193I\*, 242L\*, 364S\*                    |

​	注：`*` 和 `**` 表示 BEB 分析中的后验概率支持水平。其中，`*` 表示该位点受到正选择的后验概率 > 0.95，`**` 表示后验概率 > 0.99

- `UGT2EvolutionarySites.png`：展示 UGT2 在系统发育背景下的位点变异模式、Evo2 位点评分分布以及候选正选择位点位置。

  [`SiteView output`](../images/UGT2EvolutionarySites.png)

  图中下方轨道表示 M8 模型下各密码子位点受到正选择的后验概率，红色标记表示通过 Bayes Empirical Bayes（BEB）分析识别出的候选正选择位点，其中 `*` 表示后验概率 > 0.95，`**` 表示后验概率 > 0.99。

#### 6.2.3 环境关联分析

在识别出具有分支层选择压力异质性的候选基因后，可进一步使用 EnvAssoc 模块探索这些进化信号是否与环境因素相关。该模块首先使用 aBSREL 检测系统发育树中各分支上的 episodic diversifying selection 信号，随后结合 PGLS（phylogenetic generalized least squares）分析环境变量与分支层选择指标之间的统计关联，从而为候选基因的适应性演化提供生态环境解释。

这里以 UGT1 为例，假设对应的 CDS 序列文件存放在`/home/user/DEMO/DEMO3/DATA/UGT1.fasta`、系统发育树文件存放在`/home/user/DEMO/DEMO3/DATA/UGT1.tree`，环境数据表存放在`/home/user/DEMO/DEMO3/DATA/R_env.csv`，运行以下命令：

```bash
phyloselect envassoc -s /home/user/DEMO/DEMO3/DATA/UGT1.fasta -t /home/user/DEMO/DEMO3/DATA/UGT1.tree -e /home/user/DEMO/DEMO3/DATA/R_env.csv -o /home/user/DEMO/DEMO3/results
```
运行完成后，会得到以下结果：

1. `Table1_aBSREL_branch_summary.csv`：aBSREL 分支层结果汇总表，下表展示了该文件中的部分结果示例：

| Taxon        | ω_weighted | LRT    | P-value | Q-value | Selection status |
|--------------------|-----|----------|--------------------------|------------------|--------------------|
| R. kirilowii | 20.26      | 37.574 | <0.001  | <0.001  | Excluded         |
| R. tibetica  | 1e+10      | 0.503  | 0.329   | 1       | Excluded         |
| R. hobsonii  | 0.584      | 4.126  | 0.047   | 1       | Retained         |
| R. amabilis  | 0.3849     | 0      | 1       | 1       | Retained |

其中，`ω_weighted` 表示 aBSREL 估计的分支层加权 ω 值，LRT 表示似然比检验统计量。P-value 和 Q-value 分别表示原始显著性水平和 FDR 校正后的显著性水平。`Selection status` 综合反映统计显著性和质量控制结果，包括三种状态：`Significant` 表示该分支 Q-value < 0.05，且通过质量过滤；`Retained` 表示该分支通过质量过滤，但未达到显著水平；`Excluded` 表示该分支因极端 ω 估计、不稳定的 dN/dS 值或替换信号较弱等原因被排除，不用于后续解释。虽然 `R. kirilowii` 显示出显著的 aBSREL 信号，但由于其分支层 ω 估计不稳定，因此被排除在后续解释之外。

2. `Table2_PGLS_environment_association.csv`：PGLS 环境关联分析结果表，展示分支层选择指标与环境变量之间的统计关系，下表展示了该文件中的部分结果示例：

| Environmental factor | N    | β        | SE      | t      | P-value | Q-value | Significant after FDR |
| -------------------- | ---- | -------- | ------- | ------ | ------- | ------- | --------------------- |
| GHIL_0.05            | 16   | -3.25    | 0.8891  | -3.656 | 0.003   | 0.327   | FALSE                 |
| bio_5_0.05           | 16   | -0.03814 | 0.01414 | -2.697 | 0.017   | 0.785   | FALSE                 |
| bio_10_0.05          | 16   | -0.0446  | 0.01708 | -2.611 | 0.021   | 0.785   | FALSE                 |
| GPD_0.05             | 16   | -0.05588 | 0.0235  | -2.378 | 0.032   | 0.785   | FALSE                 |

其中，N 表示纳入该环境变量分析的分支或物种数量，β 表示回归系数，SE 表示标准误，t 表示 t 统计量，P-value 表示原始显著性水平，Q-value 表示 FDR 校正后的显著性水平。`Significant after FDR` 表示该环境因子在多重检验校正后是否仍达到显著水平。

在该示例中，所有环境变量的 `Significant after FDR` 均为 `FALSE`，说明当前 UGT1 数据集中未检测到经过 FDR 校正后仍显著的环境关联信号。

**补充案例：北柴胡柴胡皂苷相关基因**

由于 UGT1 示例在 FDR 校正后未检测到显著环境关联信号，本文进一步提供一个北柴胡（*Bupleurum chinense*）柴胡皂苷生物合成相关基因的补充案例，用于展示 EnvAssoc 在存在可检测信号时的输出形式。

运行以下命令：

```
phyloselect envassoc -s /home/user/DEMO/DEMO3/DATA/Bupleurum_chinense.fasta -t /home/user/DEMO/DEMO3/DATA/Bupleurum_chinense.tree -e /home/user/DEMO/DEMO3/DATA/B_env.csv -o /home/user/DEMO/DEMO3/results
```

运行完成后，结果如下：

1. `Table1_aBSREL_branch_summary.csv`：aBSREL 分支层结果汇总表，下表展示了该文件中的部分结果示例：

| Taxon | ω_weighted | LRT    | P-value | Q-value | Selection status |
| ----- | ---------- | ------ | ------- | ------- | ---------------- |
| S81   | 1.607      | 19.988 | <0.001  | <0.001  | Significant      |
| B51   | 65.44      | 9.303  | 0.003   | 0.133   | Excluded         |
| F61   | 0.4905     | 4.747  | 0.034   | 1       | Retained         |
| F41   | 1e+10      | 3.593  | 0.061   | 1       | Excluded         |

2. `Table2_PGLS_environment_association.csv`：PGLS 环境关联分析结果表，展示分支层选择指标与环境变量之间的统计关系，下表展示了该文件中的部分结果示例：

| Environmental factor | N    | β         | SE        | t      | P-value | Q-value | Significant after FDR |
| -------------------- | ---- | --------- | --------- | ------ | ------- | ------- | --------------------- |
| Environmental factor | N    | β         | SE        | t      | P-value | Q-value | Significant after FDR |
| bio_16               | 17   | 0.001425  | 0.0002451 | 5.813  | <0.001  | <0.001  | TRUE                  |
| bio_18               | 17   | -0.001214 | 0.0005024 | -2.417 | 0.029   | 0.099   | FALSE                 |
| bio_8                | 17   | -0.04023  | 0.01672   | -2.406 | 0.029   | 0.099   | FALSE                 |

在该补充案例中，aBSREL 分析检测到分支 S81 具有正选择信号，随后 PGLS 分析显示 precipitation of the wettest quarter（bio_16）在 FDR 校正后仍达到显著水平。该结果说明，在存在更明确分支层选择信号的数据集中，EnvAssoc 可以进一步用于筛选与候选基因进化变化相关的潜在环境因子。

#### 6.2.4 分子对接分析

Docking 模块用于从结构层面比较候选蛋白与底物、产物或辅助因子之间的相对对接评分模式。该模块以蛋白结构文件和对接配置文件为输入，结合系统发育树组织不同物种和不同基因的对接结果，用于辅助解释候选基因的结构层面差异。

这里以 UGT 基因为例，假设对接配置文件存放在 `/home/user/DEMO/DEMO4/DATA/docking_config.csv`，系统发育树文件存放在 `/home/user/DEMO/DEMO4/DATA/UGT.tree`，受体结构文件按基因存放在 `/home/user/DEMO/DEMO4/DATA/receptors` 目录下，运行以下命令：

```bash
phyloselect docking -c /home/user/DEMO/DEMO4/DATA/docking_config.csv -t /home/user/DEMO/DEMO4/DATA/UGT.tree -o /home/user/DEMO/DEMO4/results
```

运行完成后，将得到以下主要结果：

1. ``LigandBindingProfile.png``：展示候选蛋白与指定底物、产物之间的整体对接评分模式。该图可用于比较不同基因或不同物种在配体结合相关结构特征上的相对差异。

​	[`docking output`](../images/LigandBindingProfile.png)

2. `docking_summary.csv`：汇总各候选蛋白与底物、产物及辅助因子的对接评分结果，可用于查看具体数值并进行后续比较分析。




---

