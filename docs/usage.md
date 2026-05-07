# nfdata-omics/r2t2p: Usage

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

This page describes the main inputs required to run R2T2P, including the RNA-seq and Ribo-seq samplesheet,
reference annotation inputs, and the optional setup required for proteomics searches.

## Samplesheet input

You will need to create a samplesheet with information about the RNA-seq and Ribo-seq libraries you would like to
analyse before running the pipeline. Use `--input` to specify its location:

```bash
--input '[path to samplesheet file]'
```

The samplesheet must be a comma-separated file with a header row. R2T2P uses the `library_type` column to route
libraries through RNA-seq-specific and Ribo-seq-specific steps, and the `condition` column to define contrasts for
differential analyses when `--control_label` is provided.

### Multiple runs of the same sample

Use the same `sample` and `library_type` values when the same library has been sequenced more than once, for
example across multiple lanes. The pipeline will concatenate these raw reads before downstream analysis. Repeated
runs for the same `sample` and `library_type` must all be either single-end or paired-end.

Below is an example for the same paired-end RNA-seq sample sequenced across 3 lanes:

```csv title="samplesheet.csv"
sample,fastq_1,fastq_2,library_type,condition
CONTROL_RNA_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz,RNA,control
CONTROL_RNA_REP1,AEG588A1_S1_L003_R1_001.fastq.gz,AEG588A1_S1_L003_R2_001.fastq.gz,RNA,control
CONTROL_RNA_REP1,AEG588A1_S1_L004_R1_001.fastq.gz,AEG588A1_S1_L004_R2_001.fastq.gz,RNA,control
```

### Full samplesheet

The pipeline will auto-detect whether each library is single-end or paired-end from the presence of `fastq_2`.
Single-end libraries should leave `fastq_2` empty.

A samplesheet containing RNA-seq and Ribo-seq libraries may look like this:

```csv title="samplesheet.csv"
sample,fastq_1,fastq_2,library_type,condition
CONTROL_RNA_REP1,/path/to/control_rna_rep1_R1.fastq.gz,/path/to/control_rna_rep1_R2.fastq.gz,RNA,control
CONTROL_RNA_REP2,/path/to/control_rna_rep2_R1.fastq.gz,/path/to/control_rna_rep2_R2.fastq.gz,RNA,control
CONTROL_RIBO_REP1,/path/to/control_ribo_rep1.fastq.gz,,Ribo,control
CONTROL_RIBO_REP2,/path/to/control_ribo_rep2.fastq.gz,,Ribo,control
TREATED_RNA_REP1,/path/to/treated_rna_rep1_R1.fastq.gz,/path/to/treated_rna_rep1_R2.fastq.gz,RNA,treated
TREATED_RIBO_REP1,/path/to/treated_ribo_rep1.fastq.gz,,Ribo,treated
```

| Column         | Required | Description                                                                                                                                                                                                |
| -------------- | -------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `sample`       | Yes      | Custom sample or library name. This value cannot contain spaces. Use the same value for repeated sequencing runs that should be concatenated.                                                              |
| `fastq_1`      | Yes      | Full path to read 1. Files must be gzipped and end with `.fastq.gz`, `.fq.gz`, `.fasta.gz`, or `.fa.gz`.                                                                                                   |
| `fastq_2`      | No       | Full path to read 2 for paired-end libraries. Files must be gzipped and end with `.fastq.gz`, `.fq.gz`, `.fasta.gz`, or `.fa.gz`. Leave empty for single-end libraries.                                    |
| `library_type` | Yes      | Library type used to route reads through the workflow. Use `RNA` for RNA-seq libraries and `Ribo` for Ribo-seq libraries.                                                                                  |
| `condition`    | For DE   | Experimental condition label. If differential analyses are requested, `--control_label` must match one of the values in this column, and all other conditions will be compared against that control label. |

Use the examples above as templates for preparing your own samplesheet.

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run nfdata-omics/r2t2p \
    --input ./samplesheet.csv \
    --outdir ./results \
    --fasta ./genome.fa.gz \
    --gtf ./annotation.gtf.gz \
    --control_label control \
    -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.
You can provide `--gff` instead of `--gtf`, but only one of the two annotation formats should be supplied.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

> [!WARNING]
> Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/running/run-pipelines#configuring-pipelines), other infrastructural tweaks (such as output directories), or module arguments (args).

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run nfdata-omics/r2t2p -profile docker -params-file params.yaml
```

with:

```yaml title="params.yaml"
input: './samplesheet.csv'
outdir: './results/'
fasta: './genome.fa.gz'
gtf: './annotation.gtf.gz'
control_label: 'control'
<...>
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

## Optional proteomics setup

The Protein module is optional. It is enabled only when `--fragpipe_manifest` is provided. When proteomics is
enabled, R2T2P builds protein databases from the reference annotation and ORFquant results, adds decoys and
contaminants with Philosopher, and runs FragPipe searches against the generated databases.

Due to licensing restrictions, some FragPipe companion tools cannot be redistributed inside the pipeline
container. Before running the Protein module, download the required tools manually and make them available to the
pipeline with `--fragpipe_tools_folder` and `--fragpipe_diann_folder`.

The tools folder must contain the following JAR files in versioned subdirectories:

```text
tools/
  MSFragger-<version>/
    MSFragger-<version>.jar
  IonQuant-<version>/
    IonQuant-<version>.jar
  diaTracer-<version>/
    diaTracer-<version>.jar
```

R2T2P validates these paths with the following patterns:

```text
<TOOLS_DIR>/MSFragger-*/MSFragger-*.jar
<TOOLS_DIR>/IonQuant-*/IonQuant-*.jar
<TOOLS_DIR>/diaTracer-*/diaTracer-*.jar
```

The DIA-NN folder must contain the Linux executable named `diann-linux`:

```text
diann/
  diann-linux
```

If needed, make the DIA-NN binary executable before launching the pipeline:

```bash
chmod +x /path/to/diann/diann-linux
```

For academic use, the external tools can be downloaded from their official sources:

- MSFragger: <https://msfragger.nesvilab.org/>
- IonQuant: <https://ionquant.nesvilab.org/>
- diaTracer: <https://diatracer.nesvilab.org/>
- DIA-NN: <https://github.com/vdemichev/DiaNN>

To enable the Protein module, provide the FragPipe manifest, the FragPipe workflow file, the required annotation
file, and the local folders containing the external tools:

```bash
nextflow run nfdata-omics/r2t2p \
    --input ./samplesheet.csv \
    --outdir ./results \
    --fasta ./genome.fa.gz \
    --gtf ./annotation.gtf.gz \
    --fragpipe_manifest ./proteomics_manifest.tsv \
    --fragpipe_workflow ./fragpipe.workflow \
    --fragpipe_annotation ./tmt_annotation.txt \
    --fragpipe_tools_folder /path/to/tools \
    --fragpipe_diann_folder /path/to/diann \
    -profile docker
```

The FragPipe manifest is a tab-delimited file with no header row. Each row corresponds to one LC-MS/MS run and
must contain four columns in this order:

```text
<path_to_LC-MS_file>    <experiment_name>    <bioreplicate>    <data_type>
```

For example:

```tsv title="proteomics_manifest.tsv"
/storage/run_01.mzML	exp_a	1	DDA
/storage/run_02.mzML	exp_a	2	DDA
/storage/run_03.mzML	exp_b	1	DDA
/storage/run_04.mzML	exp_b	2	DDA
```

Use paths that are accessible from the machine or cluster where Nextflow is running.

### FragPipe annotation file

The `--fragpipe_annotation` parameter is currently required whenever `--fragpipe_manifest` is supplied. For TMT
workflows, provide the channel annotation file expected by FragPipe/TMT-Integrator. The file should contain two
whitespace-delimited columns: the TMT channel and the sample label assigned to that channel.

For example:

```text title="tmt_annotation.txt"
126	control_rep1
127N	control_rep2
127C	control_rep3
128N	treated_rep1
128C	treated_rep2
129N	treated_rep3
129C	NA
130N	NA
```

Use `NA` for channels that should not be assigned to a sample. During execution, R2T2P creates a copy of this file
for each `experiment_name` in the FragPipe manifest and appends the experiment name to every non-`NA` sample label.
For example, if the manifest contains `exp_a`, the label `control_rep1` becomes `control_rep1_exp_a` in the
per-experiment annotation file used by FragPipe. Because this suffix is added automatically, the labels in
`--fragpipe_annotation` should normally be the base sample labels without the experiment suffix.

For non-TMT workflows, the parameter is still required by the current pipeline validation. If the selected FragPipe
workflow does not use a TMT annotation file, provide a small placeholder file that can be read by the pipeline, for
example:

```text title="fragpipe_annotation_placeholder.txt"
126	NA
```

The FragPipe workflow file can be generated with FragPipe or downloaded from the FragPipe workflow collection.
R2T2P rewrites the `database.db-path` entry in the workflow file so that FragPipe searches use the protein databases
generated by the pipeline.

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nfdata-omics/r2t2p
```

### Reproducibility

It is a good idea to specify the pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nfdata-omics/r2t2p releases page](https://github.com/nfdata-omics/r2t2p/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducibility, you can use share and reuse [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

> [!TIP]
> If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.

## Core Nextflow arguments

> [!NOTE]
> These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen)

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

> [!IMPORTANT]
> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to check if your system is supported, please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer environment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://charliecloud.io/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `wave`
  - A generic configuration profile to enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above (requires Nextflow ` 24.03.0-edge` or later).
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the pipeline steps, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher resources request (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#set-max-resources) and [customise process resources](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#customize-process-resources) section of the nf-core website.

### Custom Containers

In some cases, you may wish to change the container or conda environment used by a pipeline steps for a particular tool. By default, nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However, in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#update-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#modifying-tool-arguments) section of the nf-core website.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
