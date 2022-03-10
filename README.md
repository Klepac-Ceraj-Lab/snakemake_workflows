# Snakemake Workflows

by  
  Lauren Tso  
  Kevin Bonham, PhD


## Initial setup

```sh
$ poetry install
$ conda install -c bioconda trimmomatic trf
```

Use `poetry shell` with this directory as the working directory every time you're ready to use.


### bioBakery databases

```sh
$ metaphlan --install --bowtie2db $METAPHLAN_DATABASE_DIR # be sure to set this in `config.yaml` as well
$ kneaddata_database --download human_genome bowtie2 $KNEADDATA_DATABASE_DIR
$ humann_databases --download uniref uniref90_diamond $HUMANN_DATABASE_DIR
$ humann_databases --download chocophlan full $HUMANN_DATABASE_DIR
$ humann_databases --download utility_mapping full $HUMANN_DATABASE_DIR
```

## Usage

1. Make a copy of `config_template.yaml` called `config.yaml`.
2. Make a copy of `cluster_config.yaml` called `cluster.yaml`.
3. Set the variables for your environment.
4. Run:

```sh
$ snakemake -s /home/vklepacc/software/repos/snakemake_workflows/run.smk \
    --configfile config.yaml --cluster-config cluster.yaml \
    --cluster "sbatch -n {cluster.processors} -N 1 -t {cluster.time} --mem {cluster.memory} -o output/logs/{rule}-%j.out -e output/logs/{rule}-%j.err -p newnodes" \
    --jobs 16 --latency-wait 15
```
