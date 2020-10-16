# Snakemake Workflows

by  
  Lauren Tso  
  Kevin Bonham, PhD

## Usage

1. Make a copy of `config_template.yaml` called `config.yaml`.
2. Make a copy of `cluster_config.yaml` called `cluster.yaml`.
3. Set the variables for your environment.
4. Run:

```sh
$ snakemake -s /home/vklepacc/software/repos/snakemake_workflows/run.snakefile \
    --configfile config.yaml --cluster-config cluster.yaml \
    --cluster "sbatch -n {cluster.processors} -N 1 -t {cluster.time} --mem {cluster.memory} -o output/logs/{rule}-%j.out -e output/logs/{rule}-%j.err -p newnodes" \
    --jobs 16 --latency-wait 15
```
