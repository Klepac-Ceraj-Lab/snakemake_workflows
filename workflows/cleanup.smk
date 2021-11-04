(samples, lanes) = glob_wildcards(os.path.join(kneadfolder, "{sample}_kneaddata_{otherstuff}.fastq.gz"))
samples = list(set(samples))
samples.sort()


rule cleanup_kneaddata:
    input:
        fwd = os.path.join(kneadfolder, "{sample}_kneaddata_paired_1.fastq"),
        rev = os.path.join(kneadfolder, "{sample}_kneaddata_paired_2.fastq"),
    output:
        fwd = os.path.join(kneadfolder, "{sample}_kneaddata_paired_1.fastq.gz"),
        rev = os.path.join(kneadfolder, "{sample}_kneaddata_paired_2.fastq.gz"),
    run: 
        shell("gzip {input.fwd}")