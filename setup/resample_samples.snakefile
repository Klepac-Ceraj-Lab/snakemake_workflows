(samples,) = glob_wildcards(os.path.join(input_folder, "{sample}.fastq.gz"))
samples = list(set(samples))
samples.sort()