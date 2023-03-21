# General settings

Adjust `config.yaml` to configure the workflow execution. The config file `samples.tsv` contain the following sample information:

- `samples.tsv` lists the samples as well as the read sets for each sample, with one set of reads for each `sample_id-unit`. Note that `sample_id-unit` combinations should be unique.

The pipeline will call all samples individually but will output the gVCF files required for joint genotyping. These files can be genotypyes with for example GLnexus.
