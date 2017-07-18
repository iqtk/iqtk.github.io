**We use [GitHub issues](https://github.com/iqtk/iqtk/issues) for
tracking requests and bugs.**

## Installation

The latest `iqtk` docker image can be obtained as follows:

```bash
docker pull gcr.io/jbei-cloud/iqtk:0.0.1
```

The toolkit can also be pip installed directly but if using this approach we suggest doing so within a virtual environment like `conda` or `virtualenv`.

```bash
pip install iqtk
```

## Examples

Most users of this system will be concerned with writing and running analytical workflows. The following examples provide a starting point for these concerns.

### Running a workflow

The core workflows can be run as simply as the following,

```bash
iqtk run expression --config=path/to/your/run_conf.json
```

provided a config file with the necessary input files. The following is an example of this for the RNA-seq analysis workflow:

```json
{
  "_meta": {
    "workflow": "core:expression"
  },
  "ref_fasta": "gs://cflow-public/data/genomes/Drosophila_melanogaster/Ensembl/BDGP5.25/Sequence/BowtieIndex/genome.fa",
  "genes_gtf": "gs://cflow-public/data/genomes/Drosophila_melanogaster/Ensembl/BDGP5.25/Annotation/Archives/archive-2015-07-17-14-30-26/Genes/genes.gtf",
  "cond_a_pairs": [
    ["gs://cflow-public/data/rnaseq/downsampled_reads/GSM794483_C1_R1_1_small.fq",
     "gs://cflow-public/data/rnaseq/downsampled_reads/GSM794483_C1_R1_2_small.fq"],
    ["gs://cflow-public/data/rnaseq/downsampled_reads/GSM794484_C1_R2_1_small.fq",
     "gs://cflow-public/data/rnaseq/downsampled_reads/GSM794484_C1_R2_2_small.fq"],
    ["gs://cflow-public/data/rnaseq/downsampled_reads/GSM794485_C1_R3_1_small.fq",
     "gs://cflow-public/data/rnaseq/downsampled_reads/GSM794485_C1_R3_2_small.fq"]
   ],
   "cond_b_pairs": [
     ["gs://cflow-public/data/rnaseq/downsampled_reads/GSM794486_C2_R1_1_small.fq",
      "gs://cflow-public/data/rnaseq/downsampled_reads/GSM794486_C2_R1_2_small.fq"],
     ["gs://cflow-public/data/rnaseq/downsampled_reads/GSM794487_C2_R2_1_small.fq",
      "gs://cflow-public/data/rnaseq/downsampled_reads/GSM794487_C2_R2_2_small.fq"],
     ["gs://cflow-public/data/rnaseq/downsampled_reads/GSM794488_C2_R3_1_small.fq",
      "gs://cflow-public/data/rnaseq/downsampled_readsGSM794488_C2_R3_2_small.fq"]
    ]
}

```

### Writing a workflow

To illustrate the structure (and hopefully simplicity) of building new workflows, per one of the core objectives of the project, the following  example (a simplified version of the full RNA-seq workflow) is provided:

```python

class TranscriptomicsWorkflow(Workflow):
    def __init__(self):
        """Initialize a workflow."""
        self.tag = 'tuxedo-transcriptomics'
        self.arg_template = [details omitted]
        super(TranscriptomicsWorkflow, self).__init__()

    def define(self):
        p, args = self.p, self.args

        # For each condition, create a PCollection to store the input read pairs.
        reads_a = util.fc_create(p, args.cond_a_pairs)

        th_a = ops.tophat(reads_a,
                          ref_fasta=args.ref_fasta,
                          args=args,
                          genes_gtf=args.genes_gtf,
                          tag='cond_b')

        # The same, but for condition B.
        th_b = ops.tophat(reads_b,
                          ref_fasta=args.ref_fasta,
                          args=args,
                          genes_gtf=args.genes_gtf,
                          tag='cond_a')

        # Subset the outputs of the tophat steps to obtain only the bam (alignment)
        # files. Then combine the collections.
        align_a = util.match(th_a, {'file_type': 'bam'})
        align_b = util.match(th_b, {'file_type': 'bam'})
        align = util.combine(p, (align_a, align_b))

        # For each set of reads, perform a transcriptome assembly with cufflinks,
        # yielding one gtf feature annotation for each input read set.
        cuff = ops.cufflinks(align, args=args)

        return cuff

```

Workflows are composed of operations, the basic building block that can be shared and remixed to build new workflows. A trivial illustrative example is the following:

```python
class TopHat(task.ContainerTask):

    def __init__(self, args, tag=None):
        container = task.ContainerTaskResources(
            disk=60, cpu_cores=4, ram=8,
            image='gcr.io/jbei-cloud/tophat:0.0.1')
        super(TopHat, self).__init__(task_label='dummy', args=args,
                                     container=container)

    def process(self, input_file):

        cmd = util.Command(['cat', localize(input_file), '>',
                            self.out_path + '/file.txt'])

        yield self.submit(cmd.txt, inputs=[input_file],
                          expected_outputs=[{'txt': 'file.txt'}])
```

As one can see significant benefit can be derived by moving from simpler domain-specific workflow languages to the more expressive language of python.

For more detailed examples of workflows and operations check out any of those provided as part of the core toolkit, e.g. [the one for RNA-seq analysis](inquiry/toolkit/rna_quantification/workflow.py).

### Data schema

A key objective of the project is to provide consistent delivery of the data resulting from workflow runs to databases according to a controlled and standardized schemas. Significant effort on the part of the Global Alliance for Genomics and Health (GA4GH) is underway in this area. Here we make use of lightly adapted versions of [those schemas](https://github.com/ga4gh/ga4gh-schemas/tree/master/src/main/proto/ga4gh). The following is an example of the schema used by the RNA-seq analysis workflow above:

```python
message DiffExpressionLevel {
  option (gen_bq_schema.table_name) = "differential_expression";
  string id = 1;
  string geneid = 2;
  string gene = 3;
  string locus = 4;
  string sample1 = 5;
  string sample2 = 6;
  string status = 7;
  float expression1 = 8;
  float expression2 = 9;
  float lnFoldChange = 10;
  float testStatistic = 11;
  float pValue = 12;
  float qValue = 13;
  bool significant = 14;
}
```

Browse the full schema [here].

## Syncing data from instruments

Included in the toolkit we have a sync daemon that currently is a simple wrapper of the `gcloud rsync ...` utility. The plan is to extend this with the ability to sync to various other providers and to be an intermediate for instrument logs as they are generated, beyond just syncing resulting data.

The uplink daemon can be initiated as follows:

```bash
iqtk uplink --local_path=/my/source/dir \
            --remote_path=/my/target/uri \
            --service_account=[your service account address] \
            --service_account_key_path=[your sa key path] \
            --sleep_time=600
```

Make sure the serivce account you create has read and write access to the target bucket and logs; documentation [here](https://cloud.google.com/storage/docs/authentication). The following is an example:

```bash
# example
uplink -l /path/to/example -r gs://[your-bucket]/[your-target-path] -r 30
```

### Acknowledgements

We would like to acknowledge the value of input received from members of the Google Genomics team (summarized this [post](https://opensource.googleblog.com/2016/11/docker-dataflow-happier-workflows.html)). See also [DockerFlow](https://github.com/googlegenomics/dockerflow) for a Java implementation of Airflow-style container workflow orchestration with Beam.

### Contact

Want to get in touch? You can [provide feedback](https://goo.gl/forms/2cOmuUrQ3n3CKpim1) regarding this or other documentation, [reach out to us](https://goo.gl/forms/j8FWdNJqABAoJvcW2) regarding collaboration, or [request a new feature or analytical capability](https://goo.gl/forms/dQm3SDcoNZsV7AAd2).

Read more about the Joint BioEnergy Institute (JBEI) at https://www.jbei.org/.

© Regents of the University of California, 2017. Licensed under a BSD-3 <a href="https://github.com/.../blob/master/LICENSE">license</a>.
