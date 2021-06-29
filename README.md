# PEACH

**P**harmacogenomic **E**valuator **A**nd **C**aller of **H**aplotypes (PEACH) is a pharmacogenomics tool developed for 
the [Hartwig Medical Foundation pipeline](https://github.com/hartwigmedical/pipeline5). 
It imports haplotypes and related variants from a curated JSON file, reports the presence of these variants in a 
germline VCF, and infers the simplest combination of haplotypes that explains the presence of these variants. 

The two output files are:
* A file that contains the determined genotype of the sample for each gene in the JSON, expressed in terms of haplotypes.
* A file that contains calls for all positions of variants in the JSON file, including annotation and filters wrt both v37 and v38 reference genomes.

## Contents

* [Installation](#installation)
* [Arguments](#arguments)
  + [Mandatory Arguments](#mandatory-arguments)
  + [Optional Arguments](#optional-arguments)
* [Input](#input)
  + [VCF](#vcf)
  + [JSON](#json)
  + [Datastore file locations](#datastore-file-locations)
* [Output](#output)
  + [Genotype TSV file](#genotype-tsv-file)
  + [Calls TSV file](#calls-tsv-file)
* [Algorithm](#algorithm)
  + [Preparation](#preparation)
  + [Get_VCF Variant Calls](#get-vcf-variant-calls)
  + [Annotate Calls with Panel Information](#annotate-calls-with-panel-information)
  + [Infer Haplotypes](#infer-haplotypes)
  + [Examples](#examples)
  + [Restrictions](#restrictions)
* [Tests](#tests)
* [Version History and Download Links](#version-history-and-download-links)

## Installation
PEACH has been designed to work with Python 3.6.

1. Download this repository. See for instance [Github instructions](https://docs.github.com/en/github/creating-cloning-and-archiving-repositories/cloning-a-repository).
2. Install the requirements in one fo the following ways:
    1. Generate a local Python 3.6 virtual environment and install the requirements:
    ```bash
    $ python3.6 -m venv [path/to/new/virtual/environment, for example: ./peach]
    $ source [path/to/new/venv, for example: ./peach/bin/activate]
    (peach) $ pip install -r requirements.txt
    ```
    2. Install the requirements directly:
    ```bash
    $ pip install -r requirements.txt
    ```
3. Run PEACH.

## Arguments
If you have installed PEACH's requirements into a venv, then remember to source the venv before running `peach`.

#### Example Usage
```
(peach) $ peach \
    --vcf input.vcf.gz \
    --sample_r_id COLO829R \
    --sample_t_id COLO829T \
    --panel /path/to/panel.json 
    --outputdir /path/to/outputdir/ \
    --version 1.0 \
```

### Mandatory Arguments
Long&nbsp;Argument | Short Argument | Description
---|---|---
--vcf | -i | Path to germline VCF file of sample. For instance the germline VCF output from PURPLE. Calls should be wrt v37. Support for v38 calls is experimental.
--panel | -p | Path to a JSON file that contains the variants and haplotypes to test on.
--outputdir | -o | Directory to write the output to.
--sample_t_id | -t | The tumor sample ID of the run. Only used for the names of the output files.
--sample_r_id | -r | The ref sample ID of the run.
--tool_version | -v | The version of PEACH. It is included in the output files.

### Optional Arguments
Currently, the only optional argument enables an experimental feature.

Long&nbsp;Argument | Short Argument | Options | Default | Description
---|---|---|---|---
--vcf_reference_assembly_version | -a | V37, V38 | V37 | The version of the reference assembly wrt which the vcf has been constructed. Support for V38 is experimental. Default is V37.

## Input
### VCF
PEACH has been designed to work with VCF files that follow the VCF Version 4.2 format, see 
[specification](https://samtools.github.io/hts-specs/VCFv4.2.pdf). In addition to the required fields, 
each row should contain an annotation field "ANN", as described in 
[annotation format specification](http://snpeff.sourceforge.net/VCFannotationformat_v1.0.pdf), that contains
the subsections "Gene Name" and "HGVS.c". One of the samples in the VCF should have a label 
equal to the `sample_r_id` argument,
and the "GT" sub-field for this sample should be included and filled in with diploid calls.

The calls in the VCF should be with respect to a v37 reference genome. Support for calls wrt v38 is experimental.

### JSON
For an example of a valid panel JSON (with fake data), see 
[example](https://github.com/hartwigmedical/scripts/blob/master/peach/src/test_resources/test_panel.json).
All fields in the example JSON are required. Additional fields are ignored. 
Relevant differences between the v37 and v38 reference sequences for a gene should be included as an entry in the "variants" field
of that gene where the "referenceAlleleV37" and "referenceAlleleV38" fields are different. The set of rs id's with such entries 
has to be equal to the set of rs id's with entries in the "refSeqDifferenceAnnotations" field of that gene.
This ensures that the variant annotation for variants at these locations can be translated to v38 properly

PEACH does not (properly) support panel JSON files that contain (partially) overlapping genes.
Variants in a panel JSON file are not allowed to (partially) overlap.

### Datastore file locations 
(Only relevant for internal use)

Panel:
* Smaller panel for DPYD with haplotypes and haplotypes restricted to those in SOC tests (`/data/common/dbs/peach/panelfiles/min_DPYD.json`).
* Panel with common DPYD haplotypes and variants (`/data/common/dbs/peach/panelfiles/DPYD.json`).

## Output
PEACH outputs two TSV files. One contains genotypes/haplotypes for each gene, the other contains calls for all variants from the panel JSON.
### Genotype TSV file
Name: `[sample_t_id].peach.genotype.tsv`

Column | Example Value | Description
---|---|---
gene | DPYD | Gene for which this haplotype is called.
haplotype | *1 | Called haplotype. If no haplotype could be called, has value "Unresolved Haplotype".
zygosity | HOM | Whether haplotype call is homozygous (HOM) or heterozygous (HET). If no haplotype could be called, has value "N/A".
function | No function | Functionality of this haplotype. Wild type has function "Normal Function". If no haplotype could be called, has value "Unknown Function".
linked_drugs | 5-Fluoracil;Capecitabine | Drugs for which this haplotype is relevant, separated by ";".
url_prescription_info | https://www.some_url.com/5-Fluoracil;https://www.some_other_url.com/Capecitabine | For each listed drug, a url with information on how to translate abnormal haplotype function into an appropriate treatment adjustment. Separated by ";".
panel_version | DPYDpanel_v1.3 | Name and version of panel JSON. Both are taken from fields in the JSON.
repo_version | 1.0 | Version of PEACH.

### Calls TSV file
Name: `[sample_t_id].peach.calls.tsv`

Column | Example Value | Description
---|---|---
gene | DPYD | Gene to which the variant is related.
chromosome | 1 | Chromosome of variant.
position_v37 | 98348885 | Position on chromosome wrt v37 reference genome. If unknown, has value "UNKNOWN".
position_v38 | 97883329 | Position on chromosome wrt v38 reference genome. If unknown, has value "UNKNOWN".
ref_v37 | G | Reference allele wrt v37. If unknown, has value "UNKNOWN".
ref_v38 | A | Reference allele wrt v38. If unknown, has value "UNKNOWN".
allele1 | A | First of the called alleles. Order of alleles is lexicographical order.
allele2 | A | Second of the called alleles. Order of alleles is lexicographical order.
rsid | rs1801265 | Rs id(s) of variant. If more than one, then they are separated by ";". Taken from VCF if available. If not, taken from matching variant in panel JSON, if match exists. If not, has value ".".
variant_annotation_v37 | 85C>T | Variant annotation wrt v37. See [Get_VCF Variant Calls](#get-vcf-variant-calls) and [Annotate Calls with Panel Information](#annotate-calls-with-panel-information) for details.
filter_v37 | PASS | Has value PASS, NO_CALL, UNKNOWN or INFERRED_PASS. See [Get_VCF Variant Calls](#get-vcf-variant-calls) and [Annotate Calls with Panel Information](#annotate-calls-with-panel-information) for details.
variant_annotation_v38 | REF_CALL | Variant annotation wrt v38. See [Get_VCF Variant Calls](#get-vcf-variant-calls) and [Annotate Calls with Panel Information](#annotate-calls-with-panel-information) for details.
filter_v38 | NO_CALL | Has value PASS, NO_CALL, UNKNOWN or INFERRED_PASS. See [Get_VCF Variant Calls](#get-vcf-variant-calls) and [Annotate Calls with Panel Information](#annotate-calls-with-panel-information) for details.
panel_version | DPYDpanel_v1.3 | Name and version of panel JSON. Both are taken from fields in the JSON.
repo_version | 1.0 | Version of PEACH.

## Algorithm
Haplotypes are commonly defined wrt a v38 reference genome. 
If the calls in the input VCF are wrt v37, then this requires a translation of v37 calls to v38 calls.
PEACH extracts the required knowledge of the differences between these reference genomes from the information in the panel JSON.
No matter wrt what reference genome version the calls in the input VCF have been defined, 
PEACH always tries to translate these calls to the other reference genome by using the information in the panel JSON, 
resulting in *full calls* with information for both versions. 
For haplotype calling we use the v38 information from the full calls. 

Define the *VCF RG version* as the reference genome version wrt which the calls in the VCF are defined. 
Define the *non-VCF RG version* as the reference genome version that is not the VCF RG version.

In broad strokes, PEACH does the following:
* Extract relevant calls from VCF, where relevance is determined by the panel JSON.
* Use information from the panel JSON to translate these calls to the other reference genome, where possible, and construct full calls.
* For each gene:
  + Determine for each variant how often each alt allele occurs (alt wrt v38).
  + Determine the unique simplest combination of haplotypes that completely explains that combination of alt alleles and counts.
If there is no unique simplest combination of haplotypes that completely explains the combination of alt alleles and counts, then declare "Unresolved Haplotype".
* Create output files.

### Preparation
First, the panel JSON is loaded and checked for consistency. 

### Get VCF Variant Calls
The calls for the sample `sample_r_id` are extracted from the VCF, and they are compared to the variants in the panel JSON file. 
Calls are ignored when none of the following are true:
* At least one of the rs id's of the call matches an rs id from the panel JSON.
* At least one of the covered positions of the call matches at least one of the covered positions of the variants in the panel JSON 
  for the VCF RG version.
In this comparison, the *covered positions* of a call or variant are the positions of the bases in the reference allele.

Let's call the remaining variants the *observed VCF calls*. 

For each variant in the panel JSON for which there is no matching observed VCF call, 
a call is created that is homozygously the reference allele, an *inferred VCF call*. 
More specifically, a call is added for a panel variant when there are no observed VCF calls such that either:
* The variant rs id matches one of the calls rs id's.
* The covered positions of the variant (partially) match the covered positions of the call (for the VCF RG version).

The observed and inferred VCF calls together form the list of calls that will be considered by PEACH, 
which we will call the *combined VCF calls*.

The annotation and filter of the combined VCF calls are determined in the following way:

Condition | Variant Annotation VCF RG Version | Filter VCF RG Version
---|---|---
Homozygous alt or heterozygous observed call. | From HGVS.c field in ANN in VCF | PASS  
Homozygous reference observed call. | REF_CALL | PASS
Inferred call. | REF_CALL | NO_CALL

### Annotate Calls with Panel Information
For each of the combined VCF calls, an attempt is made to find a variant in the panel JSON that has the same VCF position and reference allele.

If successful:
* If the rs id of the call has not been set, then it is set to the value from the matching variant.
* The reference allele and position wrt the non-VCF RG version are determined from the matching variant.

If unsuccessful:
* Set reference allele and position wrt non-VCF RG version to "UNKNOWN".

Also, the correct annotation and filter wrt the non-VCF RG version are determined according to the following table, 
where an asterisk (*) indicates that any value is allowed:

Matches variant in panel JSON | Type of call | Call is homozygous reference wrt non-VCF RG version | Reference alleles v37 and v38 are identical | All alleles are reference with v37 or v38 | Variant Annotation Non-VCF RG Version | Filter Non-VCF RG Version
---|---|---|---|---|---|---
False|*|*|*|*|Variant Annotation VCF RG Version + "?"|UNKNOWN
True|Inferred|True|*|*|REF_CALL|NO_CALL
True|Inferred|False|*|*|From "refSeqDifferenceAnnotations" field in panel JSON|INFERRED_PASS
True|Observed|True|*|*|REF_CALL|PASS
True|Observed|False|True|*|Variant Annotation VCF RG Version|PASS
True|Observed|False|False|True|From "refSeqDifferenceAnnotations" field in panel JSON|PASS
True|Observed|False|False|False|Variant Annotation VCF RG Version + "?"|UNKNOWN

Note that an asterisk does not indicate that every value is actually possible. 
For instance, calls that do not match any variants from the panel JSON can only be observed calls, not inferred calls.
Furthermore, if there is no matching variant in the panel JSON, then PEACH does not know what the correct non-VCF RG version reference allele is, 
so it would be unknown whether the reference alleles wrt v37 and v38 are identical. 

Let's call these calls with both v37 and v38 details *full calls*.

### Infer Haplotypes
For the haplotype calling PEACH only uses the v38 and general information from the full calls, not the v37 information.

The goal is to find the simplest combination of haplotypes that explains the called variants. 

Sometimes, more than one combination of haplotypes could explain the calls. 
As an example, consider the DPYD gene and two variants for that gene: c.1905+1G>A and c.1627A>G. 
Separately, these variants form the haplotypes *2A and *5, and the haplotype *2B consists of both of these variants together.
A haplotype can contain multiple variants if these variants have a tendency to be inherited together.
If each variant is called once and if all of these variants and haplotypes are included in the panel JSON, then, to take this combined inheritance into account, 
PEACH prefers to call the haplotype combination as *2B_HET/*1_HET and not as *2A_HET/*5_HET. 
If you want PEACH to call *2A_HET/*5_HET instead of *2B_HET/*1_HET in this situation, then simply don't include *2B in the panel JSON.

To make this more precise, define the *length* of a haplotype combination as the total number of non-wild-type
haplotypes in the combination, where homozygous haplotype calls are counted as 2. 
PEACH will always attempt to call the unique haplotype combination of minimum length that explain all of the variant calls.
If there are no haplotype combinations that explain all of the variant calls, 
or if there is more than one combination of the same minimum length, 
then the haplotype combination for that gene is called as "Unresolved Haplotype". 

The only valid haplotype combination of length 0 is the homozygous wild type haplotype,
valid haplotype combinations of length 1 always include precisely one heterozygous wild type haplotype call,
and valid haplotype combinations of length at least 2 do not contain any calls for the wild type haplotype.

Note that when at least one of the VCF calls overlaps with but is not identical to 
one of the variants in the panel JSON, then the haplotype combination "Unresolved Haplotype" will be called, 
because this variant will be an observed VCF call that is not part of any haplotypes in the panel JSON.

#### Haplotype Calling Algorithm 
Haplotypes are called for each gene separately. First, collect the full calls that correspond to that gene. 
Then, extract the alt alleles wrt v38 from these full calls, 
and count the number of times each combination of position (v38) and alt allele (v38) occurs.
Use recursive descent to determine all haplotype combinations that perfectly explain all of these variants.
If there are no such combinations, then no haplotype combination can be called for this gene. 
If such combinations do exist, then the next step is to determine the length of
each valid haplotype combination. Find the minimum length of the valid haplotype combinations, 
and select the haplotype combinations whose length is equal to this minimum. 
If precisely one such haplotype combination exists, then this combination will be called for this gene.
If more than one haplotype combination of minimum length exists, then no haplotype combination is called for this gene.

### Examples
The data in these examples will be the completely fictional.
The examples will focus on fairly "standard" situations, and they will exclude all information that is not necessary to understand these situations. 
For details on non-standard situations, see the more detailed subsections of the [Algorithm](#algorithm) section.

Suppose that the panel JSON contains the following variants and haplotypes for the fictional gene FAKE,
and that FAKE is the only gene in the panel JSON.

Rs Id|Reference Allele V37|Reference Allele V38|V38 Annotation for Reference Sequence Difference
---|---|---|---
rs1|A|A|N/A
rs2|TA|GC|c.6543GC>TA
rs3|GG|GG|N/A

Haplotype|Variants (Rs Id: Variant Allele wrt v38)
---|---
*1 (wild type)|None
*2|rs1: T
*3|rs2: TA
*4|rs3: G
*5|rs1: T, rs3: G

#### No Calls
If there are no calls wrt v37 in the VCF, then the full calls are:

Rs Id|Allele1|Allele2|Variant Annotation V37|Filter V37|Variant Annotation V38|Filter V38
---|---|---|---|---|---|---
rs1|A|A|REF_CALL|NO_CALL|REF_CALL|NO_CALL
rs2|TA|TA|REF_CALL|NO_CALL|c.6543GC>TA|INFERRED_PASS
rs3|GG|GG|REF_CALL|NO_CALL|REF_CALL|NO_CALL

The only valid haplotype combination that explains these variants is *3_HOM, 
so this is the haplotype combination that is called for FAKE.

#### Homozygous Wild Type
Suppose that the v37 calls are the following:

Rs Id|Allele1|Allele2|Variant Annotation V37|Filter V37
---|---|---|---|---
rs1|A|A|c.8483A>T|PASS
rs2|GC|GC|c.6543TA>GC|PASS
rs3|GG|GG|c.4838GG>G|PASS

In this case, the full calls are:

Rs Id|Allele1|Allele2|Variant Annotation V37|Filter V37|Variant Annotation V38|Filter V38
---|---|---|---|---|---|---
rs1|A|A|REF_CALL|PASS|REF_CALL|PASS
rs2|GC|GC|c.6543TA>GC|PASS|REF_CALL|PASS
rs3|GG|GG|REF_CALL|PASS|REF_CALL|PASS

The only valid haplotype combination is *1_HOM, so this haplotype combination is called for FAKE.

#### Heterozygous Wild Type
Suppose that the v37 calls are the following:

Rs Id|Allele1|Allele2|Variant Annotation V37|Filter V37
---|---|---|---|---
rs1|A|A|c.8483A>T|PASS
rs2|TA|GC|c.6543TA>GC|PASS
rs3|GG|GG|c.4838GG>G|PASS

The full calls are:

Rs Id|Allele1|Allele2|Variant Annotation V37|Filter V37|Variant Annotation V38|Filter V38
---|---|---|---|---|---|---
rs1|A|A|REF_CALL|PASS|REF_CALL|PASS
rs2|GC|TA|c.6543TA>GC|PASS|c.6543GC>TA|PASS
rs3|GG|GG|REF_CALL|PASS|REF_CALL|PASS

The only valid haplotype combination for FAKE is *3_HET/*1_HET, so this haplotype combination is called.

#### Multiple Valid Haplotypes
Suppose that the v37 calls are the following:

Rs Id|Allele1|Allele2|Variant Annotation V37|Filter V37
---|---|---|---|---
rs1|A|T|c.8483A>T|PASS
rs2|GC|GC|c.6543TA>GC|PASS
rs3|G|G|c.4838GG>G|PASS

The resulting full calls are:

Rs Id|Allele1|Allele2|Variant Annotation V37|Filter V37|Variant Annotation V38|Filter V38
---|---|---|---|---|---|---
rs1|A|T|c.8483A>T|PASS|c.8483A>T|PASS
rs2|GC|GC|c.6543TA>GC|PASS|REF_CALL|PASS
rs3|G|G|c.4838GG>G|PASS|c.4838GG>G|PASS

The valid haplotype combinations are *2_HET/*4_HOM and *4_HET/*5_HET. 
These combinations have lengths 3 and 2, respectively, so the second combination is preferred.
The called haplotype combination for FAKE is *4_HET/*5_HET.

### Restrictions
PEACH does not support calling for multiple (partially) overlapping genes.
If one wishes to attain results for (partially) overlapping genes anyway,
split them across separate panel JSON files and run PEACH multiple times.

Variants in a panel JSON file are not allowed to (partially) overlap.

## Tests
To run PEACH's test suite, including mypy, run the script `test_peach`.
If you have installed PEACH's requirements into a venv, then remember to source the venv before running `test_peach`.

## Version History and Download Links
* Upcoming
    + Nothing
* [1.1](https://github.com/hartwigmedical/peach/releases/tag/v1.1)
    + Add shell script `peach` for running PEACH.
    + Remove VCF filtering step. 
        + Remove VCFTools dependency.
        + Remove filtered vcf output file.
    + Add experimental support for input VCF's for reference genomes with version v38. 
    + Change format of arguments to PEACH. 
        + Arguments are no longer positional. 
        + Remove arguments vcftools, recreate_bed and transcript_tsv, since they are no longer needed.
        + Add optional parameter for experimental v38 VCF support.
    + Change format of panel JSON.
        + Change key "url_prescription_info" to "urlPrescriptionInfo" for consistency with other keys.
        + Add "annotationV37" key for reference sequence differences, for support of v38 reference genomes.
    + Adjust format of genotype TSV output file.
        + Split "haplotype" column into "haplotype" and "zygosity".
    + Add script for running tests.
* [1.0](https://github.com/hartwigmedical/peach/releases/tag/v1.0)
    + First release
