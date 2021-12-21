import re
from typing import Dict, Union

import pandas as pd

from base.constants import UNKNOWN_FUNCTION_STRING
from base.reference_assembly import ReferenceAssembly
from call_data import HaplotypeCall
from config.panel import Panel
from analysis.pgx_analysis import PgxAnalysis


class GenotypeReporter(object):
    GENE_COLUMN_NAME = "gene"
    CHROMOSOME_V37_COLUMN_NAME = "chromosome"
    CHROMOSOME_V38_COLUMN_NAME = "chromosome_v38"
    POSITION_V37_COLUMN_NAME = "position_v37"
    POSITION_V38_COLUMN_NAME = "position_v38"
    REF_ALLELE_V37_COLUMN_NAME = "ref_v37"
    REF_ALLELE_V38_COLUMN_NAME = "ref_v38"
    FIRST_ALLELE_COLUMN_NAME = "allele1"
    SECOND_ALLELE_COLUMN_NAME = "allele2"
    RS_IDS_COLUMN_NAME = "rsid"
    ANNOTATION_V37_COLUMN_NAME = "variant_annotation_v37"
    ANNOTATION_V38_COLUMN_NAME = "variant_annotation_v38"
    FILTER_V37_COLUMN_NAME = "filter_v37"
    FILTER_V38_COLUMN_NAME = "filter_v38"
    PANEL_VERSION_COLUMN_NAME = "panel_version"
    TOOL_VERSION_COLUMN_NAME = "repo_version"
    CALLS_TSV_COLUMNS = (
        GENE_COLUMN_NAME,
        CHROMOSOME_V37_COLUMN_NAME,
        POSITION_V37_COLUMN_NAME,
        POSITION_V38_COLUMN_NAME,
        REF_ALLELE_V37_COLUMN_NAME,
        REF_ALLELE_V38_COLUMN_NAME,
        FIRST_ALLELE_COLUMN_NAME,
        SECOND_ALLELE_COLUMN_NAME,
        RS_IDS_COLUMN_NAME,
        ANNOTATION_V37_COLUMN_NAME,
        FILTER_V37_COLUMN_NAME,
        ANNOTATION_V38_COLUMN_NAME,
        FILTER_V38_COLUMN_NAME,
        PANEL_VERSION_COLUMN_NAME,
        TOOL_VERSION_COLUMN_NAME,
        CHROMOSOME_V38_COLUMN_NAME,
    )
    CHROMOSOME_INDEX_NAME = "chromosome_index"

    UNKNOWN_CHROMOSOME_STRING = "UNKNOWN"
    UNKNOWN_POSITION_STRING = "UNKNOWN"
    UNKNOWN_REF_ALLELE_STRING = "UNKNOWN"

    TSV_SEPARATOR = "\t"
    RS_ID_SEPARATOR = ";"

    @classmethod
    def get_calls_tsv_text(
        cls, pgx_analysis: PgxAnalysis, panel_id: str, version: str, input_reference_assembly: ReferenceAssembly
    ) -> str:
        panel_calls_df = cls.__get_panel_calls_df(pgx_analysis, panel_id, version, input_reference_assembly)
        return str(panel_calls_df.to_csv(sep=cls.TSV_SEPARATOR, index=False))

    @classmethod
    def __get_panel_calls_df(
        cls, pgx_analysis: PgxAnalysis, panel_id: str, version: str, input_reference_assembly: ReferenceAssembly
    ) -> pd.DataFrame:
        data_frame = pd.DataFrame(columns=cls.CALLS_TSV_COLUMNS)
        for full_call in pgx_analysis.get_all_full_calls():
            sorted_alleles = sorted(full_call.alleles)

            position_v37: Union[str, int]
            position_v38: Union[str, int]
            if full_call.reference_site_v37 is not None:
                chromosome_v37 = full_call.reference_site_v37.start_coordinate.chromosome
                position_v37 = full_call.reference_site_v37.start_coordinate.position
                reference_allele_v37 = full_call.reference_site_v37.allele
            else:
                chromosome_v37 = cls.UNKNOWN_CHROMOSOME_STRING
                position_v37 = cls.UNKNOWN_POSITION_STRING
                reference_allele_v37 = cls.UNKNOWN_REF_ALLELE_STRING
            if full_call.reference_site_v38 is not None:
                chromosome_v38 = full_call.reference_site_v38.start_coordinate.chromosome
                position_v38 = full_call.reference_site_v38.start_coordinate.position
                reference_allele_v38 = full_call.reference_site_v38.allele
            else:
                chromosome_v38 = cls.UNKNOWN_CHROMOSOME_STRING
                position_v38 = cls.UNKNOWN_POSITION_STRING
                reference_allele_v38 = cls.UNKNOWN_REF_ALLELE_STRING

            line_data: Dict[str, Union[str, int]] = {
                cls.GENE_COLUMN_NAME: full_call.gene,
                cls.CHROMOSOME_V37_COLUMN_NAME: chromosome_v37,
                cls.POSITION_V37_COLUMN_NAME: position_v37,
                cls.CHROMOSOME_V38_COLUMN_NAME: chromosome_v38,
                cls.POSITION_V38_COLUMN_NAME: position_v38,
                cls.REF_ALLELE_V37_COLUMN_NAME: reference_allele_v37,
                cls.REF_ALLELE_V38_COLUMN_NAME: reference_allele_v38,
                cls.FIRST_ALLELE_COLUMN_NAME: sorted_alleles[0],
                cls.SECOND_ALLELE_COLUMN_NAME: sorted_alleles[1],
                cls.RS_IDS_COLUMN_NAME: cls.RS_ID_SEPARATOR.join(full_call.rs_ids),
                cls.ANNOTATION_V37_COLUMN_NAME: full_call.variant_annotation_v37,
                cls.FILTER_V37_COLUMN_NAME: full_call.filter_v37.name,
                cls.ANNOTATION_V38_COLUMN_NAME: full_call.variant_annotation_v38,
                cls.FILTER_V38_COLUMN_NAME: full_call.filter_v38.name,
                cls.PANEL_VERSION_COLUMN_NAME: panel_id,
                cls.TOOL_VERSION_COLUMN_NAME: version,
            }
            data_frame = data_frame.append(line_data, ignore_index=True)

        if pd.isna(data_frame).any(axis=None):
            raise ValueError(f"This should never happen: Unhandled NaN values:\n{data_frame}")

        data_frame = cls.__sort_call_data_frame(data_frame, input_reference_assembly)

        return data_frame

    @classmethod
    def __sort_call_data_frame(
        cls, data_frame: pd.DataFrame, input_reference_assembly: ReferenceAssembly
    ) -> pd.DataFrame:
        if input_reference_assembly == ReferenceAssembly.V37:
            column_for_chromosome_sorting = cls.CHROMOSOME_V37_COLUMN_NAME
            column_for_position_sorting = cls.POSITION_V37_COLUMN_NAME
        elif input_reference_assembly == ReferenceAssembly.V38:
            column_for_chromosome_sorting = cls.CHROMOSOME_V38_COLUMN_NAME
            column_for_position_sorting = cls.POSITION_V38_COLUMN_NAME
        else:
            error_msg = f"Unrecognized reference assembly: {input_reference_assembly}"
            raise ValueError(error_msg)

        data_frame.loc[:, cls.CHROMOSOME_INDEX_NAME] = data_frame.loc[:, column_for_chromosome_sorting].apply(
            lambda x: tuple(int(y) if y.isnumeric() else y for y in re.split(r"(\d+)", x))
        )

        data_frame = data_frame.sort_values(
            by=[
                cls.CHROMOSOME_INDEX_NAME,
                column_for_position_sorting,
                cls.GENE_COLUMN_NAME,
                cls.REF_ALLELE_V37_COLUMN_NAME,
                cls.REF_ALLELE_V38_COLUMN_NAME,
            ]
        ).reset_index(drop=True)
        data_frame = data_frame.loc[:, cls.CALLS_TSV_COLUMNS]
        return data_frame


class HaplotypeReporter(object):
    GENOTYPE_TSV_COLUMNS = (
        "gene",
        "haplotype",
        "function",
        "linked_drugs",
        "url_prescription_info",
        "panel_version",
        "repo_version",
        "haplotype_only",
        "zygosity_only",
    )

    HAPLOTYPE_HOMOZYGOUS_ZYGOSITY = "HOM"
    HAPLOTYPE_HETEROZYGOUS_ZYGOSITY = "HET"

    UNRESOLVED_HAPLOTYPE_STRING = "Unresolved Haplotype"
    NOT_APPLICABLE_ZYGOSITY_STRING = "N/A"

    TSV_SEPARATOR = "\t"
    DRUG_SEPARATOR = ";"

    @classmethod
    def get_genotype_tsv_text(cls, pgx_analysis: PgxAnalysis, panel: Panel, version: str) -> str:
        gene_to_haplotype_calls = pgx_analysis.get_gene_to_haplotype_calls()

        genes_in_analysis = set(gene_to_haplotype_calls.keys())
        assert genes_in_analysis == panel.get_genes(), (
            f"Gene lists inconsistent.\n"
            f"From analysis={sorted(list(genes_in_analysis))}\n"
            f"From panel={sorted(list(panel.get_genes()))}"
        )

        gene_to_drug_info = {}
        for gene_info in panel.get_gene_infos():
            sorted_drugs = sorted(
                [drug for drug in gene_info.drugs], key=lambda info: (info.name, info.url_prescription_info)
            )
            gene_to_drug_info[gene_info.gene] = (
                cls.DRUG_SEPARATOR.join([drug.name for drug in sorted_drugs]),
                cls.DRUG_SEPARATOR.join([drug.url_prescription_info for drug in sorted_drugs]),
            )

        header = cls.TSV_SEPARATOR.join(cls.GENOTYPE_TSV_COLUMNS)
        lines = [header]
        for gene in sorted(gene_to_haplotype_calls.keys()):
            if gene_to_haplotype_calls[gene]:
                for haplotype_call in sorted(gene_to_haplotype_calls[gene], key=lambda call: call.haplotype_name):
                    line_contents = [
                        gene,
                        f"{haplotype_call.haplotype_name}_{cls.__get_zygosity(haplotype_call)}",
                        panel.get_haplotype_function(gene, haplotype_call.haplotype_name),
                        gene_to_drug_info[gene][0],
                        gene_to_drug_info[gene][1],
                        panel.get_id(),
                        version,
                        haplotype_call.haplotype_name,
                        cls.__get_zygosity(haplotype_call),
                    ]
                    lines.append(cls.TSV_SEPARATOR.join(line_contents))
            else:
                line_contents = [
                    gene,
                    cls.UNRESOLVED_HAPLOTYPE_STRING,
                    UNKNOWN_FUNCTION_STRING,
                    gene_to_drug_info[gene][0],
                    gene_to_drug_info[gene][1],
                    panel.get_id(),
                    version,
                    cls.UNRESOLVED_HAPLOTYPE_STRING,
                    cls.NOT_APPLICABLE_ZYGOSITY_STRING,
                ]
                lines.append(cls.TSV_SEPARATOR.join(line_contents))
        text = "\n".join(lines) + "\n"
        return text

    @classmethod
    def __get_zygosity(cls, haplotype_call: HaplotypeCall) -> str:
        if haplotype_call.count == 2:
            return cls.HAPLOTYPE_HOMOZYGOUS_ZYGOSITY
        elif haplotype_call.count == 1:
            return cls.HAPLOTYPE_HETEROZYGOUS_ZYGOSITY
        else:
            error_msg = f"Invalid haplotype count: haplotype call={haplotype_call}"
            raise ValueError(error_msg)
