from typing import NamedTuple

from base.reference_assembly import ReferenceAssembly


class ToolConfig(NamedTuple):
    vcf_path: str
    panel_path: str
    output_dir: str
    sample_t_id: str
    sample_r_id: str
    tool_version: str
    vcf_reference_assembly: ReferenceAssembly

    def get_calls_output_file_path(self) -> str:
        return f"{self.output_dir}/{self.sample_t_id}.peach.calls.tsv"

    def get_genotype_output_file_path(self) -> str:
        return f"{self.output_dir}/{self.sample_t_id}.peach.genotype.tsv"
