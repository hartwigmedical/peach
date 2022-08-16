from typing import NamedTuple, Optional

from util.reference_assembly import ReferenceAssembly


class ToolConfig(NamedTuple):
    vcf_path: str
    panel_path: str
    output_dir: str
    sample_t_id: Optional[str]
    sample_r_id: str
    tool_version: str
    vcf_reference_assembly: ReferenceAssembly

    def get_calls_output_file_path(self) -> str:
        return f"{self.output_dir}/{self.get_sample_name()}.peach.calls.tsv"

    def get_genotype_output_file_path(self) -> str:
        return f"{self.output_dir}/{self.get_sample_name()}.peach.genotype.tsv"

    def get_sample_name(self) -> str:
        if self.sample_t_id is None:
            return self.sample_r_id
        else:
            return self.sample_t_id
