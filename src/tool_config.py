from pathlib import Path
from typing import NamedTuple, Optional

from util.reference_assembly import ReferenceAssembly


class ToolConfig(NamedTuple):
    vcf_path: Path
    panel_path: Path
    output_dir: Path
    sample_t_id: Optional[str]
    sample_r_id: str
    tool_version: str
    vcf_reference_assembly: ReferenceAssembly

    def get_calls_output_file_path(self) -> Path:
        return self.output_dir / f"{self.get_sample_name()}.peach.calls.tsv"

    def get_genotype_output_file_path(self) -> Path:
        return self.output_dir / f"{self.get_sample_name()}.peach.genotype.tsv"

    def get_sample_name(self) -> str:
        if self.sample_t_id is None:
            return self.sample_r_id
        else:
            return self.sample_t_id
