import argparse
from typing import List

from base.reference_assembly import ReferenceAssembly
from config.tool_config import ToolConfig


class ArgumentParser(object):
    @classmethod
    def get_tool_config(cls, sys_args: List[str]) -> ToolConfig:
        parser = argparse.ArgumentParser(
            prog="peach",
            description=(
                "Run pharmacogenomics panel on v37 germline VCF file. The pharmacogenomic annotations are done on "
                "v38, so output for both reference genomes is given where possible."
            )
        )
        parser._action_groups.pop()
        required = parser.add_argument_group("required arguments")
        required.add_argument(
            "--vcf", "-i", type=str, required=True, help="VCF file to use for pharmacogenomics analysis.")
        required.add_argument(
            "--panel", "-p", type=str, required=True, help="Json file with the panel variants.")
        required.add_argument(
            "--outputdir", "-o", type=str, required=True, help="Directory to store output of pharmacogenomic analysis.")
        required.add_argument(
            "--sample_t_id", "-t", type=str, required=True, help="The sample ID of the tumor.")
        required.add_argument(
            "--sample_r_id", "-r", type=str, required=True, help="The sample ID of the normal.")
        required.add_argument(
            "--version", "-v", type=str, required=True, help="The version of the tool.")

        args = parser.parse_args(sys_args)

        config = ToolConfig(
            args.vcf,
            args.panel,
            args.outputdir,
            args.sample_t_id,
            args.sample_r_id,
            args.version,
            ReferenceAssembly.V37,
        )
        return config
