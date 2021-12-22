import argparse
from enum import Enum
from typing import List, Optional, Any

from base.reference_assembly import ReferenceAssembly
from tool_config import ToolConfig


class EnumAction(argparse.Action):
    """
    Argparse action for handling Enums
    """

    #  https://stackoverflow.com/questions/43968006/support-for-enum-arguments-in-argparse/55500795
    def __init__(self, **kwargs: Any) -> None:
        # Pop off the type value
        enum = kwargs.pop("type", None)

        # Ensure an Enum subclass is provided
        if enum is None:
            raise ValueError("type must be assigned an Enum when using EnumAction")
        if not issubclass(enum, Enum):
            raise TypeError("type must be an Enum when using EnumAction")

        # Generate choices from the Enum
        kwargs.setdefault("choices", tuple(e.name for e in enum))

        super(EnumAction, self).__init__(**kwargs)

        self._enum = enum

    def __call__(
        self,
        parser: argparse.ArgumentParser,
        namespace: argparse.Namespace,
        values: Any,
        option_string: Optional[str] = None,
    ) -> None:
        # Convert value back into an Enum
        enum = self._enum[values]
        setattr(namespace, self.dest, enum)


class ArgumentParser(object):
    def get_tool_config(self, sys_args: List[str]) -> ToolConfig:
        parser = argparse.ArgumentParser(
            prog="peach",
            description=(
                "Run pharmacogenomics panel on v37 germline VCF file. The pharmacogenomic annotations are done on "
                "v38, so output for both reference genomes is given where possible."
            ),
        )
        parser._action_groups.pop()
        required = parser.add_argument_group("required arguments")
        required.add_argument(
            "--vcf", "-i", type=str, required=True, help="VCF file to use for pharmacogenomics analysis."
        )
        required.add_argument("--panel", "-p", type=str, required=True, help="Json file with the panel variants.")
        required.add_argument(
            "--outputdir", "-o", type=str, required=True, help="Directory to store output of pharmacogenomic analysis."
        )
        required.add_argument("--sample_t_id", "-t", type=str, required=True, help="The sample ID of the tumor.")
        required.add_argument("--sample_r_id", "-r", type=str, required=True, help="The sample ID of the normal.")
        required.add_argument("--tool_version", "-v", type=str, required=True, help="The version of the tool.")
        experimental_optional = parser.add_argument_group("experimental optional arguments")
        experimental_optional.add_argument(
            "--vcf_reference_assembly_version",
            "-a",
            type=ReferenceAssembly,
            default=ReferenceAssembly.V37,
            required=False,
            action=EnumAction,
            help=(
                "The version of the reference assembly wrt which the vcf has been constructed. "
                "Support for V38 is experimental. Default is V37."
            ),
        )

        args = parser.parse_args(sys_args)

        config = ToolConfig(
            args.vcf,
            args.panel,
            args.outputdir,
            args.sample_t_id,
            args.sample_r_id,
            args.tool_version,
            args.vcf_reference_assembly_version,
        )
        return config
