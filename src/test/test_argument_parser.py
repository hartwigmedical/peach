import unittest
from copy import deepcopy
from io import StringIO
from unittest.mock import patch

from argument_parser import ArgumentParser
from base.reference_assembly import ReferenceAssembly
from config.tool_config import ToolConfig


class TestArgumentParser(unittest.TestCase):
    def test_representative_example(self) -> None:
        arguments = [
            "--vcf", "vcf_file", "-t", "tumor_sample_id", "-r", "ref_sample_id", "-v", "script_version",
            "-o", "output_directory", "-p", "panel_location", "-a", "V37",
        ]
        actual_config = ArgumentParser.get_tool_config(arguments)

        expected_config = ToolConfig(
            "vcf_file",
            "panel_location",
            "output_directory",
            "tumor_sample_id",
            "ref_sample_id",
            "script_version",
            ReferenceAssembly.V37,
        )
        self.assertEqual(expected_config, actual_config)

    def test_all_long_args(self) -> None:
        arguments = [
            "--vcf", "vcf_file", "--sample_t_id", "tumor_sample_id", "--sample_r_id", "ref_sample_id",
            "--tool_version", "script_version", "--outputdir", "output_directory", "--panel",
            "panel_location", "--vcf_reference_assembly_version", "V38",
        ]
        actual_config = ArgumentParser.get_tool_config(arguments)

        expected_config = ToolConfig(
            "vcf_file",
            "panel_location",
            "output_directory",
            "tumor_sample_id",
            "ref_sample_id",
            "script_version",
            ReferenceAssembly.V38,
        )
        self.assertEqual(expected_config, actual_config)

    def test_all_short_args(self) -> None:
        arguments = [
            "-i", "vcf_file", "-t", "tumor_sample_id", "-r", "ref_sample_id",
            "-v", "script_version", "-o", "output_directory", "-p", "panel_location", "-a", "V38",
        ]
        actual_config = ArgumentParser.get_tool_config(arguments)

        expected_config = ToolConfig(
            "vcf_file",
            "panel_location",
            "output_directory",
            "tumor_sample_id",
            "ref_sample_id",
            "script_version",
            ReferenceAssembly.V38,
        )
        self.assertEqual(expected_config, actual_config)

    def test_defaults(self) -> None:
        arguments = [
            "--vcf", "vcf_file", "-t", "tumor_sample_id", "-r", "ref_sample_id", "-v", "script_version",
            "-o", "output_directory", "-p", "panel_location",
        ]
        actual_config = ArgumentParser.get_tool_config(arguments)

        expected_config = ToolConfig(
            "vcf_file",
            "panel_location",
            "output_directory",
            "tumor_sample_id",
            "ref_sample_id",
            "script_version",
            ReferenceAssembly.V37,
        )
        self.assertEqual(expected_config, actual_config)

    def test_missing_required_arguments(self) -> None:
        minimum_arguments = [
            "-i vcf_file", "-t tumor_sample_id", "-r ref_sample_id", "-v script_version",
            "-o output_directory", "-p panel_location",
        ]
        ArgumentParser.get_tool_config(minimum_arguments)

        for i in range(len(minimum_arguments)):
            too_few_arguments = deepcopy(minimum_arguments)
            too_few_arguments.pop(i)
            with patch("sys.stderr", new_callable=StringIO):  # silence the help message in this test
                with self.subTest(arguments=too_few_arguments):
                    with self.assertRaises(SystemExit):
                        ArgumentParser.get_tool_config(too_few_arguments)


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
