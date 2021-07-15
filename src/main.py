import json
import logging
import os
import sys

from argument_parser import ArgumentParser
from config.tool_config import ToolConfig
from config.panel import Panel
from pgx_analysis import PgxAnalyser, PgxAnalysis
from pgx_reporter import HaplotypeReporter, GenotypeReporter
from vcf_reader import VcfReader


def main(tool_config: ToolConfig) -> None:
    """ Run pharmacogenomics analysis on sample """
    set_up_logging()

    logging.info(f"PEACH STARTING FOR {tool_config.sample_t_id}")

    # Check if output dir exists, create if it does not
    if not os.path.exists(tool_config.output_dir):
        try:
            os.makedirs(tool_config.output_dir)
        except FileExistsError:
            # Directory already exists
            pass

    # Get configuration
    logging.info("Creating panel config from JSON")
    panel = load_panel(tool_config.panel_path)

    if panel.is_empty():
        raise ValueError("No panel is given, so no analysis can be performed.")

    # Get data for patient
    logging.info("Reading call data from VCF")
    vcf_call_data = VcfReader.get_call_data(tool_config, panel)

    # Compute output from input data
    pgx_analysis = PgxAnalyser.create_pgx_analysis(vcf_call_data, panel)

    # Output
    logging.info("Creating output files")
    print_calls_to_file(pgx_analysis, tool_config, panel.get_id())
    print_genotypes_to_file(pgx_analysis, panel, tool_config)

    # TODO: add genes CYP2D6, CYP3A4, CYP3A5

    logging.info(f"PEACH FINISHED FOR {tool_config.sample_t_id}")


def set_up_logging() -> None:
    logging.basicConfig(
        format='%(asctime)s - [%(levelname)-8s] - %(message)s',
        level=logging.INFO,
        datefmt='%Y-%m-%d %H:%M:%S')


def load_panel(panel_path: str) -> Panel:
    """ Load manually annotated JSON panel file """
    try:
        with open(panel_path, "r+", encoding="utf-8") as json_file:
            return Panel.from_json(json.load(json_file))
    except IOError:
        raise FileNotFoundError(f"Panel file {panel_path} not found or cannot be opened.")


def print_calls_to_file(
        pgx_analysis: PgxAnalysis, tool_config: ToolConfig, panel_id: str
) -> None:
    calls_file = tool_config.get_calls_output_file_path()
    if os.path.exists(calls_file):
        raise IOError(f"Calls output file {calls_file} already exists. Exiting.")
    with open(calls_file, "w") as f:
        text = GenotypeReporter.get_calls_tsv_text(
            pgx_analysis, panel_id, tool_config.tool_version, tool_config.vcf_reference_assembly)
        f.write(text)
    if not os.path.exists(calls_file):
        raise FileNotFoundError(f"Failed to write calls output file {calls_file}")


def print_genotypes_to_file(pgx_analysis: PgxAnalysis, panel: Panel, tool_config: ToolConfig) -> None:
    genotype_file = tool_config.get_genotype_output_file_path()
    if os.path.exists(genotype_file):
        raise IOError(f"Genotype output file {genotype_file} already exists. Exiting.")
    with open(genotype_file, "w") as f:
        f.write(HaplotypeReporter.get_genotype_tsv_text(pgx_analysis, panel, tool_config.tool_version))
    if not os.path.exists(genotype_file):
        raise FileNotFoundError(f"Failed to write calls output file {genotype_file}")


if __name__ == "__main__":
    config = ArgumentParser.get_tool_config(sys.argv[1:])
    main(config)
