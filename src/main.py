import argparse
import json
import os
import sys
from typing import List

from base.reference_assembly import ReferenceAssembly
from config.panel import Panel
from pgx_analysis import PgxAnalyser, PgxAnalysis
from pgx_reporter import HaplotypeReporter, GenotypeReporter
from vcf_reader import VcfReader


def main(vcf: str, sample_t_id: str, sample_r_id: str, version: str,
         panel_path: str, outputdir: str, vcf_reference_assembly: ReferenceAssembly) -> None:
    """ Run pharmacogenomics analysis on sample """
    print("\n[INFO] ## START PHARMACOGENOMICS ANALYSIS")

    # Check if output dir exists, create if it does not
    if not os.path.exists(outputdir):
        try:
            os.makedirs(outputdir)
        except FileExistsError:
            # Directory already exists
            pass

    # Get configuration
    panel = load_panel(panel_path)

    if panel.is_empty():
        raise ValueError("No panel is given, so no analysis can be performed.")

    # Get data for patient
    vcf_call_data = VcfReader.get_call_data(vcf, panel, sample_r_id, vcf_reference_assembly)

    if vcf_reference_assembly != ReferenceAssembly.V37:
        raise NotImplementedError()

    # Compute output from input data
    pgx_analysis = PgxAnalyser.create_pgx_analysis(vcf_call_data, panel)

    # Output
    print_calls_to_file(pgx_analysis, outputdir, sample_t_id, panel.get_id(), version)
    print_genotypes_to_file(pgx_analysis, panel, outputdir, sample_t_id, panel.get_id(), version)

    # TODO: add genes CYP2D6, CYP3A4, CYP3A5

    print("[INFO] ## PHARMACOGENOMICS ANALYSIS FINISHED\n")


def load_panel(panel_path: str) -> Panel:
    """ Load manually annotated JSON panel file """
    try:
        with open(panel_path, "r+", encoding="utf-8") as json_file:
            data = json.load(json_file)
            return Panel.from_json(data)
    except IOError:
        raise FileNotFoundError(f"Panel file {panel_path} not found or cannot be opened.")


def print_calls_to_file(pgx_analysis: PgxAnalysis, outputdir: str, sample_t_id: str,
                        panel_id: str, version: str) -> None:
    calls_file = f"{outputdir}/{sample_t_id}.peach.calls.tsv"
    if os.path.exists(calls_file):
        raise IOError(f"Calls output file {calls_file} already exists. Exiting.")
    with open(calls_file, "w") as f:
        f.write(GenotypeReporter.get_calls_tsv_text(pgx_analysis, panel_id, version))
    if not os.path.exists(calls_file):
        raise FileNotFoundError(f"Failed to write calls output file {calls_file}")


def print_genotypes_to_file(pgx_analysis: PgxAnalysis, panel: Panel, outputdir: str, sample_t_id: str,
                            panel_id: str, version: str) -> None:
    genotype_file = f"{outputdir}/{sample_t_id}.peach.genotype.tsv"
    if os.path.exists(genotype_file):
        raise IOError(f"Genotype output file {genotype_file} already exists. Exiting.")
    with open(genotype_file, "w") as f:
        f.write(HaplotypeReporter.get_genotype_tsv_text(pgx_analysis, panel, panel_id, version))
    if not os.path.exists(genotype_file):
        raise FileNotFoundError(f"Failed to write calls output file {genotype_file}")


def parse_args(sys_args: List[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="peach",
        description=("Run pharmacogenomics panel on v37 germline VCF file. The pharmacogenomic annotations are done on "
                     "v38, so output for both reference genomes is given where possible.")
    )
    parser._action_groups.pop()
    required = parser.add_argument_group("required arguments")
    required.add_argument(
        "--vcf", "-i", type=str, help="VCF file to use for pharmacogenomics analysis.", required=True)
    required.add_argument(
        "--sample_r_id", "-r", type=str, help="The sample ID of the normal.", required=True)
    required.add_argument(
        "--sample_t_id", "-t", type=str, help="The sample ID of the tumor.", required=True)
    required.add_argument(
        "--panel", "-p", type=str, help="Json file with the panel variants.", required=True)
    required.add_argument(
        "--outputdir", "-o", type=str, help="Directory to store output of pharmacogenomic analysis.", required=True)
    required.add_argument(
        "--version", "-v", type=str, help="The version of the tool.", required=True)

    return parser.parse_args(sys_args)


if __name__ == "__main__":
    args = parse_args(sys.argv[1:])
    main(args.vcf, args.sample_t_id, args.sample_r_id, args.version, args.panel, args.outputdir, ReferenceAssembly.V37)
