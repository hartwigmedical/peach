from typing import Set

from base.gene_coordinate import GeneCoordinate
from base.reference_site import ReferenceSite
from panel.annotation import Annotation
from panel.drug_info import DrugInfo
from panel.gene_info import GeneInfo
from panel.haplotype import Haplotype
from panel.panel import Panel
from panel.rs_id_info import RsIdInfo
from panel.variant import Variant


def get_wide_example_panel(include_transcript_ids: bool) -> Panel:
    dpyd_two_a_variant = Variant("rs3918290", "T")
    dpyd_two_b_variant = Variant("rs1801159", "C")
    dpyd_three_variant = Variant("rs72549303", "TG")
    fake_variant = Variant("rs1212125", "C")
    fake2_variant = Variant("rs1212127", "C")

    dpyd_haplotypes = frozenset({
        Haplotype("*2A", "No Function", frozenset({dpyd_two_a_variant})),
        Haplotype("*2B", "No Function", frozenset({dpyd_two_a_variant, dpyd_two_b_variant})),
        Haplotype("*3", "Normal Function", frozenset({dpyd_three_variant})),
    })
    dpyd_rs_id_infos = frozenset({
        RsIdInfo(
            "rs3918290", ReferenceSite(GeneCoordinate("1", 97915614), "C"),
            ReferenceSite(GeneCoordinate("chr1", 97450058), "C"), None),
        RsIdInfo(
            "rs72549309", ReferenceSite(GeneCoordinate("1", 98205966), "GATGA"),
            ReferenceSite(GeneCoordinate("chr1", 97740410), "GATGA"), None),
        RsIdInfo(
            "rs1801159", ReferenceSite(GeneCoordinate("1", 97981395), "T"),
            ReferenceSite(GeneCoordinate("chr1", 97515839), "T"), None),
        RsIdInfo(
            "rs72549303", ReferenceSite(GeneCoordinate("1", 97915621), "TG"),
            ReferenceSite(GeneCoordinate("chr1", 97450065), "TC"), Annotation("6744CA>GA", "6744GA>CA")),
    })
    dpyd_drugs = frozenset({
        DrugInfo("5-Fluorouracil", "https://www.pharmgkb.org/chemical/PA128406956/guidelineAnnotation/PA166104939"),
        DrugInfo("Capecitabine", "https://www.pharmgkb.org/chemical/PA448771/guidelineAnnotation/PA166104963"),
    })

    fake_haplotypes = frozenset({
        Haplotype("*4A", "Reduced Function", frozenset({fake_variant})),
    })
    fake_rs_id_infos = frozenset({
        RsIdInfo(
            "rs1212125", ReferenceSite(GeneCoordinate("5", 97915617), "T"),
            ReferenceSite(GeneCoordinate("chr5", 97450060), "T"), None),
    })
    fake_drugs = frozenset({
        DrugInfo("Aspirin", "https://www.pharmgkb.org/some_other_url"),
    })

    fake2_haplotypes = frozenset({
        Haplotype("*4A", "Reduced Function", frozenset({fake2_variant})),
    })
    fake2_rs_id_infos = frozenset({
        RsIdInfo(
            "rs1212127", ReferenceSite(GeneCoordinate("16", 97915617), "C"),
            ReferenceSite(GeneCoordinate("chr16", 97450060), "T"), Annotation("1324C>T", "1324T>C")),
    })
    fake2_drugs = frozenset({
        DrugInfo("Aspirin", "https://www.pharmgkb.org/some_other_url"),
    })

    dpyd_gene_info = GeneInfo(
        "DPYD", "*1", "ENST00000370192" if include_transcript_ids else None,
        dpyd_haplotypes, dpyd_rs_id_infos, dpyd_drugs,
    )
    fake_gene_info = GeneInfo(
        "FAKE", "*1", "ENST00000192883" if include_transcript_ids else None,
        fake_haplotypes, fake_rs_id_infos, fake_drugs,
    )
    fake2_gene_info = GeneInfo(
        "FAKE2", "*1", "ENST0000021918" if include_transcript_ids else None,
        fake2_haplotypes, fake2_rs_id_infos, fake2_drugs,
    )
    gene_infos = frozenset({dpyd_gene_info, fake_gene_info, fake2_gene_info})

    name = "WideTestPanel"
    version = "1.0"
    return Panel(name, version, gene_infos)


def get_narrow_example_panel(included_haplotypes: Set[str]) -> Panel:
    dpyd_two_a_variant = Variant("rs3918290", "T")
    dpyd_five_variant = Variant("rs1801159", "C")
    dpyd_three_variant = Variant("rs72549303", "TG")
    dpyd_seven_variant = Variant("rs2938101", "AGT")

    possible_dpyd_haplotypes = [
        Haplotype("*2A", "No Function", frozenset({dpyd_two_a_variant})),
        Haplotype("*2B", "No Function", frozenset({dpyd_two_a_variant, dpyd_five_variant})),
        Haplotype("*3", "Normal Function", frozenset({dpyd_three_variant})),
        Haplotype("*5", "Normal Function", frozenset({dpyd_five_variant})),
        Haplotype("*7", "Normal Function", frozenset({dpyd_seven_variant})),
        Haplotype("*9", "Normal Function", frozenset({dpyd_two_a_variant, dpyd_three_variant})),
        Haplotype("*10", "Normal Function", frozenset({dpyd_three_variant, dpyd_five_variant})),
    ]
    possible_rs_id_infos = [
        RsIdInfo(
            "rs3918290", ReferenceSite(GeneCoordinate("1", 97915614), "C"),
            ReferenceSite(GeneCoordinate("chr1", 97450058), "C"), None),
        RsIdInfo(
            "rs72549309", ReferenceSite(GeneCoordinate("1", 98205966), "GATGA"),
            ReferenceSite(GeneCoordinate("chr1", 97740410), "GATGA"), None),
        RsIdInfo(
            "rs1801159", ReferenceSite(GeneCoordinate("1", 97981395), "T"),
            ReferenceSite(GeneCoordinate("chr1", 97515839), "T"), None),
        RsIdInfo(
            "rs72549303", ReferenceSite(GeneCoordinate("1", 97915621), "TG"),
            ReferenceSite(GeneCoordinate("chr1", 97450065), "TC"), Annotation("6744CA>GA", "6744GA>CA")),
        RsIdInfo(
            "rs2938101", ReferenceSite(GeneCoordinate("1", 97912838), "A"),
            ReferenceSite(GeneCoordinate("chr1", 97453984), "A"), None),
    ]

    unknown_haplotypes = included_haplotypes.difference({haplotype.name for haplotype in possible_dpyd_haplotypes})
    if unknown_haplotypes:
        raise ValueError(f"Not all dpyd haplotype names recognized: unknown={unknown_haplotypes}")

    # use a set to cause the order of haplotypes to be random, to make sure this order will not matter.
    included_dpyd_haplotypes = frozenset(
        {haplotype for haplotype in possible_dpyd_haplotypes if haplotype.name in included_haplotypes}
    )
    included_dpyd_rs_ids = {"rs72549303"}.union(
        {variant.rs_id for haplotype in included_dpyd_haplotypes for variant in haplotype.variants}
    )
    included_dpyd_rs_id_infos = frozenset(
        {rs_id_info for rs_id_info in possible_rs_id_infos if rs_id_info.rs_id in included_dpyd_rs_ids}
    )
    dpyd_drugs = frozenset({
        DrugInfo("5-Fluorouracil", "https://www.pharmgkb.org/chemical/PA128406956/guidelineAnnotation/PA166104939"),
        DrugInfo("Capecitabine", "https://www.pharmgkb.org/chemical/PA448771/guidelineAnnotation/PA166104963"),
    })

    gene_infos = frozenset({
        GeneInfo(
            "DPYD", "*1", None, included_dpyd_haplotypes, included_dpyd_rs_id_infos,
            dpyd_drugs,
        ),
    })
    name = "NarrowTestPanel"
    version = "1.0"
    return Panel(name, version, gene_infos)


def get_empty_panel() -> Panel:
    return Panel("EmptyPanel", "0.3", frozenset())
