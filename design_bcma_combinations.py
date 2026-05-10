#!/usr/bin/env python3
"""Generate constrained BCMA scFv combination-mutation designs.

The design rules implemented here are the user-specified constraints:
- chain A residues 1-121 are VH; residues 137-244 are VL.
- heavy-chain designs contain 6 mutations: 2 in each of HCDR1/HCDR2/HCDR3.
- light-chain designs contain 5 mutations: at least 1 in each LCDR, using the provided candidate sites.
- candidate sites are restricted to the 17 provided mutations.

The output CSV deliberately leaves SWISS-MODEL/PISA result fields blank/pending because those
web tools require external job execution for custom uploaded structures in this environment.
"""
from __future__ import annotations

import csv
from collections import OrderedDict
from pathlib import Path

PDB = Path("relaxed_teclistamab_bcma.pdb")
OUT_CSV = Path("bcma_combination_designs.csv")
OUT_FASTA = Path("bcma_mutant_sequences.fasta")
OUT_LOG = Path("bcma_design_run.log")

AA3_TO_1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}

# User-supplied candidate mutation set, annotated into CDR buckets for validation.
CDR_BUCKET = {
    "S31D": "HCDR1", "F35Y": "HCDR1",
    "I58W": "HCDR2", "P63A": "HCDR2", "K66R": "HCDR2",
    "D101E": "HCDR3", "G102A": "HCDR3", "V104T": "HCDR3",
    "G159V": "LCDR1", "K166Q": "LCDR1", "S167G": "LCDR1",
    "S187N": "LCDR2", "D188M": "LCDR2",
    "S229T": "LCDR3", "S230N": "LCDR3", "H232T": "LCDR3", "V233L": "LCDR3",
}

# Twenty constrained combinations.  The highest-priority rows preserve experimentally supportive
# motifs where possible, particularly P63A/K66R/V104T plus G159V/D188M/H232T/S229T.
DESIGNS = [
    ("BCMA-C01", ["S31D", "F35Y", "P63A", "K66R", "D101E", "V104T"], ["G159V", "K166Q", "D188M", "S229T", "H232T"]),
    ("BCMA-C02", ["S31D", "F35Y", "P63A", "K66R", "G102A", "V104T"], ["G159V", "K166Q", "D188M", "S229T", "H232T"]),
    ("BCMA-C03", ["S31D", "F35Y", "I58W", "K66R", "D101E", "V104T"], ["G159V", "K166Q", "D188M", "S229T", "H232T"]),
    ("BCMA-C04", ["S31D", "F35Y", "I58W", "P63A", "D101E", "V104T"], ["G159V", "K166Q", "D188M", "S229T", "H232T"]),
    ("BCMA-C05", ["S31D", "F35Y", "P63A", "K66R", "D101E", "V104T"], ["G159V", "S167G", "D188M", "S229T", "H232T"]),
    ("BCMA-C06", ["S31D", "F35Y", "P63A", "K66R", "G102A", "V104T"], ["G159V", "S167G", "D188M", "S229T", "H232T"]),
    ("BCMA-C07", ["S31D", "F35Y", "I58W", "K66R", "G102A", "V104T"], ["G159V", "K166Q", "D188M", "S229T", "H232T"]),
    ("BCMA-C08", ["S31D", "F35Y", "I58W", "P63A", "G102A", "V104T"], ["G159V", "K166Q", "D188M", "S229T", "H232T"]),
    ("BCMA-C09", ["S31D", "F35Y", "P63A", "K66R", "D101E", "V104T"], ["G159V", "K166Q", "S187N", "S229T", "H232T"]),
    ("BCMA-C10", ["S31D", "F35Y", "P63A", "K66R", "G102A", "V104T"], ["G159V", "K166Q", "S187N", "S229T", "H232T"]),
    ("BCMA-C11", ["S31D", "F35Y", "I58W", "K66R", "D101E", "V104T"], ["G159V", "K166Q", "S187N", "S229T", "H232T"]),
    ("BCMA-C12", ["S31D", "F35Y", "I58W", "P63A", "D101E", "V104T"], ["G159V", "K166Q", "S187N", "S229T", "H232T"]),
    ("BCMA-C13", ["S31D", "F35Y", "P63A", "K66R", "D101E", "V104T"], ["G159V", "K166Q", "D188M", "S230N", "V233L"]),
    ("BCMA-C14", ["S31D", "F35Y", "P63A", "K66R", "G102A", "V104T"], ["G159V", "K166Q", "D188M", "S230N", "V233L"]),
    ("BCMA-C15", ["S31D", "F35Y", "I58W", "K66R", "D101E", "V104T"], ["G159V", "K166Q", "D188M", "S230N", "V233L"]),
    ("BCMA-C16", ["S31D", "F35Y", "I58W", "P63A", "D101E", "V104T"], ["G159V", "K166Q", "D188M", "S230N", "V233L"]),
    ("BCMA-C17", ["S31D", "F35Y", "P63A", "K66R", "D101E", "V104T"], ["G159V", "S167G", "S187N", "S229T", "H232T"]),
    ("BCMA-C18", ["S31D", "F35Y", "P63A", "K66R", "G102A", "V104T"], ["G159V", "S167G", "S187N", "S229T", "H232T"]),
    ("BCMA-C19", ["S31D", "F35Y", "I58W", "K66R", "G102A", "V104T"], ["G159V", "S167G", "D188M", "S230N", "V233L"]),
    ("BCMA-C20", ["S31D", "F35Y", "I58W", "P63A", "G102A", "V104T"], ["G159V", "S167G", "D188M", "S230N", "V233L"]),
]


def chain_a_sequence() -> str:
    residues = OrderedDict()
    with PDB.open() as handle:
        for line in handle:
            if not line.startswith("ATOM") or line[21] != "A" or line[12:16].strip() != "CA":
                continue
            resseq = int(line[22:26])
            residues[resseq] = AA3_TO_1[line[17:20].strip()]
    expected = list(range(1, 245))
    observed = list(residues)
    if observed != expected:
        raise ValueError(f"Chain A residue numbering is not contiguous 1-244: {observed[:3]}...{observed[-3:]}")
    return "".join(residues[i] for i in expected)


def apply_mutations(seq: str, muts: list[str]) -> str:
    chars = list(seq)
    for mut in muts:
        wt, new = mut[0], mut[-1]
        pos = int(mut[1:-1])
        if chars[pos - 1] != wt:
            raise ValueError(f"{mut} does not match chain A WT residue {chars[pos - 1]} at {pos}")
        chars[pos - 1] = new
    return "".join(chars)


def validate_design(heavy: list[str], light: list[str]) -> tuple[bool, str]:
    buckets = {bucket: 0 for bucket in ["HCDR1", "HCDR2", "HCDR3", "LCDR1", "LCDR2", "LCDR3"]}
    for mut in heavy + light:
        if mut not in CDR_BUCKET:
            return False, f"non-candidate mutation {mut}"
        buckets[CDR_BUCKET[mut]] += 1
    checks = [
        len(heavy) == 6,
        5 <= len(light) <= 6,
        buckets["HCDR1"] >= 2,
        buckets["HCDR2"] >= 2,
        buckets["HCDR3"] >= 2,
        buckets["LCDR1"] >= 1,
        buckets["LCDR2"] >= 1,
        buckets["LCDR3"] >= 1,
    ]
    if not all(checks):
        return False, str(buckets)
    return True, "passes user constraints"


def adjacent_notes(muts: list[str]) -> str:
    positions = sorted(int(mut[1:-1]) for mut in muts)
    near = [(a, b) for a, b in zip(positions, positions[1:]) if b - a <= 2]
    return ";".join(f"{a}-{b}" for a, b in near) if near else "none"


def rationale(heavy: list[str], light: list[str]) -> str:
    notes = []
    if {"P63A", "K66R", "V104T"}.issubset(heavy):
        notes.append("retains experimentally favorable P63A/K66R/V104T-containing heavy motif")
    if {"G159V", "D188M", "H232T"}.issubset(light):
        notes.append("retains light-chain motif present in best observed multi-mutant trend")
    if "S229T" in light:
        notes.append("tests S229T only in combinations that also include supportive G159V/D188M or H232T context")
    if "K166Q" in light:
        notes.append("uses K166Q rather than experimentally detrimental K166R-containing patterns")
    return "; ".join(notes) or "diversifies CDR coverage while remaining in the 17-site candidate set"


def main() -> None:
    wt_seq = chain_a_sequence()
    rows = []
    fasta_entries = []
    log_lines = [
        "BCMA scFv combination-design generation log",
        f"working_directory={Path.cwd()}",
        f"template_pdb={PDB}",
        f"chain_a_length={len(wt_seq)}",
        "heavy_region=chain A residues 1-121",
        "light_region=chain A residues 137-244",
        "candidate_mutation_count=17",
    ]
    for design_id, heavy, light in DESIGNS:
        ok, message = validate_design(heavy, light)
        if not ok:
            raise ValueError(f"{design_id} failed validation: {message}")
        all_muts = heavy + light
        sequence = apply_mutations(wt_seq, all_muts)
        rows.append({
            "design_id": design_id,
            "heavy_mutations": "/".join(heavy),
            "light_mutations": "/".join(light),
            "combined_mutation": "/".join(all_muts) + ":WT",
            "n_heavy_mutations": len(heavy),
            "n_light_mutations": len(light),
            "constraint_check": message,
            "near_adjacent_pairs_pos_delta_le_2": adjacent_notes(all_muts),
            "experimental_rationale": rationale(heavy, light),
            "swiss_model_template": str(PDB),
            "swiss_model_status": "pending_external_online_run_custom_template_upload",
            "pisa_status": "pending_external_online_run_after_model_download",
            "interface_area_A2": "",
            "delta_iG_kcal_per_mol": "",
        })
        fasta_entries.append(f">{design_id} {'/'.join(all_muts)}\n{sequence}\n")
        log_lines.append(
            f"{design_id}: heavy={len(heavy)} light={len(light)} "
            f"constraint={message} near_adjacent={adjacent_notes(all_muts)}"
        )

    with OUT_CSV.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    OUT_FASTA.write_text("".join(fasta_entries))
    log_lines.extend([
        f"design_rows_written={len(rows)}",
        f"csv_output={OUT_CSV}",
        f"fasta_output={OUT_FASTA}",
        "swiss_model_status=pending_external_online_run_custom_template_upload",
        "pisa_status=pending_external_online_run_after_model_download",
        "interface_area_A2=blank_until_online_pisa_completion",
        "delta_iG_kcal_per_mol=blank_until_online_pisa_completion",
    ])
    OUT_LOG.write_text("\n".join(log_lines) + "\n")


if __name__ == "__main__":
    main()
