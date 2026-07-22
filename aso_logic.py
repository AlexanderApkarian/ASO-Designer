from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal


CHEMISTRY_PRESETS = {
    "Unmodified DNA": {
        "gap_length": 20,
        "wing_length": 0,
        "wing_chemistry": "None",
        "backbone_modification": "PO",
    },
    "5-10-5 MOE/DNA": {
        "gap_length": 10,
        "wing_length": 5,
        "wing_chemistry": "MOE",
        "backbone_modification": "PS",
    },
    "KT777/valeriasen": {
        "gap_length": 10,
        "wing_length": 5,
        "wing_chemistry": "MOE",
        "backbone_modification": "MIXED",
        "linkage_pattern": "KT777",
        "methyl_c": True,
    },
    "4-10-4 LNA/DNA": {
        "gap_length": 10,
        "wing_length": 4,
        "wing_chemistry": "LNA",
        "backbone_modification": "PS",
    },
    "3-12-3 LNA/DNA": {
        "gap_length": 12,
        "wing_length": 3,
        "wing_chemistry": "LNA",
        "backbone_modification": "PS",
    },
    "3-10-3 LNA/DNA": {
        "gap_length": 10,
        "wing_length": 3,
        "wing_chemistry": "LNA",
        "backbone_modification": "PS",
    },
    "3-9-3 LNA/DNA": {
        "gap_length": 9,
        "wing_length": 3,
        "wing_chemistry": "LNA",
        "backbone_modification": "PS",
    },
    "Custom": None,
}

DEFAULT_CHEMISTRY = "3-12-3 LNA/DNA"

RIBOSE_MODIFICATION_OPTIONS = (
    "DNA",
    "LNA",
    "MOE",
    "2'OMe",
    "2'F",
)

NUCLEOBASE_MODIFICATION_OPTIONS = ("None", "5MeC")

BASE_MODIFICATION_OPTIONS = (
    "DNA",
    "DNA + 5MeC",
    "LNA",
    "LNA + 5MeC",
    "MOE",
    "MOE + 5MeC",
    "2'OMe",
    "2'OMe + 5MeC",
    "2'F",
    "2'F + 5MeC",
)

LINKAGE_OPTIONS = ("PS", "PO")

PENALTY_POSITION_MODES = (
    "Central core positions",
    "All ASO positions",
    "Selected ASO positions",
)

PENALTY_BASE_MODES = (
    "Like-for-like only",
    "All mismatches with scoring",
    "Non-wobble only",
)

ASO_BASE_OPTIONS = ("A", "C", "G", "U")


DisplaySpanKind = Literal["mutation", "gap"]


@dataclass(frozen=True)
class DisplaySpan:
    start: int
    end: int
    kind: DisplaySpanKind


@dataclass(frozen=True)
class AsoInputs:
    target_gene: str = ""
    snp_identifier: str = ""
    chemistry_number: str = ""
    aso_chemistry: str = DEFAULT_CHEMISTRY
    gap_length: int = 12
    wing_length: int = 3
    wing_chemistry: str = "LNA"
    backbone_modification: str = "PS"
    custom_base_modifications: tuple[str, ...] = ()
    custom_linkages: tuple[str, ...] = ()
    microwalk_step_size: int = 1
    mutation_type: str = "Insertion"
    mutation_length: int = 0
    mutation_start: int = 1
    rna_sequence: str = ""


@dataclass(frozen=True)
class ChemistrySettings:
    label: str
    gap_length: int
    wing_length: int
    wing_chemistry: str
    backbone_modification: str
    base_modifications: tuple[str, ...]
    linkages: tuple[str, ...]
    methyl_c: bool = False
    linkage_pattern: str = ""

    @property
    def aso_length(self) -> int:
        return len(self.base_modifications)


@dataclass(frozen=True)
class AsoRow:
    row_number: int
    aso_id: str
    idt_code: str
    clean_sequence: str
    display_sequence: str
    display_spans: tuple[DisplaySpan, ...]
    starting_position: int
    chemistry: str
    grid_cells: tuple[str, ...]
    core_start: int
    core_end: int


@dataclass(frozen=True)
class AsoResult:
    inputs: AsoInputs
    chemistry: ChemistrySettings
    clean_reversed_rna: str
    aso_length: int
    mutation_start_reversed: int
    mutation_start_reversed_options: tuple[int, ...]
    variant_bases: int
    crop_start: int
    crop_end: int
    displayed_bases: int
    required_asos: int
    complete_required_asos: int
    coverage_warning: str
    ambiguity_warning: str
    grid_width_status: str
    header_positions: tuple[int, ...]
    header_bases: tuple[str, ...]
    rows: tuple[AsoRow, ...] = field(default_factory=tuple)


@dataclass(frozen=True)
class IdtConversionRow:
    row_number: int
    input_sequence: str
    clean_sequence: str
    idt_code: str


@dataclass(frozen=True)
class ChemistryOptimizationRow:
    row_number: int
    motif_aso_positions: tuple[int, ...]
    motif_gap_positions: tuple[int, ...]
    clean_sequence: str
    idt_code: str
    base_modifications: tuple[str, ...]
    linkages: tuple[str, ...]


@dataclass(frozen=True)
class PenaltyAsoInputs:
    target_gene: str = ""
    target_identifier: str = ""
    chemistry_number: str = ""
    aso_chemistry: str = DEFAULT_CHEMISTRY
    gap_length: int = 12
    wing_length: int = 3
    wing_chemistry: str = "LNA"
    backbone_modification: str = "PS"
    custom_base_modifications: tuple[str, ...] = ()
    custom_linkages: tuple[str, ...] = ()
    parent_start: int = 0
    parent_count: int = 5
    microwalk_step_size: int = 1
    penalty_position_mode: str = "Central core positions"
    selected_penalty_positions: str = ""
    penalty_base_mode: str = "Like-for-like only"
    rna_sequence: str = ""


@dataclass(frozen=True)
class PenaltyAsoRow:
    row_number: int
    parent_aso_id: str
    penalty_aso_id: str
    parent_start: int
    penalty_aso_position: int
    target_position_3to5: int
    target_position_5to3: int
    local_rna_context: str
    target_base: str
    canonical_aso_base: str
    penalty_aso_base: str
    mismatch_pair: str
    priority: str
    score: int
    reason: str
    clean_sequence: str
    idt_code: str
    chemistry: str
    grid_cells: tuple[str, ...]
    penalty_grid_index: int


@dataclass(frozen=True)
class PenaltyAsoResult:
    inputs: PenaltyAsoInputs
    chemistry: ChemistrySettings
    clean_reversed_rna: str
    aso_length: int
    crop_start: int
    crop_end: int
    header_positions: tuple[int, ...]
    header_bases: tuple[str, ...]
    parent_starts: tuple[int, ...]
    rows: tuple[PenaltyAsoRow, ...] = field(default_factory=tuple)


class AsoInputError(ValueError):
    pass


def chemistry_key(value: str) -> str:
    s = str(value).upper().strip().replace("\u00a0", "")
    s = s.replace(" ", "").replace("_", "-")
    aliases = {
        "5-10-5MOE/DNA": "5-10-5 MOE/DNA",
        "5-10-5MOE": "5-10-5 MOE/DNA",
        "KT777": "KT777/valeriasen",
        "VALERIASEN": "KT777/valeriasen",
        "KT777/VALERIASEN": "KT777/valeriasen",
        "VALERIASEN/KT777": "KT777/valeriasen",
        "4-10-4LNA/DNA": "4-10-4 LNA/DNA",
        "4-10-4LNA": "4-10-4 LNA/DNA",
        "3-12-3LNA/DNA": "3-12-3 LNA/DNA",
        "3-12-3LNA": "3-12-3 LNA/DNA",
        "3-10-3LNA/DNA": "3-10-3 LNA/DNA",
        "3-10-3LNA": "3-10-3 LNA/DNA",
        "3-9-3LNA/DNA": "3-9-3 LNA/DNA",
        "3-9-3LNA": "3-9-3 LNA/DNA",
        "UNMODIFIED": "Unmodified DNA",
        "UNMODIFIEDDNA": "Unmodified DNA",
        "PLAIN": "Unmodified DNA",
        "PLAINDNA": "Unmodified DNA",
        "CUSTOM": "Custom",
    }
    return aliases.get(s, "")


def _coerce_string_tuple(value: tuple[str, ...] | list[str] | str) -> tuple[str, ...]:
    if isinstance(value, str):
        if not value.strip():
            return ()
        return tuple(part.strip() for part in value.split(",") if part.strip())
    return tuple(str(part) for part in value)


def resolve_chemistry(inputs: AsoInputs) -> ChemistrySettings:
    label = chemistry_key(inputs.aso_chemistry)
    if not label:
        raise AsoInputError(f"ASO chemistry not recognised: {inputs.aso_chemistry}")

    preset = CHEMISTRY_PRESETS[label]
    if preset is not None:
        gap_length = int(preset["gap_length"])
        wing_length = int(preset["wing_length"])
        wing_chemistry = str(preset["wing_chemistry"])
        backbone_modification = str(preset["backbone_modification"])
        methyl_c = bool(preset.get("methyl_c", False))
        linkage_pattern = str(preset.get("linkage_pattern", ""))
        base_modifications, linkages = build_gapmer_pattern(
            gap_length,
            wing_length,
            wing_chemistry,
            backbone_modification,
            methyl_c=methyl_c,
            linkage_pattern=linkage_pattern,
        )
    else:
        gap_length = int(inputs.gap_length)
        wing_length = int(inputs.wing_length)
        wing_chemistry = inputs.wing_chemistry
        backbone_modification = inputs.backbone_modification
        methyl_c = False
        linkage_pattern = ""
        custom_base_modifications = _coerce_string_tuple(inputs.custom_base_modifications)
        custom_linkages = _coerce_string_tuple(inputs.custom_linkages)
        if custom_base_modifications:
            base_modifications = tuple(normalise_base_modification(item) for item in custom_base_modifications)
            if any(not item for item in base_modifications):
                raise AsoInputError("Custom base modifications include an unrecognised option.")
            linkages = tuple(normalise_linkage(item) for item in custom_linkages)
            if any(not item for item in linkages):
                raise AsoInputError("Custom linkages must be PS or PO.")
            if len(linkages) != len(base_modifications) - 1:
                raise AsoInputError("Custom chemistry must have one fewer linkage than base positions.")
            if not base_modifications:
                raise AsoInputError("Custom chemistry must contain at least one base.")
            backbone_modification = summarise_linkages(linkages)
            wing_chemistry = "Custom"
            gap_length = 0
            wing_length = 0
        else:
            base_modifications, linkages = build_gapmer_pattern(
                gap_length,
                wing_length,
                wing_chemistry,
                backbone_modification,
            )

    if gap_length < 0 or wing_length < 0:
        raise AsoInputError("Gap length and wing length must be non-negative.")

    return ChemistrySettings(
        label=label,
        gap_length=gap_length,
        wing_length=wing_length,
        wing_chemistry=normalise_wing_chemistry(wing_chemistry) or str(wing_chemistry),
        backbone_modification=summarise_linkages(linkages),
        base_modifications=base_modifications,
        linkages=linkages,
        methyl_c=methyl_c,
        linkage_pattern=linkage_pattern,
    )


def chemistry_display_label(chemistry: ChemistrySettings) -> str:
    backbone_label = f"{chemistry.backbone_modification} backbone modification"
    if chemistry.label == "KT777/valeriasen":
        return f"{chemistry.label}, 5-10-5 MOE/DNA, {backbone_label}, 5MeC C bases"
    if chemistry.label == "Custom":
        if chemistry.gap_length == 0 and chemistry.wing_length == 0:
            return f"Custom, {chemistry.aso_length}-base custom pattern, {backbone_label}"
        custom_label = (
            f"{chemistry.wing_length}-{chemistry.gap_length}-{chemistry.wing_length} "
            f"{chemistry.wing_chemistry}/DNA"
        )
        return f"Custom, {custom_label}, {backbone_label}"
    return f"{chemistry.label}, {backbone_label}"


def clean_rna_for_reverse(raw: str) -> str:
    s = str(raw).lower()
    for old in ("-", "_", " ", "\t", "\r", "\n"):
        s = s.replace(old, "")
    return s


def reversed_clean_rna(raw: str) -> str:
    return clean_rna_for_reverse(raw)[::-1]


def ambiguous_insertion_start_indexes(clean_forward_rna: str, start_index: int, insertion_length: int) -> tuple[int, ...]:
    if insertion_length <= 0:
        return (start_index,)

    start = start_index
    end = start_index + insertion_length
    if start < 0 or end > len(clean_forward_rna):
        return (start_index,)

    leftmost = start
    left_end = end
    while leftmost > 0 and clean_forward_rna[leftmost - 1] == clean_forward_rna[left_end - 1]:
        leftmost -= 1
        left_end -= 1

    rightmost = start
    right_end = end
    while right_end < len(clean_forward_rna) and clean_forward_rna[rightmost] == clean_forward_rna[right_end]:
        rightmost += 1
        right_end += 1

    return tuple(range(leftmost, rightmost + 1))


def _format_position_options_one_based(indexes: tuple[int, ...]) -> str:
    positions = tuple(index + 1 for index in indexes)
    if not positions:
        return ""
    if len(positions) == 1:
        return str(positions[0])
    if positions == tuple(range(positions[0], positions[-1] + 1)):
        return f"{positions[0]}-{positions[-1]}"
    if len(positions) <= 8:
        return ", ".join(str(position) for position in positions)
    return f"{positions[0]}-{positions[-1]} ({len(positions)} possible positions)"


def normalise_mutation_type(value: str) -> str:
    s = str(value).upper().strip().replace(" ", "").replace("-", "").replace("_", "")
    if s in {"INSERTION", "INSERT", "INS", "I"}:
        return "INSERTION"
    if s in {"DELETION", "DELETE", "DEL", "D"}:
        return "DELETION"
    if s in {"SUBSTITUTION", "SUBSTITUTE", "SUB", "S"}:
        return "SUBSTITUTION"
    return ""


def normalise_wing_chemistry(value: str) -> str:
    s = str(value).upper().strip().replace(" ", "").replace("-", "").replace("_", "")
    if s in {"LNA", "AFFINITYPLUS"}:
        return "LNA"
    if s in {"MOE", "2MOE", "2OMOE", "2OMETHOXYETHYL"}:
        return "MOE"
    if s in {"2OM", "2OME", "OM", "OME", "2OMETHYL", "2OMERNA"}:
        return "2'OMe"
    if s in {"2F", "F", "FLUORO", "2FLUORO"}:
        return "2'F"
    if s in {"DNA5MEC", "DNA+5MEC", "5MEC", "5MEDC", "METHYLC"}:
        return "DNA + 5MeC"
    if s == "DNA":
        return "DNA"
    if s in {"NONE", "NO", "NOMOD", "UNMODIFIED", ""}:
        return "NONE"
    return ""


def normalise_backbone(value: str) -> str:
    s = str(value).upper().strip().replace(" ", "").replace("-", "").replace("_", "")
    if s in {"PS", "FULLPS", "PHOSPHOROTHIOATE"}:
        return "PS"
    if s in {"PO", "NONE", "NO", "NOPS", "PHOSPHODIESTER", ""}:
        return "PO"
    if s in {"MIXED", "PSPO", "POPS", "MIXEDPSPO"}:
        return "MIXED PS/PO"
    return ""


def normalise_ribose_modification(value: str) -> str:
    s = str(value).upper().strip()
    s = s.replace("\u2019", "'").replace("\u2032", "'")
    compact = s.replace(" ", "").replace("-", "").replace("_", "").replace("/", "")
    if compact in {"DNA", "NONE", "NO", "NOMOD", "UNMODIFIED", ""}:
        return "DNA"
    if compact in {"LNA", "AFFINITYPLUS"}:
        return "LNA"
    if compact in {"MOE", "2MOE", "2OMOE", "2OMETHOXYETHYL"}:
        return "MOE"
    if compact in {"2'OME", "2OME", "OME", "OM", "2OMETHYL", "2OMERNA"}:
        return "2'OMe"
    if compact in {"2'F", "2F", "F", "FLUORO", "2FLUORO"}:
        return "2'F"
    return ""


def normalise_nucleobase_modification(value: str) -> str:
    s = str(value).upper().strip()
    compact = s.replace(" ", "").replace("-", "").replace("_", "").replace("/", "")
    if compact in {"NONE", "NO", "NOMOD", "UNMODIFIED", ""}:
        return "None"
    if compact in {"5MEC", "5MEDC", "MEC", "MEDC", "METHYLC", "METHYLDC"}:
        return "5MeC"
    return ""


def split_base_modification(value: str) -> tuple[str, str]:
    s = str(value).strip()
    compact = s.upper().replace("\u2019", "'").replace("\u2032", "'")
    compact = compact.replace(" ", "").replace("-", "").replace("_", "").replace("/", "")

    nucleobase_mod = "5MeC" if any(token in compact for token in ("5MEC", "5MEDC", "METHYLC", "METHYLDC")) else "None"

    ribose_mod = ""
    if "LNA" in compact or "AFFINITYPLUS" in compact:
        ribose_mod = "LNA"
    elif "MOE" in compact or "2OMETHOXYETHYL" in compact:
        ribose_mod = "MOE"
    elif "2'OME" in compact or "2OME" in compact or "2OMETHYL" in compact or "2OMERNA" in compact:
        ribose_mod = "2'OMe"
    elif "2'F" in compact or "2F" in compact or "FLUORO" in compact or compact == "F":
        ribose_mod = "2'F"
    elif "DNA" in compact or nucleobase_mod == "5MeC" or compact in {"NONE", "NO", "NOMOD", "UNMODIFIED", ""}:
        ribose_mod = "DNA"
    else:
        ribose_mod = normalise_ribose_modification(s)

    if not ribose_mod:
        return "", ""
    return ribose_mod, nucleobase_mod


def combine_base_modification(ribose_modification: str, nucleobase_modification: str = "None") -> str:
    ribose_mod = normalise_ribose_modification(ribose_modification)
    nucleobase_mod = normalise_nucleobase_modification(nucleobase_modification)
    if not ribose_mod or not nucleobase_mod:
        return ""
    if nucleobase_mod == "5MeC":
        return f"{ribose_mod} + 5MeC"
    return ribose_mod


def normalise_base_modification(value: str) -> str:
    ribose_mod, nucleobase_mod = split_base_modification(value)
    if not ribose_mod or not nucleobase_mod:
        return ""
    return combine_base_modification(ribose_mod, nucleobase_mod)


def normalise_linkage(value: str) -> str:
    s = str(value).upper().strip().replace(" ", "").replace("-", "").replace("_", "")
    if s in {"PS", "PHOSPHOROTHIOATE"}:
        return "PS"
    if s in {"PO", "NONE", "NO", "PHOSPHODIESTER", ""}:
        return "PO"
    return ""


def summarise_linkages(linkages: tuple[str, ...]) -> str:
    unique = set(linkages)
    if unique == {"PS"}:
        return "PS"
    if unique <= {"PO"}:
        return "PO"
    if unique <= {"PS", "PO"}:
        return "mixed PS/PO"
    return "mixed"


def build_gapmer_pattern(
    gap_length: int,
    wing_length: int,
    wing_chemistry: str,
    backbone_modification: str,
    *,
    methyl_c: bool = False,
    linkage_pattern: str = "",
) -> tuple[tuple[str, ...], tuple[str, ...]]:
    gap_len = int(gap_length)
    wing_len = int(wing_length)
    if gap_len < 0 or wing_len < 0:
        raise AsoInputError("Gap length and wing length must be non-negative.")

    wing_mod = normalise_base_modification(wing_chemistry)
    if wing_mod == "DNA" and normalise_wing_chemistry(wing_chemistry) == "NONE":
        wing_mod = "DNA"
    if not wing_mod:
        raise AsoInputError("Wing chemistry must be DNA, DNA + 5MeC, LNA, MOE, 2'OMe, 2'F, or None.")

    gap_mod = "DNA + 5MeC" if methyl_c else "DNA"
    base_modifications = tuple(
        [wing_mod] * wing_len + [gap_mod] * gap_len + [wing_mod] * wing_len
    )
    if not base_modifications:
        raise AsoInputError("ASO length must be greater than zero.")

    backbone_key = normalise_backbone(backbone_modification)
    if not backbone_key:
        raise AsoInputError("Backbone must be PS, PO, None, or mixed PS/PO.")

    if linkage_pattern.upper() == "KT777":
        linkages = ["PS"] * (len(base_modifications) - 1)
        for idx in range(1, max(wing_len - 1, 1)):
            if idx < len(linkages):
                linkages[idx] = "PO"
        right_start = len(base_modifications) - wing_len
        for idx in range(right_start, len(base_modifications) - 2):
            if 0 <= idx < len(linkages):
                linkages[idx] = "PO"
        return base_modifications, tuple(linkages)

    linkage = "PS" if backbone_key == "PS" else "PO"
    return base_modifications, tuple([linkage] * (len(base_modifications) - 1))


def validate_base_sequence(seq: str) -> None:
    bad = sorted({b for b in seq.upper() if b not in "ACGTU"})
    if bad:
        raise AsoInputError("RNA sequence must contain only A/C/G/T/U, whitespace, hyphens, or underscores.")


def complement_base(base: str) -> str:
    return {
        "A": "U",
        "U": "A",
        "T": "A",
        "C": "G",
        "G": "C",
    }.get(base.upper(), "")


def clean_seq_for_idt(seq: str) -> str:
    s = str(seq).strip()
    for old in (" ", "\t", "\r", "\n", "-", "##", "5'", "3'", "5?", "3?"):
        s = s.replace(old, "")
    return s.upper()


def dna_code(base: str) -> str:
    b = base.upper()
    if b in {"A", "C", "G"}:
        return b
    if b in {"T", "U"}:
        return "T"
    raise AsoInputError(f"Invalid DNA base: {base}")


def lna_code(base: str) -> str:
    return "+" + dna_code(base)


def rna_code(base: str) -> str:
    b = base.upper()
    if b in {"A", "C", "G"}:
        return b
    if b in {"T", "U"}:
        return "U"
    raise AsoInputError(f"Invalid RNA base: {base}")


def moe_code(base: str, position: int, total_length: int) -> str:
    if position == 1:
        prefix = "/52MOEr"
    elif position == total_length:
        prefix = "/32MOEr"
    else:
        prefix = "/i2MOEr"
    return prefix + dna_code(base) + "/"


def ome_code(base: str) -> str:
    return "m" + rna_code(base)


def fluoro_code(base: str, position: int, total_length: int) -> str:
    if position == 1:
        prefix = "/52F"
    elif position == total_length:
        prefix = "/32F"
    else:
        prefix = "/i2F"
    return prefix + rna_code(base) + "/"


def methyl_dna_c_code(base: str, position: int, total_length: int) -> str:
    if base.upper() != "C":
        return dna_code(base)
    if position == 1:
        return "/5Me-dC/"
    if position == total_length:
        return "/3Me-dC/"
    return "/iMe-dC/"


def modified_base_code(base: str, modification: str, position: int, total_length: int) -> str:
    mod = normalise_base_modification(modification)
    ribose_mod, nucleobase_mod = split_base_modification(mod)
    if not ribose_mod:
        raise AsoInputError(f"Unrecognised base modification: {modification}")
    if ribose_mod == "DNA" and nucleobase_mod == "5MeC":
        return methyl_dna_c_code(base, position, total_length)
    if ribose_mod == "DNA":
        return dna_code(base)
    if ribose_mod == "LNA":
        return lna_code(base)
    if ribose_mod == "MOE":
        return moe_code(base, position, total_length)
    if ribose_mod == "2'OMe":
        return ome_code(base)
    if ribose_mod == "2'F":
        return fluoro_code(base, position, total_length)
    raise AsoInputError(f"Unrecognised base modification: {modification}")


def idt_aso_per_position(
    raw_seq: str,
    base_modifications: tuple[str, ...] | list[str],
    linkages: tuple[str, ...] | list[str],
) -> str:
    try:
        seq = clean_seq_for_idt(raw_seq)
        if not seq:
            return ""

        if any(base not in "ACGTU" for base in seq.upper()):
            return "#ERROR: sequence must contain only A/C/G/T/U"

        base_mods = tuple(normalise_base_modification(mod) for mod in base_modifications)
        if any(not mod for mod in base_mods):
            return "#ERROR: unrecognised base modification"
        linkage_mods = tuple(normalise_linkage(linkage) for linkage in linkages)
        if any(not linkage for linkage in linkage_mods):
            return "#ERROR: linkages must be PS or PO"

        if len(seq) != len(base_mods):
            return f"#ERROR: sequence length is {len(seq)}, expected {len(base_mods)}"
        if len(linkage_mods) != max(0, len(base_mods) - 1):
            return "#ERROR: chemistry must have one fewer linkage than base positions"

        tokens = [
            modified_base_code(base, mod, position, len(seq))
            for position, (base, mod) in enumerate(zip(seq, base_mods), start=1)
        ]
        pieces = [tokens[0]]
        for linkage, token in zip(linkage_mods, tokens[1:]):
            pieces.append("*" if linkage == "PS" else "")
            pieces.append(token)
        return "".join(pieces)
    except Exception as exc:
        return f"#ERROR: {exc}"


def idt_aso_custom(
    raw_seq: str,
    gap_length: int,
    wing_length: int,
    wing_chemistry: str,
    backbone_modification: str,
) -> str:
    try:
        seq = clean_seq_for_idt(raw_seq)
        if not seq:
            return ""

        gap_len = int(gap_length)
        wing_len = int(wing_length)
        if gap_len < 0 or wing_len < 0:
            return "#ERROR: gap and wing lengths must be non-negative"

        expected_len = gap_len + 2 * wing_len
        if len(seq) != expected_len:
            return f"#ERROR: sequence length is {len(seq)}, expected {expected_len}"

        if any(base not in "ACGTU" for base in seq.upper()):
            return "#ERROR: sequence must contain only A/C/G/T/U"

        base_modifications, linkages = build_gapmer_pattern(
            gap_len,
            wing_len,
            wing_chemistry,
            backbone_modification,
        )
        return idt_aso_per_position(seq, base_modifications, linkages)
    except Exception as exc:
        return f"#ERROR: {exc}"


def convert_sequences_to_idt(
    raw_sequences: str | list[str] | tuple[str, ...],
    chemistry: ChemistrySettings,
) -> tuple[IdtConversionRow, ...]:
    if isinstance(raw_sequences, str):
        lines = raw_sequences.splitlines()
    else:
        lines = list(raw_sequences)

    rows: list[IdtConversionRow] = []
    for line in lines:
        input_sequence = str(line).strip()
        if not input_sequence:
            continue
        clean_sequence = clean_seq_for_idt(input_sequence)
        rows.append(
            IdtConversionRow(
                row_number=len(rows) + 1,
                input_sequence=input_sequence,
                clean_sequence=clean_sequence,
                idt_code=idt_aso_per_position(
                    clean_sequence,
                    chemistry.base_modifications,
                    chemistry.linkages,
                ),
            )
        )
    return tuple(rows)


def _most_common_modification(modifications: list[str]) -> str:
    if not modifications:
        return "DNA"
    counts: dict[str, int] = {}
    for modification in modifications:
        counts[modification] = counts.get(modification, 0) + 1
    return max(counts, key=lambda item: (counts[item], -modifications.index(item)))


def chemistry_optimization_walk(
    raw_seq: str,
    base_modifications: tuple[str, ...] | list[str],
    linkages: tuple[str, ...] | list[str],
    wing_length: int,
    gap_length: int,
    motif_positions: tuple[int, ...] | list[int],
    step_size: int = 1,
) -> tuple[ChemistryOptimizationRow, ...]:
    seq = clean_seq_for_idt(raw_seq)
    if not seq:
        raise AsoInputError("ASO sequence is empty after removing whitespace and hyphens.")
    if any(base not in "ACGTU" for base in seq.upper()):
        raise AsoInputError("ASO sequence must contain only A/C/G/T/U, whitespace, or hyphens.")

    base_mods = tuple(normalise_base_modification(mod) for mod in base_modifications)
    if any(not mod for mod in base_mods):
        raise AsoInputError("Chemistry pattern includes an unrecognised base modification.")
    linkage_mods = tuple(normalise_linkage(linkage) for linkage in linkages)
    if any(not linkage for linkage in linkage_mods):
        raise AsoInputError("Linkages must be PS or PO.")
    if len(seq) != len(base_mods):
        raise AsoInputError(f"ASO sequence length is {len(seq)}, expected {len(base_mods)} from the chemistry pattern.")
    if len(linkage_mods) != max(0, len(base_mods) - 1):
        raise AsoInputError("Chemistry pattern must have one fewer linkage than base positions.")

    motif = tuple(sorted(set(int(position) for position in motif_positions)))
    if not motif:
        raise AsoInputError("Select at least one base modification to iterate.")
    if motif[0] < 0 or motif[-1] >= len(base_mods):
        raise AsoInputError("Selected motif positions are outside the ASO chemistry pattern.")

    wing_len = max(0, int(wing_length))
    gap_len = max(0, int(gap_length))
    if gap_len > 0 and wing_len * 2 + gap_len <= len(base_mods):
        core_start = wing_len
        core_end = wing_len + gap_len - 1
    else:
        core_start = 0
        core_end = len(base_mods) - 1

    if any(position < core_start or position > core_end for position in motif):
        raise AsoInputError("Selected motif positions must be within the central core.")

    motif_start = motif[0]
    motif_width = motif[-1] - motif_start + 1
    if motif_width > core_end - core_start + 1:
        raise AsoInputError("Selected motif is wider than the central core.")
    step = int(step_size)
    if step < 1:
        raise AsoInputError("Step size must be at least 1.")

    motif_offsets = tuple(position - motif_start for position in motif)
    motif_modifications = tuple(base_mods[position] for position in motif)
    core_background_candidates = [
        base_mods[position]
        for position in range(core_start, core_end + 1)
        if position not in motif
    ]
    core_background = _most_common_modification(core_background_candidates)
    background_mods = list(base_mods)
    for position in motif:
        background_mods[position] = core_background

    rows: list[ChemistryOptimizationRow] = []
    for start in range(core_start, core_end - motif_width + 2, step):
        shifted_positions = tuple(start + offset for offset in motif_offsets)
        if shifted_positions[-1] > core_end:
            continue
        shifted_mods = list(background_mods)
        for position, modification in zip(shifted_positions, motif_modifications):
            shifted_mods[position] = modification
        rows.append(
            ChemistryOptimizationRow(
                row_number=len(rows) + 1,
                motif_aso_positions=tuple(position + 1 for position in shifted_positions),
                motif_gap_positions=tuple(position - core_start + 1 for position in shifted_positions),
                clean_sequence=seq,
                idt_code=idt_aso_per_position(seq, shifted_mods, linkage_mods),
                base_modifications=tuple(shifted_mods),
                linkages=linkage_mods,
            )
        )

    if not rows:
        raise AsoInputError("No chemistry walk variants could be generated.")
    return tuple(rows)


def normalise_penalty_position_mode(value: str) -> str:
    s = str(value).upper().strip().replace("_", " ").replace("-", " ")
    s = " ".join(s.split())
    if s in {"CENTRAL CORE POSITIONS", "CORE", "CORE POSITIONS", "CENTRAL CORE"}:
        return "Central core positions"
    if s in {"ALL ASO POSITIONS", "ALL POSITIONS", "ALL"}:
        return "All ASO positions"
    if s in {"SELECTED ASO POSITIONS", "SELECTED POSITIONS", "SELECTED"}:
        return "Selected ASO positions"
    return ""


def normalise_penalty_base_mode(value: str) -> str:
    s = str(value).upper().strip().replace("_", " ").replace("-", " ")
    s = " ".join(s.split())
    if s in {
        "LIKE FOR LIKE ONLY",
        "LIKE FOR LIKE",
        "LIKEFORLIKE",
        "SAME AS TARGET",
        "TARGET BASE",
        "SAME BASE",
        "DEFAULT",
    }:
        return "Like-for-like only"
    if s in {
        "ALL MISMATCHES WITH SCORING",
        "ALL MISMATCHES SCORED",
        "AUTO RANKED",
        "AUTO",
        "RANKED",
        "ALL MISMATCHES",
        "ALL",
        "ALL BASES",
    }:
        return "All mismatches with scoring"
    if s in {"NON WOBBLE ONLY", "NONWOBBLE ONLY", "NON WOBBLE", "NO WOBBLE"}:
        return "Non-wobble only"
    return ""


def parse_position_ranges(value: str, *, minimum: int, maximum: int) -> tuple[int, ...]:
    text = str(value).replace(";", ",")
    positions: set[int] = set()
    for raw_part in text.split(","):
        part = raw_part.strip()
        if not part:
            continue
        if "-" in part:
            start_text, end_text = part.split("-", 1)
            start = int(start_text.strip())
            end = int(end_text.strip())
            if end < start:
                start, end = end, start
            positions.update(range(start, end + 1))
        else:
            positions.add(int(part))
    if not positions:
        return ()
    if min(positions) < minimum or max(positions) > maximum:
        raise AsoInputError(f"Selected positions must be between {minimum} and {maximum}.")
    return tuple(sorted(positions))


def _display_rna_base(base: str) -> str:
    b = base.upper()
    return "U" if b == "T" else b


def _local_rna_context(clean_reversed_rna: str, target_index: int) -> str:
    left = _display_rna_base(clean_reversed_rna[target_index - 1]) if target_index > 0 else "-"
    base = _display_rna_base(clean_reversed_rna[target_index])
    right = _display_rna_base(clean_reversed_rna[target_index + 1]) if target_index + 1 < len(clean_reversed_rna) else "-"
    return f"{left}-{base}-{right}"


def _is_wobble_like(aso_base: str, rna_base: str) -> bool:
    aso = _display_rna_base(aso_base)
    rna = _display_rna_base(rna_base)
    return (aso == "G" and rna == "U") or (aso == "U" and rna == "G")


def _is_ga_outlier(aso_base: str, rna_base: str) -> bool:
    pair = {_display_rna_base(aso_base), _display_rna_base(rna_base)}
    return pair == {"G", "A"}


def _same_base_class(aso_base: str, rna_base: str) -> bool:
    purines = {"A", "G"}
    pyrimidines = {"C", "U"}
    aso = _display_rna_base(aso_base)
    rna = _display_rna_base(rna_base)
    return (aso in purines and rna in purines) or (aso in pyrimidines and rna in pyrimidines)


def score_penalty_candidate(
    *,
    target_base: str,
    penalty_base: str,
    left_base: str,
    right_base: str,
    aso_position_index: int,
    aso_length: int,
) -> tuple[str, int, str]:
    rna_base = _display_rna_base(target_base)
    aso_base = _display_rna_base(penalty_base)
    score = 0
    reasons: list[str] = []

    if _is_wobble_like(aso_base, rna_base):
        score -= 3
        reasons.append("G:U wobble-like mismatch may be less disruptive")
    else:
        score += 3
        reasons.append("non-wobble mismatch")

    if _same_base_class(aso_base, rna_base):
        score += 1
        reasons.append("same purine/pyrimidine class may increase geometric penalty")

    if _is_ga_outlier(aso_base, rna_base):
        score -= 1
        reasons.append("G:A/A:G mismatches can be context-dependent outliers")

    flank_bases = [_display_rna_base(base) for base in (left_base, right_base) if base and base != "-"]
    if len(flank_bases) < 2:
        score -= 2
        reasons.append("terminal/edge context; mismatch effect is less predictable")
    else:
        reasons.append("nearest-neighbor context recorded; exact penalty effect is mismatch-specific")

    if aso_position_index in {0, aso_length - 1}:
        score -= 2
        reasons.append("terminal ASO position")
    elif aso_position_index in {1, aso_length - 2}:
        score -= 1
        reasons.append("near-terminal ASO position")

    if score >= 4:
        priority = "Recommended"
    elif score >= 2:
        priority = "Alternative"
    else:
        priority = "Lower priority"
    return priority, score, "; ".join(reasons)


def _penalty_position_indexes(
    chemistry: ChemistrySettings,
    mode: str,
    selected_positions: str,
) -> tuple[int, ...]:
    aso_length = chemistry.aso_length
    if mode == "All ASO positions":
        return tuple(range(aso_length))
    if mode == "Selected ASO positions":
        selected = parse_position_ranges(selected_positions, minimum=1, maximum=aso_length)
        if not selected:
            raise AsoInputError("Enter at least one selected ASO position.")
        return tuple(position - 1 for position in selected)

    if chemistry.gap_length > 0 and chemistry.wing_length * 2 + chemistry.gap_length <= aso_length:
        start = chemistry.wing_length
        end = chemistry.wing_length + chemistry.gap_length
        return tuple(range(start, end))
    return tuple(range(aso_length))


def _penalty_candidate_bases(target_base: str) -> tuple[str, ...]:
    canonical = _display_rna_base(complement_base(target_base))
    return tuple(base for base in ASO_BASE_OPTIONS if base != canonical)


def _like_for_like_penalty_base(target_base: str) -> str:
    return _display_rna_base(target_base)


def generate_penalty_design(inputs: PenaltyAsoInputs) -> PenaltyAsoResult:
    chemistry = resolve_chemistry(
        AsoInputs(
            aso_chemistry=inputs.aso_chemistry,
            gap_length=inputs.gap_length,
            wing_length=inputs.wing_length,
            wing_chemistry=inputs.wing_chemistry,
            backbone_modification=inputs.backbone_modification,
            custom_base_modifications=inputs.custom_base_modifications,
            custom_linkages=inputs.custom_linkages,
        )
    )

    clean_reversed = reversed_clean_rna(inputs.rna_sequence)
    if not clean_reversed:
        raise AsoInputError("RNA sequence is empty after removing whitespace, hyphens, and underscores.")
    validate_base_sequence(clean_reversed)

    aso_length = chemistry.aso_length
    if aso_length <= 0:
        raise AsoInputError("ASO length must be greater than zero.")
    if len(clean_reversed) < aso_length:
        raise AsoInputError("RNA sequence is shorter than the selected ASO chemistry length.")

    try:
        parent_start = int(inputs.parent_start)
        parent_count = int(inputs.parent_count)
        step = int(inputs.microwalk_step_size)
    except Exception as exc:
        raise AsoInputError("Parent start, number of parent ASOs, and step size must be whole numbers.") from exc
    if parent_start < 0:
        raise AsoInputError("First parent ASO start position must be non-negative.")
    if parent_count < 1:
        raise AsoInputError("Number of parent ASOs must be at least 1.")
    if step < 1:
        raise AsoInputError("Step size must be at least 1.")

    max_start = len(clean_reversed) - aso_length
    if parent_start > max_start:
        raise AsoInputError("First parent ASO start position is outside the cleaned RNA sequence.")

    position_mode = normalise_penalty_position_mode(inputs.penalty_position_mode)
    if not position_mode:
        raise AsoInputError("Penalty position mode is not recognised.")
    base_mode = normalise_penalty_base_mode(inputs.penalty_base_mode)
    if not base_mode:
        raise AsoInputError("Penalty base mode is not recognised.")

    parent_starts = tuple(range(parent_start, max_start + 1, step))[:parent_count]
    if not parent_starts:
        raise AsoInputError("No parent ASOs could be generated.")

    penalty_indexes = _penalty_position_indexes(chemistry, position_mode, inputs.selected_penalty_positions)
    if not penalty_indexes:
        raise AsoInputError("No penalty positions could be generated.")

    crop_start = min(parent_starts)
    crop_end = max(start + aso_length - 1 for start in parent_starts)
    header_positions = tuple(range(crop_start, crop_end + 1))
    header_bases = tuple(_display_rna_base(clean_reversed[pos]).lower() for pos in header_positions)

    prefix_parts = [inputs.target_gene.strip(), inputs.target_identifier.strip(), inputs.chemistry_number.strip()]
    prefix = "_".join(part for part in prefix_parts if part) or "Penalty_ASO"
    chemistry_label = chemistry_display_label(chemistry)
    rows: list[PenaltyAsoRow] = []

    for parent_number, start in enumerate(parent_starts, start=1):
        target_segment = clean_reversed[start : start + aso_length]
        parent_sequence = "".join(_display_rna_base(complement_base(base)) for base in target_segment)
        parent_id = f"{prefix}_ASO_{parent_number}"
        for local_index in penalty_indexes:
            target_index = start + local_index
            target_base = _display_rna_base(clean_reversed[target_index])
            canonical_base = _display_rna_base(parent_sequence[local_index])
            left_base = _display_rna_base(clean_reversed[target_index - 1]) if target_index > 0 else "-"
            right_base = (
                _display_rna_base(clean_reversed[target_index + 1])
                if target_index + 1 < len(clean_reversed)
                else "-"
            )
            candidate_rows = []
            if base_mode == "Like-for-like only":
                candidate_bases = (_like_for_like_penalty_base(target_base),)
            else:
                candidate_bases = _penalty_candidate_bases(target_base)
            for penalty_base in candidate_bases:
                if base_mode == "Non-wobble only" and _is_wobble_like(penalty_base, target_base):
                    continue
                priority, score, reason = score_penalty_candidate(
                    target_base=target_base,
                    penalty_base=penalty_base,
                    left_base=left_base,
                    right_base=right_base,
                    aso_position_index=local_index,
                    aso_length=aso_length,
                )
                candidate_rows.append((score, priority, penalty_base, reason))
            candidate_rows.sort(key=lambda item: (-item[0], item[2]))

            for score, priority, penalty_base, reason in candidate_rows:
                clean_sequence_list = list(parent_sequence)
                clean_sequence_list[local_index] = penalty_base
                clean_sequence = "".join(clean_sequence_list)
                grid_cells = []
                for position in header_positions:
                    if start <= position <= start + aso_length - 1:
                        grid_cells.append(clean_sequence[position - start])
                    else:
                        grid_cells.append("##")
                rows.append(
                    PenaltyAsoRow(
                        row_number=len(rows) + 1,
                        parent_aso_id=parent_id,
                        penalty_aso_id=f"{parent_id}_P{local_index + 1}{penalty_base}",
                        parent_start=start,
                        penalty_aso_position=local_index + 1,
                        target_position_3to5=target_index,
                        target_position_5to3=len(clean_reversed) - target_index - 1,
                        local_rna_context=_local_rna_context(clean_reversed, target_index),
                        target_base=target_base,
                        canonical_aso_base=canonical_base,
                        penalty_aso_base=penalty_base,
                        mismatch_pair=f"{penalty_base}:{target_base}",
                        priority=priority,
                        score=score,
                        reason=reason,
                        clean_sequence=clean_sequence,
                        idt_code=idt_aso_per_position(
                            clean_sequence,
                            chemistry.base_modifications,
                            chemistry.linkages,
                        ),
                        chemistry=chemistry_label,
                        grid_cells=tuple(grid_cells),
                        penalty_grid_index=target_index - crop_start,
                    )
                )

    if not rows:
        raise AsoInputError("No penalty ASOs could be generated with the selected settings.")

    return PenaltyAsoResult(
        inputs=inputs,
        chemistry=chemistry,
        clean_reversed_rna=clean_reversed,
        aso_length=aso_length,
        crop_start=crop_start,
        crop_end=crop_end,
        header_positions=header_positions,
        header_bases=header_bases,
        parent_starts=parent_starts,
        rows=tuple(rows),
    )


def _display_sequence(
    seq: str,
    mutation_type: str,
    mutation_start_reversed: int,
    mutation_length: int,
    crop_start: int,
    row_start: int,
    aso_length: int,
    mutation_start_reversed_options: tuple[int, ...] = (),
) -> tuple[str, tuple[DisplaySpan, ...]]:
    core_start = crop_start + row_start
    core_end = core_start + aso_length - 1

    if mutation_type == "DELETION":
        gap_pos = -1
        if core_start <= mutation_start_reversed <= core_end + 1:
            gap_pos = mutation_start_reversed - core_start
        if 0 <= gap_pos <= aso_length:
            gap_width = max(mutation_length, 4)
            text = seq[:gap_pos] + (" " * gap_width) + seq[gap_pos:]
            return text, (DisplaySpan(gap_pos, gap_pos + gap_width, "gap"),)
        return seq, ()

    spans: list[DisplaySpan] = []
    for mut_start in mutation_start_reversed_options or (mutation_start_reversed,):
        mut_end = mut_start + mutation_length - 1
        overlap_start = max(core_start, mut_start)
        overlap_end = min(core_end, mut_end)
        if overlap_start <= overlap_end:
            start = overlap_start - core_start
            end = start + (overlap_end - overlap_start + 1)
            spans.append(DisplaySpan(start, end, "mutation"))
    return seq, _merge_display_spans(spans)


def _merge_display_spans(spans: list[DisplaySpan]) -> tuple[DisplaySpan, ...]:
    if not spans:
        return ()

    merged: list[DisplaySpan] = []
    for span in sorted(spans, key=lambda item: (item.start, item.end)):
        if not merged or span.start > merged[-1].end:
            merged.append(span)
            continue
        previous = merged[-1]
        merged[-1] = DisplaySpan(previous.start, max(previous.end, span.end), previous.kind)
    return tuple(merged)


def generate_design(inputs: AsoInputs, grid_warning_width: int = 123) -> AsoResult:
    chemistry = resolve_chemistry(inputs)
    mutation_type = normalise_mutation_type(inputs.mutation_type)
    if not mutation_type:
        raise AsoInputError("Mutation type must be Insertion, Deletion, or Substitution.")

    try:
        mutation_length = int(inputs.mutation_length)
        mutation_start = int(inputs.mutation_start)
    except Exception as exc:
        raise AsoInputError("Mutation length and start position (first base = 1) must be whole numbers.") from exc

    if mutation_length < 0:
        raise AsoInputError("Mutation length must be non-negative.")
    if mutation_start < 1:
        raise AsoInputError("Mutation start position (first base = 1) must be at least 1.")
    mutation_start_index = mutation_start - 1

    clean_forward = clean_rna_for_reverse(inputs.rna_sequence)
    clean_reversed = clean_forward[::-1]
    if not clean_reversed:
        raise AsoInputError("RNA sequence is empty after removing whitespace, hyphens, and underscores.")
    validate_base_sequence(clean_reversed)

    insertion_start_options_forward: tuple[int, ...] = (mutation_start_index,)
    ambiguity_warning = ""
    if mutation_type == "DELETION":
        if mutation_start_index > len(clean_reversed):
            raise AsoInputError("Deletion start is outside the cleaned RNA sequence.")
        variant_bases = 0
        mutation_start_reversed = len(clean_reversed) - mutation_start_index
        mutation_start_reversed_options = (mutation_start_reversed,)
    else:
        if mutation_length == 0:
            raise AsoInputError("Insertion/substitution length must be greater than zero.")
        if mutation_start_index + mutation_length > len(clean_reversed):
            raise AsoInputError("Insertion/substitution span is outside the cleaned RNA sequence.")
        variant_bases = mutation_length
        mutation_start_reversed = len(clean_reversed) - (mutation_start_index + mutation_length)
        if mutation_type == "INSERTION":
            insertion_start_options_forward = ambiguous_insertion_start_indexes(
                clean_forward,
                mutation_start_index,
                mutation_length,
            )
            if len(insertion_start_options_forward) > 1:
                ambiguity_warning = (
                    "Ambiguous insertion placement: the inserted bases sit within a repeat, so the exact "
                    "inserted copy cannot be uniquely numbered. The output includes the union of all possible "
                    f"insertion start positions (first base = 1): "
                    f"{_format_position_options_one_based(insertion_start_options_forward)}."
                )
        mutation_start_reversed_options = tuple(
            sorted(
                len(clean_reversed) - (start_index + mutation_length)
                for start_index in insertion_start_options_forward
            )
        )

    aso_length = chemistry.aso_length
    if aso_length <= 0:
        raise AsoInputError("ASO length must be greater than zero.")
    try:
        microwalk_step = int(inputs.microwalk_step_size)
    except Exception as exc:
        raise AsoInputError("Step size must be a whole number.") from exc
    if microwalk_step < 1:
        raise AsoInputError("Step size must be at least 1.")

    ideal_windows = tuple(
        (
            start_reversed - aso_length,
            start_reversed + variant_bases + aso_length - 1,
        )
        for start_reversed in mutation_start_reversed_options
    )
    ideal_crop_start = min(start for start, _end in ideal_windows)
    ideal_crop_end = max(end for _start, end in ideal_windows)
    crop_start = max(0, ideal_crop_start)
    crop_end = min(len(clean_reversed) - 1, ideal_crop_end)
    displayed_bases = max(0, crop_end - crop_start + 1)

    absolute_row_starts: set[int] = set()
    complete_absolute_row_starts: set[int] = set()
    for ideal_start, ideal_end in ideal_windows:
        complete_displayed = max(0, ideal_end - ideal_start + 1)
        complete_possible = max(0, complete_displayed - aso_length + 1)
        complete_absolute_row_starts.update(
            ideal_start + offset for offset in range(0, complete_possible, microwalk_step)
        )

        clipped_start = max(0, ideal_start)
        clipped_end = min(len(clean_reversed) - 1, ideal_end)
        clipped_displayed = max(0, clipped_end - clipped_start + 1)
        clipped_possible = max(0, clipped_displayed - aso_length + 1)
        absolute_row_starts.update(clipped_start + offset for offset in range(0, clipped_possible, microwalk_step))

    row_starts = tuple(start - crop_start for start in sorted(absolute_row_starts))
    required_asos = len(row_starts)
    complete_required_asos = len(complete_absolute_row_starts)
    coverage_warning = ""
    if (crop_start != ideal_crop_start or crop_end != ideal_crop_end) and required_asos < complete_required_asos:
        coverage_warning = (
            "Partial walk: the supplied RNA sequence does not include enough flanking context to generate "
            f"the complete set. Generated {required_asos} of {complete_required_asos} possible "
            f"{aso_length}-mer ASOs. Add more RNA sequence on either side of the variant to generate the full walk."
        )
    grid_width_status = "Grid width OK" if displayed_bases <= grid_warning_width else "EXTEND GRID FURTHER RIGHT"

    header_positions = tuple(range(1, displayed_bases + 1))
    header_bases = tuple(clean_reversed[crop_start + pos - 1].upper() for pos in header_positions)
    rows: list[AsoRow] = []

    for row_number, row_start in enumerate(row_starts, start=1):
        grid_cells: list[str] = []
        for pos, base in zip(header_positions, header_bases):
            if pos <= row_start or pos > row_start + aso_length:
                grid_cells.append("##")
            else:
                grid_cells.append(complement_base(base))

        clean_sequence = "".join(cell for cell in grid_cells if cell != "##")
        display_sequence, display_spans = _display_sequence(
            clean_sequence,
            mutation_type,
            mutation_start_reversed,
            mutation_length,
            crop_start,
            row_start,
            aso_length,
            mutation_start_reversed_options,
        )
        aso_id = (
            f"{inputs.target_gene}_{inputs.snp_identifier}_{inputs.chemistry_number}_ASO_{row_number}"
        )
        rows.append(
            AsoRow(
                row_number=row_number,
                aso_id=aso_id,
                idt_code=idt_aso_per_position(
                    clean_sequence,
                    chemistry.base_modifications,
                    chemistry.linkages,
                ),
                clean_sequence=clean_sequence,
                display_sequence=display_sequence,
                display_spans=display_spans,
                starting_position=row_start,
                chemistry=chemistry_display_label(chemistry),
                grid_cells=tuple(grid_cells),
                core_start=crop_start + row_start,
                core_end=crop_start + row_start + aso_length - 1,
            )
        )

    return AsoResult(
        inputs=inputs,
        chemistry=chemistry,
        clean_reversed_rna=clean_reversed,
        aso_length=aso_length,
        mutation_start_reversed=mutation_start_reversed,
        mutation_start_reversed_options=mutation_start_reversed_options,
        variant_bases=variant_bases,
        crop_start=crop_start,
        crop_end=crop_end,
        displayed_bases=displayed_bases,
        required_asos=required_asos,
        complete_required_asos=complete_required_asos,
        coverage_warning=coverage_warning,
        ambiguity_warning=ambiguity_warning,
        grid_width_status=grid_width_status,
        header_positions=header_positions,
        header_bases=header_bases,
        rows=tuple(rows),
    )


def mutation_header_indexes(result: AsoResult) -> set[int]:
    mutation_type = normalise_mutation_type(result.inputs.mutation_type)
    if mutation_type not in {"INSERTION", "SUBSTITUTION"}:
        return set()

    indexes: set[int] = set()
    mutation_length = int(result.inputs.mutation_length)
    for start in result.mutation_start_reversed_options or (result.mutation_start_reversed,):
        end = start + mutation_length - 1
        indexes.update(
            idx
            for idx, position in enumerate(result.header_positions)
            if start <= result.crop_start + position - 1 <= end
        )
    return indexes


def header_display_positions_5to3(result: AsoResult) -> tuple[int, ...]:
    return tuple(
        len(result.clean_reversed_rna) - (result.crop_start + idx)
        for idx, _position in enumerate(result.header_positions)
    )


def variant_warning_text(result: AsoResult) -> str:
    return "\n".join(message for message in (result.ambiguity_warning, result.coverage_warning) if message)


def inputs_from_strings(**kwargs: str) -> AsoInputs:
    data = dict(kwargs)
    for key in ("gap_length", "wing_length", "mutation_length", "mutation_start"):
        if key in data and data[key] != "":
            data[key] = int(data[key])
    return AsoInputs(**data)
