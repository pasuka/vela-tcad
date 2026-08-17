#!/usr/bin/env python3
"""Regression coverage for the Sentaurus ``.par`` parameter IR.

The shipped ``.par`` corpus is not sufficient on its own to pin this contract:
it contains no ``Region`` block, no repeated section, and only a single level
of ``#include``.  These tests therefore pair the real corpus with hand-built
minimal fixtures that exercise the syntax the corpus cannot reach, so that
losslessness and fail-closed behaviour are both verified rather than assumed.
"""

from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts"))

from sentaurus_parameter_ir import (  # noqa: E402
    PARAMETER_IR_SCHEMA,
    ParParseError,
    active_blocks,
    parameter_count,
    parse_parameter_ir,
    section_index,
)

PARAMETER_IR_SCHEMA_PATH = (
    REPO / "schemas" / "vela.sentaurus_parameter_ir.v1.schema.json")

PN2D_MODELS = (
    REPO / "reference_tcad" / "pn2d_sentaurus2018" / "source" / "models.par")
SINGLEDEVICE_SILICON = (
    REPO / "reference_tcad" / "singledevice_sentaurus2018" / "source"
    / "Silicon.par")
SINGLEDEVICE_SDEVICE = (
    REPO / "reference_tcad" / "singledevice_sentaurus2018" / "source"
    / "sdevice.par")
BVMETHODS_MODELS = (
    REPO / "reference_tcad" / "bvmethods_sentaurus2018" / "source"
    / "models.par")
BVMETHODS_SRH = (
    REPO / "reference_tcad" / "bvmethods_sentaurus2018" / "source"
    / "full_physics_constant_srh.par")


def validate_against_schema(testcase: unittest.TestCase,
                            document: dict,
                            schema_path: Path) -> None:
    """Validate ``document`` against ``schema_path``.

    ``jsonschema`` is not a declared dependency of this repository, so the
    check degrades to a structural assertion when the module is unavailable,
    matching ``test_sentaurus_sde_frontend.py``.
    """
    schema = json.loads(schema_path.read_text(encoding="utf-8"))
    try:
        import jsonschema  # type: ignore[import-not-found]
    except ImportError:
        for key in schema.get("required", []):
            testcase.assertIn(key, document, f"missing required key {key!r}")
        if schema.get("additionalProperties") is False:
            unexpected = set(document) - set(schema.get("properties", {}))
            testcase.assertEqual(
                set(), unexpected,
                f"unexpected top-level keys {sorted(unexpected)}")
        return
    jsonschema.validate(document, schema)


def write_par(root: Path, name: str, body: str) -> Path:
    path = root / name
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(body.strip() + "\n", encoding="utf-8")
    return path


class ParameterIrCorpusTest(unittest.TestCase):
    """The shipped corpus must round-trip with nothing dropped."""

    def test_every_shipped_par_file_parses(self) -> None:
        for path in (PN2D_MODELS, SINGLEDEVICE_SILICON, SINGLEDEVICE_SDEVICE,
                     BVMETHODS_MODELS, BVMETHODS_SRH):
            with self.subTest(par=path.name):
                ir = parse_parameter_ir(path, path.parent)
                self.assertEqual(PARAMETER_IR_SCHEMA, ir["schema"])
                self.assertEqual([], ir["parse_errors"])
                validate_against_schema(self, ir, PARAMETER_IR_SCHEMA_PATH)

    def test_models_par_accounts_for_every_code_statement(self) -> None:
        """Losslessness: statements in == blocks + parameters + scope headers.

        This is the property that makes the IR trustworthy.  If a future
        grammar change silently swallowed a line, this count would drift.
        """
        ir = parse_parameter_ir(PN2D_MODELS, PN2D_MODELS.parent)
        self.assertEqual(105, len(ir["blocks"]))
        self.assertEqual(1301, parameter_count(ir))
        # One `Material = "Silicon" {` scope header produces neither a block
        # nor a parameter, so the corpus total is 1 + 105 + 1301.
        self.assertEqual(1407, 1 + len(ir["blocks"]) + parameter_count(ir))

    def test_bandgap_keeps_variant_unit_and_zero_kelvin_reference(self) -> None:
        """The Bandgap block carries the three facts an emitter must not guess.

        ``dEg0`` is variant-selected, the values are in eV, and ``Tpar`` is
        0 K -- so ``Eg0`` is a 0 K reference.  Treating it as a 300 K value
        shifts the gap by roughly 0.045 eV and ``ni`` by nearly a decade.
        """
        ir = parse_parameter_ir(PN2D_MODELS, PN2D_MODELS.parent)
        bandgap = section_index(ir)["Bandgap"][0]
        self.assertEqual({"kind": "material", "name": "Silicon"},
                         bandgap["scope"])
        by_name = {item["raw_name"]: item for item in bandgap["parameters"]}

        old_slotboom = by_name["dEg0(OldSlotboom)"]
        self.assertEqual("OldSlotboom", old_slotboom["variant"])
        self.assertEqual("dEg0", old_slotboom["base_name"])
        self.assertEqual([-1.595e-2], old_slotboom["values"])
        self.assertEqual("eV", old_slotboom["raw_unit"])
        # Sibling variants must stay separate rather than collapse onto dEg0.
        self.assertEqual([-4.795e-3], by_name["dEg0(Slotboom)"]["values"])
        self.assertEqual([0.0], by_name["dEg0(Bennett)"]["values"])

        self.assertEqual([1.16964], by_name["Eg0"]["values"])
        self.assertEqual([0.0], by_name["Tpar"]["values"])
        self.assertEqual("K", by_name["Tpar"]["raw_unit"])

    def test_constant_mobility_keeps_split_temperature_exponent(self) -> None:
        """``Exponent`` is a carrier pair, not a single number.

        Vela applies a hard-coded ``-2.2`` to both carriers, so the electron
        exponent has no exact target.  The IR must preserve both values for
        the mapping matrix to be able to record that as a lossy mapping.
        """
        ir = parse_parameter_ir(PN2D_MODELS, PN2D_MODELS.parent)
        constant = section_index(ir)["ConstantMobility"][0]
        by_name = {item["raw_name"]: item for item in constant["parameters"]}
        self.assertEqual("carrier_pair", by_name["Exponent"]["value_kind"])
        self.assertEqual([2.5, 2.2], by_name["Exponent"]["values"])
        self.assertEqual([1417.0, 470.5], by_name["mumax"]["values"])
        self.assertEqual("cm^2/(Vs)", by_name["mumax"]["raw_unit"])

    def test_units_are_recorded_and_never_converted(self) -> None:
        """The IR is a syntax layer: cm-based values stay cm-based."""
        ir = parse_parameter_ir(SINGLEDEVICE_SILICON,
                                SINGLEDEVICE_SILICON.parent)
        scharfetter = section_index(ir)["Scharfetter"][0]
        by_name = {item["raw_name"]: item for item in scharfetter["parameters"]}
        self.assertEqual([1.0e16, 1.0e16], by_name["Nref"]["values"])
        self.assertEqual("cm^(-3)", by_name["Nref"]["raw_unit"])
        self.assertEqual([3.0e-8, 3.0e-6], by_name["taumax"]["values"])
        self.assertEqual([-1.5, -1.5], by_name["Talpha"]["values"])

    def test_exotic_corpus_name_and_value_shapes_are_classified(self) -> None:
        """Shapes that only appear once or twice must still be understood."""
        ir = parse_parameter_ir(PN2D_MODELS, PN2D_MODELS.parent)
        parameters = [item for block in ir["blocks"]
                      for item in block["parameters"]]
        kinds = {item["value_kind"] for item in parameters}
        for expected in ("scalar", "carrier_pair", "number_list",
                         "numerical_table", "record", "empty_list",
                         "string_list"):
            self.assertIn(expected, kinds)

        by_raw = {item["raw_name"]: item for item in parameters}
        # Bracketed subscripts, e.g. Dop[1][1].
        self.assertEqual([1, 1], by_raw["Dop[1][1]"]["indices"])
        self.assertEqual("Dop", by_raw["Dop[1][1]"]["base_name"])
        # Parenthesised numeric index, e.g. efit( 0), whose value is a
        # whitespace-separated pair rather than a comma-separated one.
        self.assertEqual([0], by_raw["efit( 0)"]["indices"])
        self.assertEqual("number_list", by_raw["efit( 0)"]["value_kind"])
        # A reciprocal-spelled name must not be mistaken for an expression.
        self.assertEqual("1/kappa", by_raw["1/kappa"]["base_name"])

    def test_unsupported_sections_are_kept_and_flagged(self) -> None:
        """Unsupported sections stay in the IR so coverage stays honest."""
        ir = parse_parameter_ir(PN2D_MODELS, PN2D_MODELS.parent)
        flagged = {item["section"] for item in ir["unsupported_sections"]}
        self.assertIn("OkutoCrowell", flagged)
        self.assertIn("Lackner", flagged)
        self.assertIn("UniBo", flagged)
        sections = {block["section"] for block in ir["blocks"]}
        self.assertTrue(flagged.issubset(sections),
                        "flagged sections must still be present as blocks")


class ParameterIrIncludeTest(unittest.TestCase):
    """``#include`` semantics the shipped corpus only covers shallowly."""

    def test_fragment_included_into_material_inherits_that_scope(self) -> None:
        """This is the real ``sdevice.par`` shape, pinned explicitly."""
        ir = parse_parameter_ir(SINGLEDEVICE_SDEVICE,
                                SINGLEDEVICE_SDEVICE.parent)
        self.assertEqual(["sdevice.par", "Silicon.par"], ir["files"])
        for block in ir["blocks"]:
            self.assertEqual({"kind": "material", "name": "Silicon"},
                             block["scope"])
            self.assertEqual("Silicon.par", block["source_span"]["file"])
            self.assertEqual([{"file": "sdevice.par", "line": 2}],
                             block["include_stack"])

    def test_nested_includes_record_the_full_stack(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_par_",
                                         dir=REPO / "build") as tmp:
            root = Path(tmp)
            write_par(root, "inner.par", """
Bandgap
{
	Eg0	= 1.12	# [eV]
}
""")
            write_par(root, "middle.par", '#include "inner.par"')
            top = write_par(root, "top.par", """
Material = "Silicon" {
#include "middle.par"
}
""")
            ir = parse_parameter_ir(top, root)
            self.assertEqual(1, len(ir["blocks"]))
            block = ir["blocks"][0]
            self.assertEqual({"kind": "material", "name": "Silicon"},
                             block["scope"])
            self.assertEqual("inner.par", block["source_span"]["file"])
            self.assertEqual(
                [{"file": "top.par", "line": 2},
                 {"file": "middle.par", "line": 1}],
                block["include_stack"])

    def test_include_cycle_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_par_",
                                         dir=REPO / "build") as tmp:
            root = Path(tmp)
            write_par(root, "b.par", '#include "a.par"')
            top = write_par(root, "a.par", '#include "b.par"')
            with self.assertRaises(ParParseError) as ctx:
                parse_parameter_ir(top, root)
            self.assertIn("cycle", str(ctx.exception))

    def test_include_escaping_the_root_is_rejected(self) -> None:
        """A third-party ``.par`` must not pull in arbitrary host files."""
        with tempfile.TemporaryDirectory(prefix="vela_par_",
                                         dir=REPO / "build") as tmp:
            root = Path(tmp)
            outside = root / "outside.par"
            outside.write_text("Bandgap\n{\n\tEg0\t= 1.12\n}\n",
                               encoding="utf-8")
            inner = root / "nested"
            top = write_par(inner, "top.par", '#include "../outside.par"')
            with self.assertRaises(ParParseError) as ctx:
                parse_parameter_ir(top, inner)
            self.assertIn("escapes the import root", str(ctx.exception))

    def test_missing_include_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_par_",
                                         dir=REPO / "build") as tmp:
            root = Path(tmp)
            top = write_par(root, "top.par", '#include "absent.par"')
            with self.assertRaises(ParParseError) as ctx:
                parse_parameter_ir(top, root)
            self.assertIn("does not exist", str(ctx.exception))


class ParameterIrOverrideTest(unittest.TestCase):
    """Repeated definitions and Region scope, neither present in the corpus."""

    def test_repeated_section_keeps_both_and_marks_the_earlier(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_par_",
                                         dir=REPO / "build") as tmp:
            root = Path(tmp)
            top = write_par(root, "top.par", """
Material = "Silicon" {
Bandgap
{
	Eg0	= 1.12	# [eV]
	Chi0	= 4.05	# [eV]
}
Bandgap
{
	Eg0	= 1.16964	# [eV]
}
}
""")
            ir = parse_parameter_ir(top, root)
            self.assertEqual(2, len(ir["blocks"]),
                             "a repeated section must not overwrite the first")
            first, second = ir["blocks"]
            self.assertEqual(1, first["shadowed_by"])
            self.assertNotIn("shadowed_by", second)

            by_name = {item["raw_name"]: item for item in first["parameters"]}
            # Eg0 is redefined later, Chi0 is not.
            self.assertEqual(1, by_name["Eg0"]["shadowed_by"])
            self.assertNotIn("shadowed_by", by_name["Chi0"])

            self.assertEqual([second], active_blocks(ir))

    def test_region_scope_is_first_class(self) -> None:
        """Two regions of the same material must stay distinguishable.

        Dispatching parameters by material name alone cannot express this,
        which is why the IR keys scope by kind and name.
        """
        with tempfile.TemporaryDirectory(prefix="vela_par_",
                                         dir=REPO / "build") as tmp:
            root = Path(tmp)
            top = write_par(root, "top.par", """
Material = "Silicon" {
ConstantMobility:
{
	mumax	= 1417 ,	470.5	# [cm^2/(Vs)]
}
}

Region = "drift" {
ConstantMobility:
{
	mumax	= 1200 ,	400	# [cm^2/(Vs)]
}
}

Region = "buffer" {
ConstantMobility:
{
	mumax	= 900 ,	300	# [cm^2/(Vs)]
}
}
""")
            ir = parse_parameter_ir(top, root)
            scopes = [block["scope"] for block in ir["blocks"]]
            self.assertEqual(
                [{"kind": "material", "name": "Silicon"},
                 {"kind": "region", "name": "drift"},
                 {"kind": "region", "name": "buffer"}],
                scopes)
            # Same section name in three scopes: none may shadow another.
            for block in ir["blocks"]:
                self.assertNotIn("shadowed_by", block)
            self.assertEqual(3, len(active_blocks(ir)))

    def test_named_section_variant_is_distinct_from_the_bare_section(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_par_",
                                         dir=REPO / "build") as tmp:
            root = Path(tmp)
            top = write_par(root, "top.par", """
Material = "Silicon" {
QuantumPotentialParameters {
	gamma	= 3.6 ,	3.6
}
QuantumPotentialParameters "theta_zero" {
	theta	= 0.0 ,	0.0
}
}
""")
            ir = parse_parameter_ir(top, root)
            self.assertEqual(2, len(ir["blocks"]))
            self.assertNotIn("variant", ir["blocks"][0])
            self.assertEqual("theta_zero", ir["blocks"][1]["variant"])
            for block in ir["blocks"]:
                self.assertNotIn("shadowed_by", block)


class ParameterIrFailClosedTest(unittest.TestCase):
    """Unclassified valid syntax must abort instead of degrading to a warning."""

    def test_unknown_statement_shape_produces_no_ir(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_par_",
                                         dir=REPO / "build") as tmp:
            root = Path(tmp)
            top = write_par(root, "top.par", """
Material = "Silicon" {
Bandgap
{
	Eg0	= 1.12	# [eV]
	@@ this is valid punctuation but not a known statement @@
}
}
""")
            with self.assertRaises(ParParseError) as ctx:
                parse_parameter_ir(top, root)
            message = str(ctx.exception)
            self.assertIn("top.par:5", message)
            self.assertIn("not classified", message)

    def test_assignment_outside_a_section_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_par_",
                                         dir=REPO / "build") as tmp:
            root = Path(tmp)
            top = write_par(root, "top.par", """
Material = "Silicon" {
	Eg0	= 1.12	# [eV]
}
""")
            with self.assertRaises(ParParseError) as ctx:
                parse_parameter_ir(top, root)
            self.assertIn("outside any section block", str(ctx.exception))

    def test_unterminated_section_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_par_",
                                         dir=REPO / "build") as tmp:
            root = Path(tmp)
            top = write_par(root, "top.par", """
Material = "Silicon" {
Bandgap
{
	Eg0	= 1.12	# [eV]
}
""")
            with self.assertRaises(ParParseError) as ctx:
                parse_parameter_ir(top, root)
            self.assertIn("never closed", str(ctx.exception))

    def test_unterminated_parenthesis_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_par_",
                                         dir=REPO / "build") as tmp:
            root = Path(tmp)
            top = write_par(root, "top.par", """
Material = "Silicon" {
Table
{
	NumericalTable (
	  0.1	0.2;
}
}
""")
            with self.assertRaises(ParParseError) as ctx:
                parse_parameter_ir(top, root)
            self.assertIn("unterminated", str(ctx.exception))

    def test_error_message_identifies_file_and_line(self) -> None:
        """Diagnostics must point at the include target, not the includer."""
        with tempfile.TemporaryDirectory(prefix="vela_par_",
                                         dir=REPO / "build") as tmp:
            root = Path(tmp)
            write_par(root, "frag.par", """
Bandgap
{
	Eg0	= 1.12	# [eV]
	%%% unclassifiable %%%
}
""")
            top = write_par(root, "top.par", """
Material = "Silicon" {
#include "frag.par"
}
""")
            with self.assertRaises(ParParseError) as ctx:
                parse_parameter_ir(top, root)
            self.assertIn("frag.par:4", str(ctx.exception))


class ParameterIrCommentTest(unittest.TestCase):
    """Comment handling: prose is droppable, units are data."""

    def test_prose_comments_are_ignored_and_units_retained(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_par_",
                                         dir=REPO / "build") as tmp:
            root = Path(tmp)
            top = write_par(root, "top.par", """
* a file-level prose comment
Material = "Silicon" {
Bandgap
{ * an inline prose comment introducing the section
  * Eg = Eg0 + alpha * Tpar^2
	Eg0	= 1.12	# [eV]
}
}
""")
            ir = parse_parameter_ir(top, root)
            self.assertEqual(1, len(ir["blocks"]))
            parameters = ir["blocks"][0]["parameters"]
            self.assertEqual(1, len(parameters),
                             "prose lines must not become parameters")
            self.assertEqual("eV", parameters[0]["raw_unit"])
            self.assertEqual("1.12", parameters[0]["raw_lexeme"])


if __name__ == "__main__":
    unittest.main()
