"""Resolve Sentaurus dopant fields before reducing them to Vela node pairs."""

from dataclasses import dataclass
from decimal import localcontext
from hashlib import sha256
import math
from pathlib import Path
import re


_BUILTIN = {
    **{name + "Concentration": ("donor", name + "ActiveConcentration", ())
       for name in ("Antimony", "Arsenic", "Phosphorus", "NDopant")},
    **{name + "Concentration": ("acceptor", name + "ActiveConcentration", ())
       for name in ("Boron", "Indium", "PDopant")},
}
_NET = "DopingConcentration"
_TOTAL = "TotalConcentration"
_AGGREGATES = {"DonorConcentration", "AcceptorConcentration"}
_OTHER_OUTPUTS = {"DonorActiveConcentration", "AcceptorActiveConcentration",
                  "IonizedDonorConcentration", "IonizedAcceptorConcentration",
                  "TotalActiveConcentration", "NetActiveConcentration",
                  "AbsDopingConcentration", "Log10AbsDopingConcentration"}


def _material_key(name):
    name = name.lower()
    return {"si": "silicon", "polysilicon": "polysi"}.get(name, name)


def _registry(path):
    data = path.read_bytes()
    contents = re.sub(r'/\*.*?\*/|(?m:^[ \t]*[#*][^\n]*)|//[^\n]*', "", data.decode("utf-8-sig"), flags=re.S)
    entries = {}
    for field, body in re.findall(r'([\w]+)\s*\{([^{}]*)\}', contents):
        declaration = re.search(r'\bdoping\s*=\s*(\w+)(?:\s*\(([^()]*)\))?', body)
        if declaration is None:
            continue
        polarity, attributes = declaration.groups()
        if polarity not in ("donor", "acceptor"):
            raise ValueError(f"unsupported datexcodes doping declaration for {field}: {polarity}")
        active = re.search(r'\bactive\s*=\s*"?([\w]+)', attributes or "")
        material = re.search(r'\bmaterial\s*=\s*([^\n;]+)', attributes or "")
        materials = tuple(re.findall(r'[\w]+', material[1])) if material else ()
        value = (polarity, active[1] if active else None, materials)
        if field in entries and entries[field] != value:
            raise ValueError(f"conflicting datexcodes declarations for {field}")
        entries[field] = value
    if not entries:
        raise ValueError(f"{path.name} contains no supported dopant declarations")
    return entries, sha256(data).hexdigest()


@dataclass(frozen=True)
class DopingPolicy:
    mode: str
    compute_net: bool
    species: dict
    registry_sources: tuple

    @classmethod
    def from_inputs(cls, mode, compute_net, paths):
        species = dict(_BUILTIN)
        sources = []
        for path in map(Path, paths):
            entries, digest = _registry(path)
            species.update(entries)
            sources.append({"file": path.name, "sha256": digest})
        return cls(mode, compute_net, species, tuple(sources))

    def is_candidate(self, name):
        return (name in {_NET, _TOTAL, "NetActive"} or name.endswith("Concentration")
                or any(name in (chemical, active) for chemical, (_, active, _) in self.species.items()))

    def resolve(self, datasets, nodes, region, material, zero):
        if self.mode == "net-only":
            if _NET not in datasets:
                raise ValueError(f"region '{region}' has no {_NET} for --doping-mode net-only")
            if any(not math.isfinite(value) for value in datasets[_NET]):
                raise ValueError(f"region '{region}' has nonfinite net doping")
            pairs = [(max(value, zero), max(-value, zero)) for value in datasets[_NET]]
            return pairs, "net_split", [_NET]

        known = _AGGREGATES | _OTHER_OUTPUTS | {_NET, _TOTAL, "NetActive"}
        selected = {"donor": [], "acceptor": []}
        species_by_identity = {}
        material = _material_key(material)
        for chemical, (polarity, active, materials) in self.species.items():
            known.add(chemical)
            known.add(active)
            eligible = not materials or any(_material_key(item) in (material, "all") for item in materials)
            field = active if active in datasets else chemical if chemical in datasets else None
            if eligible and field is not None:
                identity = active or chemical
                previous = species_by_identity.get(identity)
                if previous is not None:
                    previous_polarity, previous_field = previous
                    if polarity != previous_polarity or datasets[field] != datasets[previous_field]:
                        raise ValueError(f"region '{region}' has conflicting aliases for species '{identity}'")
                    continue
                species_by_identity[identity] = (polarity, field)
                selected[polarity].append(field)

        unknown = sorted(set(datasets) - known)
        if unknown and not self.registry_sources:
            raise ValueError(f"region '{region}' has unclassified concentration fields {unknown}; "
                             "supply --datexcodes or an explicit --doping-mode")

        net = datasets.get(_NET)
        total = datasets.get(_TOTAL)
        species_fields = selected["donor"] + selected["acceptor"]
        if self.compute_net and not species_fields:
            raise ValueError(f"region '{region}' enables ComputeDopingConcentration but has no recognized dopants")

        values = [value for field in species_fields + [key for key in (_NET, _TOTAL) if key in datasets]
                  for value in datasets[field]]
        if any(not math.isfinite(value) for value in values):
            raise ValueError(f"region '{region}' has nonfinite doping")
        if any(value < zero for field in species_fields for value in datasets[field]):
            raise ValueError(f"region '{region}' has negative independent species concentrations")

        def matches(left, right, operand_scale):
            # Exported donor and acceptor fields have independent rounding errors.
            # Their errors add even when subtraction largely cancels the operands.
            return abs(left - right) <= max(float(operand_scale), float(abs(right)), 1.0) * 8e-15

        with localcontext() as context:
            decimals = [value for value in values if hasattr(value, "as_tuple") and value]
            if decimals:
                context.prec = max(32, max(value.adjusted() for value in decimals)
                                   - min(value.as_tuple().exponent for value in decimals)
                                   + len(str(len(species_fields))) + 4)
            if species_fields:
                pairs = [(sum((datasets[field][i] for field in selected["donor"]), zero),
                          sum((datasets[field][i] for field in selected["acceptor"]), zero))
                         for i in range(len(nodes))]
                for i, (donors, acceptors) in enumerate(pairs):
                    scale = donors + acceptors
                    if net is not None and not self.compute_net and not matches(donors - acceptors, net[i], scale):
                        raise ValueError(f"region '{region}' node {nodes[i]} has inconsistent {_NET} and species")
                    if total is not None and not matches(scale, total[i], scale):
                        raise ValueError(f"region '{region}' node {nodes[i]} has inconsistent {_TOTAL} and species")
                used = species_fields + ([_NET] if net is not None and not self.compute_net else [])
                if total is not None:
                    used.append(_TOTAL)
                return pairs, "species", used
            if net is None:
                raise ValueError(f"region '{region}' has no recognized dopants or {_NET}; "
                                 "use --doping-mode reported for explicit aggregate concentrations")
            if total is None:
                return [(max(value, zero), max(-value, zero)) for value in net], "net_split", [_NET]
            for i, value in enumerate(net):
                if total[i] < abs(value):
                    raise ValueError(f"region '{region}' node {nodes[i]} has {_TOTAL} smaller than |{_NET}|")
            return [((s + n) / 2, (s - n) / 2) for n, s in zip(net, total)], "net_and_total", [_NET, _TOTAL]
