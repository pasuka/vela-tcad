"""Derive the four JD inverse-series coefficients from the Fermi fugacity law.

Independent Decimal series reversion; no device or native-tool data is used.
Run from the repository root with --output <ignored evidence path> to save JSON.
"""
from __future__ import annotations

import argparse
import json
from decimal import Decimal, localcontext
from pathlib import Path


def derive(precision: int = 90) -> dict:
    if precision < 40:
        raise ValueError("Use at least 40 decimal digits for the cancelling coefficients")
    with localcontext() as context:
        context.prec = precision
        count = 6  # z through c^5 is required for log(z/c) through c^4.

        def multiply(a, b):
            return [sum((a[j] * b[i - j] for j in range(i + 1)), Decimal(0))
                    for i in range(count)]

        def power(a, exponent):
            result = [Decimal(1)] + [Decimal(0)] * (count - 1)
            for _ in range(exponent):
                result = multiply(result, a)
            return result

        fermi = [Decimal(0)] + [Decimal((-1) ** (k + 1)) /
                    (Decimal(k) * Decimal(k).sqrt()) for k in range(1, count)]
        inverse = [Decimal(0), Decimal(1)] + [Decimal(0)] * (count - 2)
        for order in range(2, count):
            inverse[order] = -sum((fermi[k] * power(inverse, k)[order]
                                  for k in range(2, count)), Decimal(0))
        remainder = inverse[1:] + [Decimal(0)]
        remainder[0] -= 1
        logarithm = [sum((Decimal((-1) ** (k + 1)) / Decimal(k) *
                         power(remainder, k)[order] for k in range(1, count)),
                        Decimal(0)) for order in range(count)]
        return {"precision_digits": precision,
                "fermi_series": [str(x) for x in fermi[1:]],
                "inverse_log_coefficients": [str(x) for x in logarithm[1:5]],
                "inverse_fugacity": [str(x) for x in inverse[1:]]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--precision", type=int, default=90)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    payload = json.dumps(derive(args.precision), indent=2) + "\n"
    if args.output:
        args.output.write_text(payload, encoding="utf-8")
    else:
        print(payload, end="")


if __name__ == "__main__":
    main()
