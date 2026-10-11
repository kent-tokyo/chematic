"""Stream an Excel compound table through chematic and write an XLSX report.

Excel I/O is delegated to elixcee; chematic only handles chemistry. The first
worksheet row must contain a SMILES column. Invalid and blank inputs remain in
the output with their source row and an error instead of being silently lost.

Run:
    pip install "chematic[excel]"
    python examples/xlsx_batch_analysis.py compounds.xlsx results.xlsx
    python examples/xlsx_batch_analysis.py compounds.xlsx results.xlsx --sheet Compounds
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import sys
import uuid
from typing import Any, Iterable, Iterator

import chematic

try:
    import elixcee
except ImportError as exc:  # pragma: no cover - depends on the optional extra
    raise SystemExit(
        'Excel support requires ElixCee: pip install "chematic[excel]"'
    ) from exc


SMILES_COLUMN_CANDIDATES = ("smiles", "canonical_smiles", "structure")
NAME_COLUMN_CANDIDATES = ("name", "compound", "compound_name", "id")
OUTPUT_COLUMNS = (
    "source_row",
    "name",
    "input_smiles",
    "canonical_smiles",
    "formula",
    "mw",
    "logp",
    "tpsa",
    "hbd",
    "hba",
    "qed",
    "lipinski_passes",
    "pains_passes",
    "parse_status",
    "error",
)
MAX_EXCEL_DATA_ROWS = 1_048_575  # Excel row limit minus the header row.
MAX_EXCEL_ROWS = 1_048_576
MAX_EXCEL_COLUMNS = 16_384


def normalized_headers(values: Iterable[Any]) -> list[str]:
    headers = ["" if value is None else str(value).strip() for value in values]
    if not headers or not any(headers):
        raise ValueError("the first worksheet row must contain column names")

    seen: dict[str, str] = {}
    for header in headers:
        if not header:
            continue
        key = header.casefold()
        if key in seen:
            raise ValueError(
                f"duplicate worksheet columns are ambiguous: {seen[key]!r} and {header!r}"
            )
        seen[key] = header
    return headers


def resolve_column(
    headers: list[str], requested: str | None, candidates: tuple[str, ...], label: str
) -> int | None:
    by_name = {
        header.casefold(): index for index, header in enumerate(headers) if header
    }
    if requested is not None:
        index = by_name.get(requested.strip().casefold())
        if index is None:
            raise ValueError(f"{label} column {requested!r} was not found")
        return index
    for candidate in candidates:
        if candidate in by_name:
            return by_name[candidate]
    return None


def value_at(row: list[Any], index: int | None) -> Any:
    return None if index is None or index >= len(row) else row[index]


def analyze_row(source_row: int, name_value: Any, smiles_value: Any) -> dict[str, Any]:
    name = "" if name_value is None else str(name_value).strip()
    if not name:
        name = f"row_{source_row}"

    result = {
        "source_row": source_row,
        "name": name,
        "input_smiles": "" if smiles_value is None else str(smiles_value),
        "canonical_smiles": "",
        "formula": "",
        "mw": None,
        "logp": None,
        "tpsa": None,
        "hbd": None,
        "hba": None,
        "qed": None,
        "lipinski_passes": None,
        "pains_passes": None,
        "parse_status": "error",
        "error": "",
    }
    if not isinstance(smiles_value, str) or not smiles_value.strip():
        result["error"] = "SMILES cell must be a non-empty text value"
        return result

    smiles = smiles_value.strip()
    result["input_smiles"] = smiles
    try:
        mol = chematic.from_smiles(smiles)
        if mol is None:
            raise ValueError("SMILES parser returned no molecule")
    except (TypeError, ValueError, RuntimeError) as exc:
        result["error"] = str(exc) or type(exc).__name__
        return result

    result.update(
        {
            "canonical_smiles": mol.smiles,
            "formula": mol.formula,
            "mw": round(mol.mw, 4),
            "logp": round(mol.logp, 4),
            "tpsa": round(mol.tpsa, 4),
            "hbd": mol.hbd,
            "hba": mol.hba,
            "qed": round(mol.qed, 6),
            "lipinski_passes": mol.lipinski_passes,
            "pains_passes": mol.pains_passes,
            "parse_status": "ok",
            "error": "",
        }
    )
    return result


def iter_analysis_rows(
    input_path: Path,
    *,
    sheet: str | None,
    smiles_column: str | None,
    name_column: str | None,
    max_rows: int,
) -> Iterator[dict[str, Any]]:
    with elixcee.open_stream(
        str(input_path),
        sheet=sheet,
        include_row_numbers=True,
        # Read one extra data row so --max-rows never truncates silently.
        max_rows=min(max_rows + 2, MAX_EXCEL_ROWS),
        max_columns=MAX_EXCEL_COLUMNS,
    ) as reader:
        try:
            _header_row, header_values = next(reader)
        except StopIteration as exc:
            raise ValueError("workbook has no rows") from exc

        headers = normalized_headers(header_values)
        smiles_index = resolve_column(
            headers, smiles_column, SMILES_COLUMN_CANDIDATES, "SMILES"
        )
        if smiles_index is None:
            choices = ", ".join(SMILES_COLUMN_CANDIDATES)
            raise ValueError(
                f"no SMILES column found; expected one of {choices}, or use --smiles-column"
            )
        name_index = resolve_column(
            headers, name_column, NAME_COLUMN_CANDIDATES, "name"
        )

        data_rows = 0
        for source_row, values in reader:
            data_rows += 1
            if data_rows > max_rows:
                raise ValueError(f"workbook exceeds --max-rows ({max_rows})")
            yield analyze_row(
                source_row,
                value_at(values, name_index),
                value_at(values, smiles_index),
            )


def write_report(
    output_path: Path, rows: Iterable[dict[str, Any]], max_rows: int
) -> tuple[int, int]:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    temporary = output_path.with_name(
        f".{output_path.name}.{uuid.uuid4().hex}.tmp.xlsx"
    )
    accepted = 0
    rejected = 0
    try:
        with elixcee.create_stream_bounded(
            str(temporary),
            max_rows=max_rows + 1,
            max_columns=len(OUTPUT_COLUMNS),
        ) as writer:
            writer.append(OUTPUT_COLUMNS)
            for row in rows:
                writer.append([row[column] for column in OUTPUT_COLUMNS])
                if row["parse_status"] == "ok":
                    accepted += 1
                else:
                    rejected += 1
        os.replace(temporary, output_path)
    except Exception:
        temporary.unlink(missing_ok=True)
        raise
    return accepted, rejected


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="input .xlsx or .xlsm workbook")
    parser.add_argument("output", type=Path, help="output .xlsx report")
    parser.add_argument("--sheet", help="worksheet name; defaults to the first sheet")
    parser.add_argument(
        "--smiles-column", help="SMILES header when auto-detection is insufficient"
    )
    parser.add_argument("--name-column", help="optional compound-name header")
    parser.add_argument(
        "--max-rows",
        type=int,
        default=100_000,
        help="maximum data rows to read and write (default: 100000)",
    )
    args = parser.parse_args(argv)
    if args.max_rows < 1 or args.max_rows > MAX_EXCEL_DATA_ROWS:
        parser.error(f"--max-rows must be between 1 and {MAX_EXCEL_DATA_ROWS}")
    if args.output.suffix.casefold() != ".xlsx":
        parser.error("output must use the .xlsx extension")
    return args


def main(argv: list[str] | None = None) -> int:
    args = parse_args(sys.argv[1:] if argv is None else argv)
    try:
        rows = iter_analysis_rows(
            args.input,
            sheet=args.sheet,
            smiles_column=args.smiles_column,
            name_column=args.name_column,
            max_rows=args.max_rows,
        )
        accepted, rejected = write_report(args.output, rows, args.max_rows)
    except (OSError, RuntimeError, TypeError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1

    print(f"Wrote {accepted + rejected} rows to {args.output}")
    print(f"  analysed: {accepted}")
    print(f"  rejected: {rejected}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
