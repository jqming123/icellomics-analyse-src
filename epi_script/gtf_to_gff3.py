#!/usr/bin/env python3
"""Convert a sorted GTF file into browser-only GFF3 for JBrowse 2."""

import argparse
import gzip
import io
import sys

GFF3_ESCAPE_TABLE = str.maketrans({
    "%": "%25", "\t": "%09", "\n": "%0A", "\r": "%0D",
    ";": "%3B", "=": "%3D", "&": "%26", ",": "%2C",
})

SUBFEATURE_TYPES = (
    "exon", "CDS", "five_prime_utr", "three_prime_utr",
    "start_codon", "stop_codon", "Selenocysteine",
)


def escape_value(value):
    return value.translate(GFF3_ESCAPE_TABLE)


def split_attribute_fields(text):
    fields = []
    current = []
    in_quotes = False
    for char in text.strip():
        if char == '"':
            in_quotes = not in_quotes
            current.append(char)
        elif char == ";" and not in_quotes:
            fields.append("".join(current))
            current = []
        else:
            current.append(char)
    if current:
        fields.append("".join(current))
    return fields


def parse_attributes(text):
    attributes = {}
    for field in split_attribute_fields(text):
        field = field.strip()
        if not field:
            continue
        key, separator, raw = field.partition(" ")
        if not separator:
            continue
        raw = raw.strip()
        if len(raw) >= 2 and raw.startswith('"') and raw.endswith('"'):
            raw = raw[1:-1]
        attributes[key] = raw
    return attributes


def format_attributes(pairs):
    return ";".join(f"{key}={value}" for key, value in pairs if value)


def convert_record(fields, counters):
    seqid, source, feature, start, end, score, strand, frame, attribute_text = fields
    attributes = parse_attributes(attribute_text)
    gene_id = attributes.get("gene_id", "")
    transcript_id = attributes.get("transcript_id", "")
    counters[feature] = counters.get(feature, 0) + 1
    index = counters[feature]

    pairs = []
    if feature == "gene":
        if not gene_id:
            return None
        pairs.append(("ID", f"gene:{gene_id}"))
        if attributes.get("gene_name"):
            pairs.append(("Name", escape_value(attributes["gene_name"])))
        if attributes.get("gene_biotype"):
            pairs.append(("biotype", escape_value(attributes["gene_biotype"])))
    elif feature == "transcript":
        if not transcript_id:
            return None
        pairs.append(("ID", f"transcript:{transcript_id}"))
        if gene_id:
            pairs.append(("Parent", f"gene:{gene_id}"))
        pairs.append(("Name", escape_value(attributes.get("transcript_name") or transcript_id)))
        if attributes.get("transcript_biotype"):
            pairs.append(("biotype", escape_value(attributes["transcript_biotype"])))
    elif feature in SUBFEATURE_TYPES:
        if not transcript_id:
            return None
        pairs.append(("ID", f"{feature}:{transcript_id}:{index}"))
        pairs.append(("Parent", f"transcript:{transcript_id}"))
        if feature == "exon" and attributes.get("exon_number"):
            pairs.append(("exon_number", escape_value(attributes["exon_number"])))
    else:
        pairs.append(("ID", f"{feature}:{index}"))

    phase = frame if feature == "CDS" and frame in ("0", "1", "2") else "."
    return "\t".join([
        seqid, source, feature, start, end, score, strand, phase,
        format_attributes(pairs),
    ])


def convert_stream(src, dst):
    counters = {}
    for raw_line in src:
        line = raw_line.rstrip("\r\n")
        if not line:
            continue
        if line.startswith("#"):
            dst.write(line + "\n")
            continue
        fields = line.split("\t")
        if len(fields) != 9:
            continue
        record = convert_record(fields, counters)
        if record:
            dst.write(record + "\n")


def open_input(path):
    with open(path, "rb") as probe:
        magic = probe.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", help="GTF path, or - for stdin (must be sorted by seqid/start)")
    parser.add_argument("output", nargs="?", help="GFF3 output path (default: stdout)")
    args = parser.parse_args(argv)

    if args.input == "-":
        src = io.TextIOWrapper(sys.stdin.buffer, encoding="utf-8", errors="replace")
    else:
        src = open_input(args.input)
    try:
        if args.output:
            with open(args.output, "w", encoding="utf-8", newline="\n") as dst:
                convert_stream(src, dst)
        else:
            convert_stream(src, sys.stdout)
    finally:
        src.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
