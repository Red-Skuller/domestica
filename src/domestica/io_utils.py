from pathlib import Path
import pandas as pd
from Bio import SeqIO
import re

_AA_RE = re.compile(r"[^ACDEFGHIKLMNPQRSTVWY]")


def clean_seq(seq: str) -> str:
    """Return an uppercase protein sequence containing only the 20 standard AAs."""
    if seq is None: return ""
    s = str(seq).upper()
    s = re.sub(r"\s+", "", s)
    s = s.replace("-", "").replace("*", "")
    return _AA_RE.sub("", s)


def read_input(path: Path, name_col="Name", seq_col="Sequence") -> list[dict]:
    """Reads FASTA or XLSX and returns a list of dicts: [{'id': ..., 'sequence': ...}]"""
    ext = path.suffix.lower()
    records = []

    if ext in [".fasta", ".fa"]:
        for rec in SeqIO.parse(path, "fasta"):
            records.append({"id": rec.id, "sequence": clean_seq(str(rec.seq))})

    elif ext in [".xlsx"]:
        df = pd.read_excel(path)
        # Assuming the standard columns are Name and Sequence or fallback to indexes
        name_col = name_col if name_col in df.columns else df.columns[0]
        seq_col = seq_col if seq_col in df.columns else df.columns[1]

        for _, row in df.iterrows():
            records.append({
                "id": str(row[name_col]),
                "sequence": clean_seq(str(row[seq_col]))
            })
    else:
        raise ValueError(f"Unsupported file extension: {ext}")

    return records