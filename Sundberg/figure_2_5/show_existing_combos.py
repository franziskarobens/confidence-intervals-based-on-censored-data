import argparse
import re
from pathlib import Path

import pandas as pd

CSV_DIR = Path(__file__).resolve().parent / "non_coverage_probs"

# ---------------- set your filter here ----------------
THETA  = None      # e.g. 1, [0.5, 1], or None for all
ALPHA  = 0.1     # e.g. 0.05, [0.05, 0.1], or None
NU     = 0.5     # e.g. 0.5, [1, 2], float("inf"), or None
DIGITS = None     # e.g. 100, or None
KIND   = "two_sided"     # e.g. "two_sided", ["lower", "upper"], or None
# -------------------------------------------------------

COLUMNS = ["theta", "alpha", "nu", "Digits", "two/onesided"]

PATTERN = re.compile(
    r"^non_coverage_(?P<kind>.+?)"
    r"_nu(?P<nu>[^_]+)"
    r"_alpha(?P<alpha>[^_]+)"
    r"_theta(?P<theta>[^_]+)"
    r"_Digits(?P<Digits>\d+)\.csv$"
)


def build_table(csv_dir: Path = CSV_DIR) -> pd.DataFrame:
    rows = []
    for f in sorted(csv_dir.glob("non_coverage_*.csv")):
        m = PATTERN.match(f.name)
        if m is None:
            print(f"warning: could not parse {f.name}")
            continue
        rows.append({
            "theta": float(m["theta"]),
            "alpha": float(m["alpha"]),
            "nu": float(m["nu"]),
            "Digits": int(m["Digits"]),
            "two/onesided": m["kind"],
        })
    df = pd.DataFrame(rows, columns=COLUMNS)
    return df.sort_values(COLUMNS).reset_index(drop=True)


def filter_table(df: pd.DataFrame, **filters) -> pd.DataFrame:
    """
    Filter by column values (single value or list); None = no filter.
    'kind' is a shortcut for the column 'two/onesided'.
    """
    mask = pd.Series(True, index=df.index)
    for key, value in filters.items():
        if value is None:
            continue
        col = "two/onesided" if key == "kind" else key
        if col not in df.columns:
            raise KeyError(f"Unknown column '{key}'. Available: {COLUMNS}")
        values = value if isinstance(value, (list, tuple, set)) else [value]
        mask &= df[col].isin(values)
    return df[mask].reset_index(drop=True)


def main(argv=None):
    p = argparse.ArgumentParser(description="List existing non-coverage CSV combinations")
    p.add_argument("--dir", type=Path, default=CSV_DIR)
    p.add_argument("--theta", type=float, nargs="+")
    p.add_argument("--alpha", type=float, nargs="+")
    p.add_argument("--nu", type=float, nargs="+")
    p.add_argument("--Digits", type=int, nargs="+")
    p.add_argument("--kind", type=str, nargs="+")
    args = p.parse_args(argv)

    df = build_table(args.dir)
    # command-line arguments override the filter values set at the top
    sub = filter_table(
        df,
        theta=args.theta or THETA,
        alpha=args.alpha or ALPHA,
        nu=args.nu or NU,
        Digits=args.Digits or DIGITS,
        kind=args.kind or KIND,
    )
    print(sub.to_string(index=False))
    print(f"\n{len(sub)} of {len(df)} combinations")
    return df, sub


if __name__ == "__main__":
    df, sub = main([])  # [] = ignore IPython/Positron arguments, use the settings above