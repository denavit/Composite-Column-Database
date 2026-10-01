from pathlib import Path
import csv
import io

import requests


def download_sheets(
    spreadsheet_id: str,
    sheets: dict[str, int],
    output_dir: str | Path = ".",
) -> None:
    """Download worksheets from a public Google Sheet as CSV files.

    Parameters
    ----------
    spreadsheet_id : str
        Google Sheets spreadsheet ID.
    sheets : dict[str, int]
        Mapping of output filename to worksheet gid.
        For example: {"k_series": 0, "lh_series": 123456789}.
    output_dir : str or Path, optional
        Directory in which to save the CSV files.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    for name, gid in sheets.items():
        url = (
            f"https://docs.google.com/spreadsheets/d/{spreadsheet_id}/"
            f"export?format=csv&gid={gid}"
        )

        response = requests.get(url, timeout=30)
        response.raise_for_status()

        # Parse Google's CSV and rewrite with minimal quoting.
        text = response.content.decode("utf-8")
        rows = csv.reader(io.StringIO(text))
        
        output_file = output_dir / f"{name}.csv"

        with output_file.open(
            "w", newline="", encoding="utf-8"
        ) as f:
            writer = csv.writer(f, quoting=csv.QUOTE_MINIMAL)
            writer.writerows(rows)

        print(f"Downloaded: {output_file}")


if __name__ == "__main__":
    download_sheets(
        spreadsheet_id="1B-o75HJTP9XQo7Jq5ii_OzN86Ddwo2XWRBEYxZy1Q20",
        sheets={
            "CCFT_C+PBC": 0,
            "RCFT_C+PBC": 1,
            "SRC_C+PBC":  2,
            "CCFT_Beams": 3,
            "RCFT_Beams": 6,
            "CCFT_Other": 4,
            "RCFT_Other": 5,
        },
    )
