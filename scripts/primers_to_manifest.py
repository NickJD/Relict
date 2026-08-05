import pandas as pd
import os
import shutil
import argparse
import csv
from pathlib import Path


def process_files(input_excel, input_dir, output_dir, manifest_name):
    # Convert strings to Path objects
    input_path = Path(input_excel)
    source_dir = Path(input_dir)
    dest_dir = Path(output_dir)
    manifest_path = dest_dir / manifest_name

    # Create output directory if it doesn't exist
    dest_dir.mkdir(parents=True, exist_ok=True)

    if not source_dir.exists():
        print(f"Error: Source directory '{source_dir}' does not exist.")
        return

    print(f"Reading Excel: {input_path}")

    try:
        # Read the specific sheet
        df = pd.read_excel(input_path, sheet_name='Primer file mapping')

        # Clean up column names (remove newlines and extra whitespace)
        df.columns = [str(c).replace('\n', ' ').replace('  ', ' ').strip() for c in df.columns]

        id_col = 'Partner culture collection isolate ID'
        if id_col not in df.columns:
            print(f"Error: Could not find column '{id_col}' in the sheet.")
            print(f"Found columns: {df.columns.tolist()}")
            return

        manifest_data = []

        print("Processing rows...")
        for _, row in df.iterrows():
            sample_id = row[id_col]

            # --- SKIP INSTRUCTION ROWS ---
            # Skip rows where the ID column contains the "Ensure that the ID..." instruction
            if pd.isna(sample_id) or "Ensure that the ID" in str(sample_id):
                continue

            # Iterate through all columns to find primer file entries
            for col in df.columns:
                if col == id_col:
                    continue

                filename = row[col]

                # --- SKIP EMPTY OR INSTRUCTION CELLS ---
                # Skip if cell is empty, NaN, or contains the "Please, leave blank" instruction
                if pd.isna(filename) or "Please, leave blank" in str(filename) or str(filename).strip() == "":
                    continue

                filename = str(filename).strip()

                # Extract a "Clean" primer name (e.g., "27F" instead of the whole header)
                # Logic: Find text inside parentheses, if none, use the whole header.
                if '(' in col and ')' in col:
                    primer_name = col.split('(')[-1].replace(')', '').strip()
                else:
                    primer_name = col

                # 1. Extract the clean identifier from the header
                # Example: "Filename of primer sequence (27F)" -> "27F"
                if '(' in col and ')' in col:
                    # Split by the last opening parenthesis and remove the trailing closing one
                    primer_identifier = col.split('(')[-1].replace(')', '').strip()
                else:
                    primer_identifier = col.strip()

                # 2. Determine direction based ONLY on the extracted identifier
                if 'F' in primer_identifier.upper() and 'R' not in primer_identifier.upper():
                    direction = "forward"
                elif 'R' in primer_identifier.upper():
                    direction = "reverse"
                else:
                    # Fallback for non-standard primers like ITS1, DS1, 18S, etc.
                    direction = primer_identifier


                src_file = source_dir / filename

                if src_file.exists():
                    dest_file = dest_dir / filename
                    shutil.copy2(src_file, dest_file)

                    # Append to manifest list: [ID, filename, direction, clean_primer_name]
                    manifest_data.append([sample_id, filename, direction, primer_name])
                else:
                    print(f"Warning: File not found: {filename} (Sample: {sample_id})")

        # Write the manifest file as a true CSV (comma delimited)
        with open(manifest_path, mode='w', newline='', encoding='utf-8') as m_file:
            writer = csv.writer(m_file)
            writer.writerow(['sequence_id', 'read_file', 'direction', 'primer'])
            writer.writerows(manifest_data)

        print(f"Success! Files copied to: {dest_dir}")
        print(f"Manifest created at: {manifest_path}")

    except Exception as e:
        print(f"An error occurred: {e}")


def main():
    parser = argparse.ArgumentParser(description="Process .ab1 files based on 'Primer file mapping' Excel sheet.")

    parser.add_argument("-i", "--input", required=True, help="Path to the Excel (.xlsx) file")
    parser.add_argument("-s", "--source", required=True, help="Directory containing the source .ab1 files")
    parser.add_argument("-o", "--output", required=True, help="Directory to copy processed files to")
    parser.add_argument("-m", "--manifest", default="manifest.csv",
                        help="Name of the output manifest (default: manifest.csv)")

    args = parser.parse_args()

    process_files(args.input, args.source, args.output, args.manifest)


if __name__ == "__main__":
    main()
