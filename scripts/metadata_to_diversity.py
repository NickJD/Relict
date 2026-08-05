import pandas as pd
import argparse
import re
import sys

def extract_bracket_content(value):
    """
    Extracts text inside brackets [ ] or ( ).
    If no brackets are found, returns the original value.
    """
    if pd.isna(value) or not isinstance(value, str):
        return value

    # Regex to find content inside () or []
    match = re.search(r'[\(\[](.*?)[\)\]]', value)
    if match:
        return match.group(1).strip()
    return value

def process_excel(input_file, output_file, literal_string, sheet_name=None, skiprows=4):
    try:
        # Define the specific column indices:
        # A=0, B=1, DJ=113, EG=136, DW=126, EK=140
        target_indices = [0, 1, 113, 136, 126, 140]

        # Quick check to see if the file is wide enough
        df_headers = pd.read_excel(input_file, nrows=0, sheet_name=sheet_name, engine='openpyxl')
        max_cols = len(df_headers.columns)

        for idx in target_indices:
            if idx >= max_cols:
                print(f"Error: The file does not have enough columns. "
                      f"Index {idx} is out of bounds (File only has {max_cols} columns).")
                sys.exit(1)

        print(f"Loading file and extracting columns from sheet: '{sheet_name if sheet_name else 'First Sheet'} '...")

        # Load the 6 physical columns, skipping the metadata header rows
        # Adjust 'skiprows' if your actual data starts on a different line
        df = pd.read_excel(input_file, usecols=target_indices, sheet_name=sheet_name,
                            engine='openpyxl', skiprows=skiprows)

        # Rename the loaded columns
        df.columns = ['Col_A', 'Col_B_Brackets', 'Col_DJ', 'Col_EG', 'Col_DW', 'Col_EK']

        # CLEANUP: Remove rows where the primary ID (Col_A) is empty
        # This removes the trailing "UoG_01" junk rows from your previous output
        df = df.dropna(subset=['Col_A'])

        # INSERT the literal string into the 3rd position (index 2)
        df.insert(2, 'User_Literal_Col', literal_string)

        # Process Column B to extract content from brackets
        df['Col_B_Brackets'] = df['Col_B_Brackets'].apply(extract_bracket_content)

        # Save to TSV
        final_column_order = ['Col_A', 'Col_B_Brackets', 'User_Literal_Col', 'Col_DJ', 'Col_EG', 'Col_DW', 'Col_EK']
        df = df[final_column_order]

        df.to_csv(output_file, sep='\t', index=False)
        print(f"Success! File saved to: {output_file}")
        print(f"Inserted literal string: '{literal_string}' into Column 3.")
        print(f"Note: Skipped the first {skiprows} rows of the Excel sheet.")

    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Convert specific columns from an Excel file to a TSV with a literal string insertion.")

    parser.add_argument("--input", required=True, help="Path to the source Excel file (.xlsx, .xlsm)")
    parser.add_argument("--output", required=True, help="Path to the output TSV file")
    parser.add_argument("--col3", required=True, help="The exact string to place in the 3rd column of the TSV")
    parser.add_argument("--sheet", default=None, help="Optional: Name of the sheet to process")
    # Added skiprows as an optional argument so you can change it easily
    parser.add_argument("--skip", type=int, default=4, help="Number of header rows to skip at the top of the Excel file")

    args = parser.parse_args()

    process_excel(args.input, args.output, args.col3, args.sheet, skiprows=args.skip)
