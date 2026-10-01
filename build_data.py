import sys
sys.path.insert(0, r"C:\Users\mdenavit\Dropbox\GitHub\libdenavit-py\src")

import csv
import json
import re
from math import nan
from pathlib import Path
from libdenavit import unit_convert

def process_specimen(row,database,unit_system):
    """Process one row from the CSV file.

    Parameters
    ----------
    row : dict
        Dictionary containing the values from one CSV row.

    Returns
    -------
    dict
        Processed data to be written to the JSON file.
    """
    
    if unit_system == 'US':
        length_units = 'in'
        force_units = 'kips'
        stress_units = 'ksi'
    else:
        raise ValueError(f'Invalid unit_system: {unit_system}')
    
    specimen = dict()
    specimen['specimen_name'] = row['Specimen'].strip()
    specimen['reference'] = f'{row['Author'].strip()} {row['Year'].strip()}'
    specimen['author'] = row['Author'].strip()
    specimen['year'] = get_plain_year(row['Year'])

    section_type, member_type = database.split("_")
    specimen['section_type'] = section_type
    specimen['member_type'] = member_type
    
    try:
        specimen['Fy'] = unit_convert(float(row['Fy']), row['Fy_units'], stress_units)
    except:
        raise ValueError(f'Invalid Fy for {specimen['reference']}, specimen {specimen['specimen_name']}: {row['Fy']}')
    
    if row['Fu'].strip() == "":
        specimen['Fu'] = nan
    else:
        try:
            specimen['Fu'] = unit_convert(float(row['Fu']), row['Fu_units'], stress_units)
        except:
            raise ValueError(f'Invalid Fu for {specimen['reference']}, specimen {specimen['specimen_name']}: {row['Fu']}')
            
    specimen['fc'] = unit_convert(float(row['fc']), row['fc_units'], stress_units)  # @todo - Implement fc_type conversions
    
    if section_type == "CCFT":
        pass
    elif section_type == "RCFT":
        pass
    elif section_type == "SRC":
        pass
    else:
        raise ValueError(f'Invalid section_type: {section_type}')

    if member_type == "C+PBC":
        pass
    elif member_type == "Beams":
        pass
    elif member_type == "Other":
        pass
    else:
        raise ValueError(f'Invalid member_type: {member_type}')

    specimen['tags'] = row['Tags'].strip()  # @todo - break into list
    specimen['notes'] = row['Notes'].strip()
    
    return specimen


def get_plain_year(value: str) -> int:
    """Extract the leading integer from a string."""
    match = re.match(r"\d+", value)
    if not match:
        raise ValueError(f"No leading number found in {value!r}")
    return int(match.group())


def build_data(database,unit_system):
    """Read a CSV file, process each row, and save the results as JSON."""

    csv_file = Path(f"{database}.csv")
    json_file = Path(f"{database}.json")

    results = []

    print(f"Processing: {csv_file}")
    with open(csv_file, "r", encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f)

        for row in reader:
            processed_specimen = process_specimen(row,database,unit_system)
            results.append(processed_specimen)

    with open(json_file, "w", encoding="utf-8") as f:
        json.dump(results, f, indent=4, ensure_ascii=False)


if __name__ == "__main__":
    unit_system = 'US'
    build_data('CCFT_C+PBC',unit_system)
    build_data('RCFT_C+PBC',unit_system)
    build_data('SRC_C+PBC',unit_system)
    build_data('CCFT_Beams',unit_system)
    build_data('RCFT_Beams',unit_system)
    build_data('CCFT_Other',unit_system)
    build_data('RCFT_Other',unit_system)
