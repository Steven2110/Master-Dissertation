import os
import pandas as pd
import re

def parse_data(subdirectory, base_path='data'):
    # Construct the full file path
    file_path = os.path.join(base_path, subdirectory, f"{subdirectory}.dat")
    
    # Initialize lists to hold the parsed data
    t0_values, interval_values, date_values = [], [], []
    x_values, y_values, z_values, megno_avg_values = [], [], [], []
    vx_values, vy_values, vz_values, megno_non_avg_values = [], [], [], []

    # Open and read the file line by line
    with open(file_path, 'r') as file:
        lines = file.readlines()

        # Iterate over lines in blocks of 3 (based on the structure)
        for i in range(0, len(lines), 3):
            # Parse the first line: t0, interval, and date
            t0, interval, date_raw = lines[i].strip().split(maxsplit=2)
            t0_values.append(float(t0))  # Store as float
            interval_values.append(float(interval))  # Store as float
            
            # Extract the date part within parentheses
            date_match = re.search(r'\((.*?)\)', date_raw)
            date = date_match.group(1) if date_match else ""
            date_values.append(date)

            # Parse the second line: x, y, z, and MEGNO (avg)
            _, x, y, z, megno_avg = lines[i + 1].strip().split(maxsplit=4)
            x_values.append(float(x))
            y_values.append(float(y))
            z_values.append(float(z))
            megno_avg_values.append(float(megno_avg))

            # Parse the third line: vx, vy, vz, and MEGNO (non-avg)
            vx, vy, vz, megno_non_avg = lines[i + 2].strip().split(maxsplit=3)
            vx_values.append(float(vx))
            vy_values.append(float(vy))
            vz_values.append(float(vz))
            megno_non_avg_values.append(float(megno_non_avg))

    # Create a DataFrame from the parsed data
    data = pd.DataFrame({
        't0': t0_values,
        'interval': interval_values,
        'date': date_values,
        'x': x_values,
        'y': y_values,
        'z': z_values,
        'MEGNO (avg)': megno_avg_values,
        'vx': vx_values,
        'vy': vy_values,
        'vz': vz_values,
        'MEGNO (non-avg)': megno_non_avg_values
    })

    # Export to CSV
    output_file_path = os.path.join(base_path, subdirectory, f"parsed_data_{subdirectory}.csv")
    data.to_csv(output_file_path, index=False, float_format='%.20f')

    print(f"Data from {file_path} has been successfully exported to {output_file_path}.")
