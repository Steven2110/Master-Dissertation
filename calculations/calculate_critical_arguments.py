import os
import pandas as pd
import numpy as np

from .functions import calculate_critical_arguments, calculate_sidereal_time

def convert_to_degree(radians):
    degrees = np.rad2deg(radians)  # Convert radians to degrees
    normalized_degrees = (degrees % 360 + 360) % 360  # Normalize to 0-360
    return normalized_degrees

def process_critical_arguments(subdirectory, base_path = 'data'):
    # Construct the input and output file paths
    input_file_path = os.path.join(base_path, subdirectory, f"calculated_orbital_elements_{subdirectory}.csv")
    output_file_path = os.path.join(base_path, subdirectory, f"calculated_orbital_elements_with_Phi_{subdirectory}.csv")

    # Load the dataset
    df = pd.read_csv(input_file_path)

    # Use the 'interval' column to compute the final Julian date and sidereal time
    final_jd = df['t0'] + df['interval'] / 86400
    df['theta'] = calculate_sidereal_time(final_jd)

    # Calculate the critical arguments and append the results in radians
    critical_arguments = df.apply(
        lambda row: calculate_critical_arguments(
            1,
            2,
            row['mean anomaly (rads)'],
            row['longitude of ascending node (Ω) (rads)'],
            row['argument of pericenter (ω) (rads)'],
            row['theta']
        ),
        axis=1
    )

    # Add the critical arguments in radians as new columns
    df[['Phi_1', 'Phi_2', 'Phi_3', 'Phi_4', 'Phi_5']] = pd.DataFrame(critical_arguments.tolist(), index=df.index)

    # Convert critical arguments to degrees and add as new columns
    df[['Phi_1_deg', 'Phi_2_deg', 'Phi_3_deg', 'Phi_4_deg', 'Phi_5_deg']] = df[['Phi_1', 'Phi_2', 'Phi_3', 'Phi_4', 'Phi_5']].applymap(convert_to_degree)

    # Save the updated dataframe to a new CSV file
    df.to_csv(output_file_path, index=False)

    print(f"The updated file with critical arguments is saved at: {output_file_path}")

