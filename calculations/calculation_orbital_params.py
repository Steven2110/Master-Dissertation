import os
import pandas as pd
import numpy as np
from astropy.time import Time
from .functions import orbital_params

def process_orbital_elements(subdirectory, base_path='data'):
    # Construct input and output file paths
    input_file_path = os.path.join(base_path, subdirectory, f"parsed_data_{subdirectory}.csv")
    output_file_path = os.path.join(base_path, subdirectory, f"calculated_orbital_elements_{subdirectory}.csv")

    # Load the data from the existing CSV file
    data = pd.read_csv(input_file_path)

    # Initialize lists to store calculated values
    t0_values = []
    interval_values = []
    date_starts = []
    date_ends = []

    # Orbital parameters
    semi_major_axis_values = []
    eccentricity_values = []
    inclination_values = []
    Omega_rads_values = []
    Omega_degs_values = []
    omega_rads_values = []
    omega_degs_values = []
    mean_anomaly_values = []
    mean_anomaly_degs_values = []

    # Loop over each row in the DataFrame to calculate the orbital elements
    for index, row in data.iterrows():
        try:
            # Calculate the year as the difference between t1 and t0
            t0 = row['t0']
            interval = row['interval']
            date_obj_t0 = Time(t0, format='jd', scale='tt')
            date_obj_end = date_obj_t0 + (interval / 86400.0)

            t0_values.append(t0)
            interval_values.append(interval)
            date_starts.append(date_obj_t0.iso)
            date_ends.append(date_obj_end.iso)

            # Position and velocity vectors
            initial_coordinates = np.array([row['x'], row['y'], row['z']])
            initial_velocities = np.array([row['vx'], row['vy'], row['vz']])

            (semi_major_axis,
                eccentricity,
                inclination,
                Omega_rads,
                Omega_degs,
                omega_rads,
                omega_degs,
                mean_anomaly_rads,
                mean_anomaly_degs) = orbital_params(initial_coordinates, initial_velocities)

            # Orbital parameters append to list
            semi_major_axis_values.append(semi_major_axis)
            eccentricity_values.append(eccentricity)
            inclination_values.append(inclination)
            Omega_rads_values.append(Omega_rads)
            Omega_degs_values.append(Omega_degs)
            omega_rads_values.append(omega_rads)
            omega_degs_values.append(omega_degs)
            mean_anomaly_values.append(mean_anomaly_rads)
            mean_anomaly_degs_values.append(mean_anomaly_degs)
        except Exception as e:
            print(f"Error processing row {index} in {input_file_path}: {e}")
            continue

    # Create a new DataFrame with the required columns
    result_df = pd.DataFrame({
        't0': t0_values,
        'interval': interval_values,
        'date_start': date_starts,
        'date_end': date_ends,
        'semi-major axis': semi_major_axis_values,
        'eccentricity': eccentricity_values,
        'inclination': inclination_values,
        'longitude of ascending node (Ω) (rads)': Omega_rads_values,
        'longitude of ascending node (Ω) (degs)': Omega_degs_values,
        'argument of pericenter (ω) (rads)': omega_rads_values,
        'argument of pericenter (ω) (degs)': omega_degs_values,
        'mean anomaly (rads)': mean_anomaly_values,
        'mean anomaly (degs)': mean_anomaly_degs_values
    })

    # Export the new DataFrame to a CSV file
    result_df.to_csv(output_file_path, index=False, float_format='%.20f')

    print(f"Data with calculated orbital elements has been saved to {output_file_path}.")
