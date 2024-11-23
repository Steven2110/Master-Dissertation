import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def plot_with_gaps_and_save(x, y, threshold, label, xlabel, ylabel, title, save_path):
        x = np.array(x)
        y = np.array(y)

        # Find gaps in data where the difference exceeds the threshold
        gaps = np.abs(np.diff(y)) > threshold
        segments = np.split(np.arange(len(y)), np.where(gaps)[0] + 1)

        plt.figure(figsize=(10, 6))
        for segment in segments:
            plt.plot(x[segment], y[segment], color='black', label=label if segment[0] == 0 else "")  # Label only first line
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.title(title)
        plt.grid(True)
        plt.legend()

        # Save the plot instead of showing it
        plt.savefig(save_path)
        plt.close()

def process_graph(subdirectory, base_path='data'):
     # Construct file path for input data
    file_path = os.path.join(base_path, subdirectory, f"calculated_orbital_elements_with_Phi_{subdirectory}.csv")
    
    # Check if the file exists
    if not os.path.exists(file_path):
        print(f"File not found: {file_path}. Skipping...")
        return

    # Load the CSV file
    df = pd.read_csv(file_path)
    
    # Convert interval from seconds to years
    seconds_in_a_year = 365.25 * 24 * 60 * 60
    df['time_years'] = df['interval'] / seconds_in_a_year

    # Define save directory
    save_dir = os.path.join(base_path, subdirectory, "critical_arguments")
    os.makedirs(save_dir, exist_ok=True)

    phi_columns = ['Phi_1_deg', 'Phi_2_deg', 'Phi_3_deg', 'Phi_4_deg', 'Phi_5_deg']

    # Plotting for all Phi columns with gaps
    for phi_col in phi_columns:
        save_path = os.path.join(save_dir, f"{phi_col}_vs_Time.png")
        plot_with_gaps_and_save(
            df['time_years'],
            df[phi_col],
            threshold=100,  # Adjust threshold as needed for gaps
            label=phi_col,
            xlabel='Time (years)',
            ylabel=f'{phi_col} (degrees)',
            title=f'{phi_col} vs Time for {subdirectory}',
            save_path=save_path
        )
        print(f"Saved plot for {phi_col} in {save_path}.")
