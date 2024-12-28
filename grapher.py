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
    file_path = os.path.join(base_path, subdirectory, f"calculated_orbital_elements_with_All_Resonance_{subdirectory}.csv")

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
    save_dir = os.path.join(base_path, subdirectory, "resonance_graphs")
    os.makedirs(save_dir, exist_ok=True)

    for i in range(1, 6):
        x = np.array(df['time_years'])
        # Resonance Relations
        y = np.array(df[f'Phi_dot_{i}'])

        # Find gaps in data where the difference exceeds the threshold
        gaps = np.abs(np.diff(y)) > 100
        segments = np.split(np.arange(len(y)), np.where(gaps)[0] + 1)

        # Create a figure with subplots for stacked visualization
        fig, axes = plt.subplots(2, 1, figsize=(14, 12), sharex=True)
        fig.suptitle(f'Resonance Analysis for {subdirectory}', fontsize=18, y=0.95)
        
        for segment in segments:
            axes[0].plot(x[segment], y[segment], color='black', label=f"Phi dot {i}" if segment[0] == 0 else "")

        axes[0].set_xlabel('Time (years)', fontsize=14)
        axes[0].set_ylabel(f'Phi Dot {i} (rad/second)', fontsize=14)
        axes[0].set_title(f'Phi Dot {i}', fontsize=16)
        axes[0].grid(True)
        axes[0].legend()

        # Critical arguments
        y = np.array(df[f'Phi_{i}_deg'])

        # Find gaps in data where the difference exceeds the threshold
        gaps = np.abs(np.diff(y)) > 100
        segments = np.split(np.arange(len(y)), np.where(gaps)[0] + 1)

        # Plot Phi values on the first subplot
        for segment in segments:
            axes[1].plot(x[segment], y[segment], color='black', label=f'Phi {i}' if segment[0] == 0 else "")

        axes[1].set_ylabel(f'Phi {i} (degrees)', fontsize=14)
        axes[1].set_title(f'Phi {i}', fontsize=16)
        axes[1].grid(True)
        axes[1].legend()

        # Save the combined plot
        save_path = os.path.join(save_dir, f"{subdirectory}_resonance_{i}.png")
        plt.savefig(save_path)
        plt.close()

        print(f"Saved stacked plot for Phi and Phi Dot in {save_path}.")

