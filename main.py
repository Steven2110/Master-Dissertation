from data.data_parser import parse_data
from calculations.calculation_orbital_params import process_orbital_elements, process_differential_orbital_elements
from calculations.calculate_critical_arguments import process_critical_arguments, process_resonance_relations
from grapher import process_graph

def main():
    subdirectories = ["EPH_0001", "EPH_2422", "EPH_2570"]

    for subdir in subdirectories:
        # Parsing
        print(f"Processing for {subdir}!\n")
        print(f"Parsing data for {subdir}!")
        parse_data(subdir)
        print(f"Parsing data for {subdir} finished!\n")

        # Calculate orbital parameters
        print(f"Calculating orbital parameters for {subdir}!\n")
        process_orbital_elements(subdir)
        print(f"Calculation for {subdir} finished!\n")

        # Calculate differential orbital parameters
        print(f"Calculating differential orbital parameters for {subdir}!\n")
        process_differential_orbital_elements(subdir)
        print(f"Calculation for {subdir} finished!\n")

        # Calculate critical arguments
        print(f"Calculating critical arguments for {subdir}!\n")
        process_critical_arguments(subdir)
        print(f"Calculation for {subdir} finished!\n")

        # Calculate resonance relations
        print(f"Calculating resonance relations for {subdir}!\n")
        process_resonance_relations(subdir)
        print(f"Calculation for {subdir} finished!\n")

        # Draw and save graph
        print(f"Drawing critical arguments graph for {subdir}!\n")
        process_graph(subdir)
        print(f"Drawing and saving critical arguments graph for {subdir} finished!\n")

main()