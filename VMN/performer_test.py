import os
import random
import csv
from decimal import Decimal, getcontext
from VMN.utils import pridict_ms, common_ion_find

getcontext().prec = 10

def read_molecules(base_dir):
    molecule_file = os.path.join(base_dir, "test_molecule.txt")
    molecules = []
    with open(molecule_file, "r") as f:
        lines = f.readlines()
        for line in lines:
            parts = line.split()
            if len(parts) > 1:
                molecules.append(parts[1].strip())
    return molecules


def simulate_molecules(base_dir, smiles_list, energy_level, num_molecules, energy_level_test, num_test):
    simulator = pridict_ms.CFMIDSimulator(base_dir)
    simulator.simulate_fragments(smiles_list, energy_level, num_molecules, energy_level_test, num_test)


def calculate_scores(base_dir, num_molecules, energy_level, num_test):
    common_ions_file = os.path.join(base_dir, f"{num_molecules}-{energy_level}-{num_test}-common_ions.csv")
    common_neutral_losses_file = os.path.join(base_dir,
                                              f"{num_molecules}-{energy_level}-{num_test}-common_neutral_losses.csv")
    target_ion = Decimal("84.08")
    target_neutral_loss = Decimal("134.04")
    tolerance = Decimal("0.02")

    def read_common_file(file_path, target):
        common_data = []
        with open(file_path, "r") as f:
            reader = csv.DictReader(f)
            for row in reader:
                value = Decimal(row['ion'] if 'ion' in row else row['neutral_loss'])
                average_intensity = Decimal(row['average_intensity'])
                common_data.append((value, average_intensity))
        common_data.sort(key=lambda x: x[1], reverse=True)
        for idx, (value, avg_intensity) in enumerate(common_data):
            if abs(value - target) <= tolerance:
                return idx + 1
        return None

    ion_rank = read_common_file(common_ions_file, target_ion)
    neutral_loss_rank = read_common_file(common_neutral_losses_file, target_neutral_loss)

    def calculate_rank_score(rank):
        if rank is None:
            return Decimal(0)
        if rank <= 5:
            return Decimal(3)
        if rank <= 10:
            return Decimal(2)
        if rank <= 15:
            return Decimal(1)
        return Decimal(0)

    ion_score = calculate_rank_score(ion_rank)
    neutral_loss_score = calculate_rank_score(neutral_loss_rank)

    return ion_score, neutral_loss_score


def main():
    base_dir = "/path/to/config-config"
    molecules = read_molecules(base_dir)
    num_tests = 3

    results = []

    for num_molecules in range(2, 9):
        for energy_level in range(3):
            total_ion_score = 0
            total_neutral_loss_score = 0
            previous_combinations = []

            for num_test in range(num_tests):
                while True:
                    selected_molecules = random.sample(molecules, num_molecules)
                    if selected_molecules not in previous_combinations:
                        previous_combinations.append(selected_molecules)
                        break

                simulate_molecules(base_dir, selected_molecules, energy_level, num_molecules, energy_level,
                                   num_test + 1)

                energy_level_str = f"energy{energy_level}"
                min_neutral_loss = 50.0

                analyzer = common_ion_find.CommonIonsAnalyzer(base_dir, energy_level_str, min_neutral_loss)
                analyzer.analyze_neutral_losses_for_each_molecule(num_molecules, energy_level, num_test + 1)
                analyzer.find_common_neutral_losses(num_molecules, energy_level, num_test + 1)
                analyzer.find_common_ions(num_molecules, energy_level, num_test + 1)

                ion_score, neutral_loss_score = calculate_scores(base_dir, num_molecules, energy_level, num_test + 1)
                total_ion_score += ion_score
                total_neutral_loss_score += neutral_loss_score

                print(f"Test {num_test + 1}: Ion Score = {ion_score}, Neutral Loss Score = {neutral_loss_score}")

            average_ion_score = total_ion_score / num_tests
            average_neutral_loss_score = total_neutral_loss_score / num_tests

            print(f"Num Molecules: {num_molecules}, Energy Level: {energy_level}")
            print(f"Average Ion Score: {average_ion_score}")
            print(f"Average Neutral Loss Score: {average_neutral_loss_score}")

            results.append((num_molecules, energy_level, average_ion_score, average_neutral_loss_score))

    results_file = os.path.join(base_dir, "test_results.csv")
    with open(results_file, mode='w', newline='', encoding='utf-8') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['num_molecules', 'energy_level', 'average_ion_score', 'average_neutral_loss_score'])
        for result in results:
            writer.writerow(result)


if __name__ == "__main__":
    main()
