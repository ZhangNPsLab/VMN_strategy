import os
import csv
from decimal import Decimal, getcontext, ROUND_DOWN

# 设置精度，可以根据需要进行调整
getcontext().prec = 10


class CommonIonsAnalyzer:
    def __init__(self, base_dir, energy_level, min_neutral_loss=0):
        self.base_dir = base_dir
        self.energy_level = energy_level
        self.min_neutral_loss = Decimal(min_neutral_loss)

    def get_molecule_files(self, num_molecules=None, energy_level=None, num_test=None):
        if num_molecules is None or energy_level is None or num_test is None:
            return [os.path.join(self.base_dir, f) for f in os.listdir(self.base_dir) if
                    f.startswith("Molecule") and f.endswith(".log")]
        else:
            return [os.path.join(self.base_dir, f) for f in os.listdir(self.base_dir) if
                    f.startswith(f"{num_molecules}-{energy_level}-{num_test}-Molecule") and f.endswith(".log")]

    def extract_top_ions(self, file_path, top_n=30):
        ions = []
        with open(file_path, 'r') as f:
            lines = f.readlines()
            start_collecting = False
            for line in lines:
                if line.strip() == self.energy_level:
                    start_collecting = True
                    continue
                if start_collecting:
                    if not line.strip():
                        break
                    parts = line.split()
                    if len(parts) >= 2:
                        mz = float(parts[0])
                        intensity = float(parts[1])
                        ions.append((mz, intensity))
        ions.sort(key=lambda x: x[1], reverse=True)
        return ions[:top_n]

    def generate_neutral_losses(self, top_ions, charge=1):
        neutral_losses = {}
        for i in range(len(top_ions)):
            mz1, intensity1 = top_ions[i]
            for j in range(len(top_ions)):
                if i != j:
                    mz2, intensity2 = top_ions[j]
                    neutral_loss = (Decimal(mz1) - Decimal(mz2)) * charge
                    if neutral_loss >= self.min_neutral_loss:
                        average_intensity = (Decimal(intensity1) + Decimal(intensity2)) / Decimal(2.0)
                        neutral_loss = neutral_loss.quantize(Decimal('1.00000'), rounding=ROUND_DOWN)
                        mz1 = Decimal(mz1).quantize(Decimal('1.00000'), rounding=ROUND_DOWN)
                        mz2 = Decimal(mz2).quantize(Decimal('1.00000'), rounding=ROUND_DOWN)
                        average_intensity = average_intensity.quantize(Decimal('1.00000'), rounding=ROUND_DOWN)

                        found = False
                        for nl in neutral_losses.keys():
                            if abs(neutral_loss - nl) < Decimal('1e-5'):
                                neutral_losses[nl].append((mz1, mz2, average_intensity))
                                found = True
                                break
                        if not found:
                            neutral_losses[neutral_loss] = [(mz1, mz2, average_intensity)]
        return neutral_losses

    def write_neutral_losses_to_csv(self, neutral_losses, file_name):
        csv_file_path = os.path.join(self.base_dir, file_name)

        total_intensities = {}
        for nl, mz_pairs in neutral_losses.items():
            intensities = [pair[2] for pair in mz_pairs]
            total_intensities[nl] = sum(intensities)

        max_total_intensity = max(total_intensities.values())

        with open(csv_file_path, mode='w', newline='', encoding='utf-8') as csv_file:
            writer = csv.writer(csv_file)
            writer.writerow(['neutral_loss', 'Intensity', 'm1', 'm2'])

            for nl, mz_pairs in neutral_losses.items():
                total_intensity = total_intensities[nl]

                normalized_intensity = (total_intensity / max_total_intensity) * 100

                row = [nl, normalized_intensity]
                for mz1, mz2, _ in mz_pairs:
                    row.extend([mz1, mz2])

                writer.writerow(row)

    def analyze_neutral_losses_for_each_molecule(self, num_molecules=None, energy_level=None, num_test=None):
        molecule_files = self.get_molecule_files(num_molecules, energy_level, num_test)
        for file in molecule_files:
            top_ions = self.extract_top_ions(file)
            neutral_losses = self.generate_neutral_losses(top_ions)
            file_name = os.path.basename(file).replace('.log', '_neutral_losses.csv')
            self.write_neutral_losses_to_csv(neutral_losses, file_name)

    def find_common_neutral_losses(self, num_molecules=None, energy_level=None, num_test=None):
        if num_molecules is None or energy_level is None or num_test is None:
            neutral_loss_files = [f for f in os.listdir(self.base_dir)
                                  if f.startswith("Molecule") and f.endswith(".csv")]
        else:
            neutral_loss_files = [f for f in os.listdir(self.base_dir)
                                  if f.startswith(f"{num_molecules}-{energy_level}-{num_test}-Molecule") and f.endswith(
                    ".csv")]
        all_neutral_losses = {}

        for file in neutral_loss_files:
            file_path = os.path.join(self.base_dir, file)
            with open(file_path, 'r') as csvfile:
                reader = csv.DictReader(csvfile)
                for row in reader:
                    neutral_loss = Decimal(row['neutral_loss'])
                    average_intensity = Decimal(row['average_intensity'])
                    found = False
                    for nl in all_neutral_losses.keys():
                        if abs(neutral_loss - nl) < Decimal('1e-5'):
                            all_neutral_losses[nl].append(average_intensity)
                            found = True
                            break
                    if not found:
                        all_neutral_losses[neutral_loss] = [average_intensity]

        common_neutral_losses = {nl: intensities for nl, intensities in all_neutral_losses.items() if
                                 len(intensities) == len(neutral_loss_files)}

        common_neutral_losses_avg = {nl: sum(intensities) / len(intensities) for nl, intensities in
                                     common_neutral_losses.items()}

        if num_molecules is None or energy_level is None or num_test is None:
            common_csv_path = os.path.join(self.base_dir, "common_neutral_losses.csv")
        else:
            common_csv_path = os.path.join(self.base_dir,
                                           f"{num_molecules}-{energy_level}-{num_test}-common_neutral_losses.csv"
                                           )

        common_neutral_losses_sorted = sorted(common_neutral_losses_avg.items(), key=lambda x: x[1], reverse=True)

        with open(common_csv_path, mode='w', newline='', encoding='utf-8') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(['neutral_loss', 'average_intensity'])
            for nl, avg_intensity in common_neutral_losses_sorted:
                writer.writerow([float(nl), float(avg_intensity)])

    def find_common_ions(self, num_molecules=None, energy_level=None, num_test=None):
        if num_molecules is None or energy_level is None or num_test is None:
            molecule_files = [f for f in os.listdir(self.base_dir) if f.startswith("Molecule") and f.endswith(".log")]
        else:
            molecule_files = [f for f in os.listdir(self.base_dir) if
                              f.startswith(f"{num_molecules}-{energy_level}-{num_test}-Molecule") and f.endswith(
                                  ".log")]

        all_ions = {}

        for file in molecule_files:
            file_path = os.path.join(self.base_dir, file)
            with open(file_path, 'r') as logfile:
                lines = logfile.readlines()
                energy_section = False
                for line in lines:
                    line = line.strip()
                    if line.startswith("energy"):
                        energy_section = True
                        continue
                    if energy_section:
                        if line:
                            parts = line.split()
                            if len(parts) >= 2:
                                ion = Decimal(parts[0])
                                intensity = Decimal(parts[1])
                                found = False
                                for i in all_ions.keys():
                                    if abs(ion - i) < Decimal('1e-5'):
                                        all_ions[i].append(intensity)
                                        found = True
                                        break
                                if not found:
                                    all_ions[ion] = [intensity]

        common_ions = {ion: intensities for ion, intensities in all_ions.items() if
                       len(intensities) == len(molecule_files)}

        common_ions_avg = {ion: sum(intensities) / len(intensities) for ion, intensities in common_ions.items()}

        if num_molecules is None or energy_level is None or num_test is None:
            common_csv_path = os.path.join(self.base_dir, "common_ions.csv")
        else:
            common_csv_path = os.path.join(self.base_dir,
                                           f"{num_molecules}-{energy_level}-{num_test}-common_ions.csv")

        common_ions_sorted = sorted(common_ions_avg.items(), key=lambda x: x[1], reverse=True)

        with open(common_csv_path, mode='w', newline='', encoding='utf-8') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(['ion', 'average_intensity'])
            for ion, avg_intensity in common_ions_sorted:
                writer.writerow([float(ion), float(avg_intensity)])


def main():
    base_dir = "/path/to/config-config"
    energy_level = input("Please enter the predicted energy level (0 for 10eV, 1 for 20eV, 2 for 40eV):").strip()
    energy_level = f"energy{energy_level}"
    min_neutral_loss = input("Please enter the minimum value for neutral loss (default is 50):").strip()
    min_neutral_loss = float(min_neutral_loss) if min_neutral_loss else 50.0

    analyzer = CommonIonsAnalyzer(base_dir, energy_level, min_neutral_loss)
    analyzer.analyze_neutral_losses_for_each_molecule()
    analyzer.find_common_neutral_losses()
    analyzer.find_common_ions()


if __name__ == "__main__":
    main()
