import docker
import os

class CFMIDSimulator:
    ENERGY_LEVELS = {
        0: "energy0",
        1: "energy1",
        2: "energy2"
    }

    def __init__(self, base_dir):
        self.base_dir = base_dir
        self.client = docker.from_env()
        self.molecule_file = os.path.join(base_dir, "molecule.txt")
        self.param_output_file = os.path.join(base_dir, "param_output.log")
        self.param_config_file = os.path.join(base_dir, "param_config.txt")
        self.output_file = os.path.join(base_dir, "output.log")

    def write_molecules_to_file(self, molecules):
        with open(self.molecule_file, "w") as f:
            for idx, smiles in enumerate(molecules):
                f.write(f"Molecule{idx+1} {smiles}\n")

    def run_simulation(self):
        try:
            container = self.client.containers.run(
                image="wishartlab/cfmid",
                command=f"cfm-predict /data/molecule.txt 0.001 /data/param_output.log /data/param_config.txt 1 /data/output.log 0 0",
                volumes={
                    self.base_dir: {"bind": "/data", "mode": "rw"}
                },
                platform="linux/amd64",
                detach=True
            )

            container.wait()

            logs = container.logs().decode("utf-8")
            print(logs)

        except docker.errors.ContainerError as e:
            print(f"Error running container: {e}")
        except Exception as e:
            print(f"An error occurred: {e}")
        finally:
            if 'container' in locals():
                container.remove()

    def read_and_filter_output_file(self, energy_level, num_molecules=None, energy_level_test=None, num_test=None):
        fragments = {}
        current_molecule_info = {}
        current_molecule = None
        write_data = False

        try:
            with open(self.output_file, "r") as f:
                output_data = f.read()

            for line in output_data.splitlines():
                if line.startswith("#"):
                    if line.startswith("#ID="):
                        current_molecule = line.split("=")[1].strip()
                        fragments[current_molecule] = []
                        current_molecule_info[current_molecule] = [line]
                        write_data = False
                    elif current_molecule is not None and not line.startswith("#In-silico") and not line.startswith("#PREDICTED"):
                        current_molecule_info[current_molecule].append(line)
                elif line.startswith("energy"):
                    current_energy_level = line.strip()
                    write_data = (current_energy_level == self.ENERGY_LEVELS[energy_level])
                    if write_data and current_molecule is not None:
                        current_molecule_info[current_molecule].append(current_energy_level)
                elif write_data:
                    if not line.strip():
                        write_data = False
                    else:
                        parts = line.split()
                        if len(parts) >= 2:
                            mass = parts[0]
                            intensity = parts[1]
                            fragments[current_molecule].append(f"{mass} {intensity}")

            # Write each molecule's data to its own file
            for molecule, fragment_list in fragments.items():
                if num_molecules is not None and energy_level_test is not None and num_test is not None:
                    log_file_name = f"{num_molecules}-{energy_level_test}-{num_test}-{molecule}.log"
                else:
                    log_file_name = f"{molecule}.log"
                log_file_path = os.path.join(self.base_dir, log_file_name)
                with open(log_file_path, "w") as log_file:
                    for info in current_molecule_info[molecule]:
                        log_file.write(info + "\n")
                    log_file.write("\n".join(fragment_list))

        except FileNotFoundError as e:
            print(f"Output file not found: {e}")
        except Exception as e:
            print(f"An error occurred while reading the output file: {e}")

        return fragments

    def simulate_fragments(self, smiles_list, energy_level, num_molecules=None, energy_level_test=None, num_test=None):
        self.write_molecules_to_file(smiles_list)
        self.run_simulation()
        self.read_and_filter_output_file(energy_level, num_molecules, energy_level_test, num_test)



def main():
    base_dir = "/path/to/config-config"
    simulator = CFMIDSimulator(base_dir)

    smiles_list = input("Please enter the SMILES representation of the molecule(s) (separated by commas):").split(',')

    energy_level = int(input("Please enter the predicted energy level (0 for 10eV, 1 for 20eV, 2 for 40eV):"))

    simulator.simulate_fragments(smiles_list, energy_level)


if __name__ == "__main__":
    main()
