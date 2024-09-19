import re
from django.http import JsonResponse
import os
from VMN.auto_neutral_losses import generate_neutral_losses, calculate_neutral_loss_percentages


def show_feature(request):
    if request.method == 'POST':
        analyze = request.POST.get('analyze')
        energy_level = request.POST.get('energyLevel')
        tolerance = float(request.POST.get('tolerance', 0.00001))

        if not analyze or not energy_level:
            return JsonResponse({'status': 'error', 'message': 'Missing required parameters.'}, status=400)

        response = {'status': 'error', 'message': 'Invalid processing path'}

        if energy_level == 'expt':
            selected_molecules = None
            if analyze == 'select':
                selected_molecules = request.POST.get('selectedMolecules', None)
                if selected_molecules:
                    selected_molecules = list(map(int, selected_molecules.split(' ')))

            response = handle_mgf_file_analysis(request, tolerance, selected_molecules)

        elif energy_level in ['10eV', '20eV', '40eV']:
            selected_molecules = None
            if analyze == 'select':
                selected_molecules = request.POST.get('selectedMolecules', None)
                if selected_molecules:
                    selected_molecules = list(map(int, selected_molecules.split(' ')))

            response = handle_simulation_file_analysis(request, energy_level, tolerance, selected_molecules)

        return JsonResponse({'status': 'success', 'data': response}, status=200)

    return JsonResponse({'status': 'error', 'message': 'Invalid request method.'}, status=405)


def handle_simulation_file_analysis(request, energy_level, tolerance, selected_molecules=None):
    user_directory = request.session.get('user_directory')
    output_file = os.path.join(user_directory, "output.log")

    if not os.path.exists(output_file):
        return {'status': 'error', 'message': 'Output log file not found'}

    spectrum_data = {}
    collecting_data = False
    current_energy = None
    molecule_count = None

    with open(output_file, 'r') as file:
        for line in file:
            line = line.strip()

            if line.startswith("#ID="):
                molecule_count = int(line.split("=")[-1].replace("Molecule", ""))
                collecting_data = False

                if selected_molecules and molecule_count not in selected_molecules:
                    continue

            # 解析能量层级
            if line.startswith("energy0"):
                current_energy = "10eV"
            elif line.startswith("energy1"):
                current_energy = "20eV"
            elif line.startswith("energy2"):
                current_energy = "40eV"

            if current_energy == energy_level:
                collecting_data = True
            elif line.startswith("energy") and collecting_data:
                collecting_data = False

            if collecting_data and line and line[0].isdigit():
                mz, intensity = map(float, line.split())
                if molecule_count not in spectrum_data:
                    spectrum_data[molecule_count] = []
                spectrum_data[molecule_count].append([mz, intensity])

            if collecting_data and (not line or not line[0].isdigit()):
                collecting_data = False

    if not spectrum_data:
        return {'status': 'error', 'message': 'No spectrum data found for selected energy level'}

    return process_spectrum_data(spectrum_data, tolerance)


def process_spectrum_data(spectrum_data, tolerance):
    nl_data = {}
    for molecule_id, data in spectrum_data.items():
        max_intensity = max([intensity for _, intensity in data])
        spectrum_data[molecule_id] = [[mz, (intensity / max_intensity) * 100] for mz, intensity in data]

        top_ions = sorted(spectrum_data[molecule_id], key=lambda x: x[1], reverse=True)[:30]
        neutral_losses = generate_neutral_losses(top_ions, 50)
        nl_percentages = calculate_neutral_loss_percentages(neutral_losses)
        nl_data[molecule_id] = nl_percentages

    mz_intensities = {}
    nl_intensities = {}

    first_molecule = list(spectrum_data.keys())[0]
    for mz, intensity in spectrum_data[first_molecule]:
        mz_intensities[mz] = intensity

    for nl, nl_intensity in nl_data[first_molecule].items():
        nl_intensities[nl] = nl_intensity

    for molecule_id, data in spectrum_data.items():
        if molecule_id == first_molecule:
            continue

        current_mz_intensities = {}
        for mz_2, intensity_2 in data:
            for mz in list(mz_intensities.keys()):
                if abs(mz - mz_2) <= tolerance:
                    mz_intensities[mz] += intensity_2
                    current_mz_intensities[mz] = mz_intensities[mz]
        mz_intensities = current_mz_intensities

        current_nl_intensities = {}
        for nl_2, nl_intensity_2 in nl_data[molecule_id].items():
            for nl in list(nl_intensities.keys()):
                if abs(nl - nl_2) <= tolerance:
                    nl_intensities[nl] += nl_intensity_2
                    current_nl_intensities[nl] = nl_intensities[nl]
        nl_intensities = current_nl_intensities

    common_mz = {mz: intensity for mz, intensity in mz_intensities.items()}
    common_nl = {nl: intensity for nl, intensity in nl_intensities.items()}

    return {
        'status': 'success',
        'characteristic_ions': common_mz,
        'characteristic_nl': common_nl
    }


def parse_mgf_for_neutral_loss(file_content, top_n):
    molecules = []
    current_molecule = None
    molecule_counter = 0

    lines = file_content.splitlines()

    for line in lines:
        line = line.strip()

        if line == "BEGIN IONS":
            if current_molecule:
                molecules.append(current_molecule)
            molecule_counter += 1
            current_molecule = {
                "id": molecule_counter,
                "pepmass": "",
                "scans": "",
                "rtinseconds": "",
                "charge": "",
                "spectrum": ""
            }

        elif line.startswith("PEPMASS="):
            current_molecule["pepmass"] = line.split("=")[-1]

        elif line.startswith("SCANS="):
            current_molecule["scans"] = line.split("=")[-1]

        elif line.startswith("RTINSECONDS="):
            current_molecule["rtinseconds"] = line.split("=")[-1]

        elif line.startswith("CHARGE="):
            current_molecule["charge"] = line.split("=")[-1]

        elif line == "END IONS":
            if current_molecule:
                current_molecule[
                    "spectrum"] = f"<a href='#' onclick=\"showNLSpectrum('SCANS-{current_molecule['scans']}', '0', {top_n})\">NL Spectrum</a>"
                molecules.append(current_molecule)
                current_molecule = None

    if current_molecule:
        current_molecule[
            "spectrum"] = f"<a href='#' onclick=\"showNLSpectrum('SCANS-{current_molecule['scans']}', '0', {top_n})\">NL Spectrum</a>"
        molecules.append(current_molecule)

    return molecules


def handle_mgf_file_analysis(request, tolerance, selected_molecules=None):
    mgf_file_content = request.session.get('mgf_file_content')
    min_neutral_loss = request.session.get('minNeutralLoss')
    top_n = request.session.get('topN')

    if not mgf_file_content:
        return {'status': 'error', 'message': 'MGF file not found'}

    spectrum_data, nl_data = parse_mgf_file(mgf_file_content, min_neutral_loss, top_n)

    common_mz, common_nl = process_common_mz_nl(spectrum_data, nl_data, tolerance, selected_molecules)

    return {
        'status': 'success',
        'characteristic_ions': common_mz,
        'characteristic_nl': common_nl
    }

def parse_mgf_file(mgf_file_content, min_neutral_loss, top_n):
    spectrum_data = {}
    nl_data = {}
    molecule_count = 1

    lines = mgf_file_content.splitlines()
    for line in lines:
        line = line.strip()
        if re.match(r'^\d+\.\d+\s+\d+', line):
            mz, intensity = map(float, line.split())
            if molecule_count not in spectrum_data:
                spectrum_data[molecule_count] = []
            spectrum_data[molecule_count].append([mz, intensity])
        elif line.startswith("END IONS"):
            molecule_count += 1

    for molecule_id, data in spectrum_data.items():
        max_intensity = max([intensity for _, intensity in data])
        spectrum_data[molecule_id] = [[mz, (intensity / max_intensity) * 100] for mz, intensity in data]

        top_ions = sorted(spectrum_data[molecule_id], key=lambda x: x[1], reverse=True)[:top_n]
        neutral_losses = generate_neutral_losses(top_ions, min_neutral_loss)
        nl_percentages = calculate_neutral_loss_percentages(neutral_losses)
        nl_data[molecule_id] = nl_percentages

    return spectrum_data, nl_data


def process_common_mz_nl(spectrum_data, nl_data, tolerance, selected_molecules=None):
    mz_intensities = {}
    nl_intensities = {}

    molecules_to_process = selected_molecules if selected_molecules else list(spectrum_data.keys())

    first_molecule = molecules_to_process[0]
    for mz, intensity in spectrum_data[first_molecule]:
        mz_intensities[mz] = intensity

    for nl, nl_intensity in nl_data[first_molecule].items():
        nl_intensities[nl] = nl_intensity

    for molecule_id in molecules_to_process[1:]:
        current_mz_intensities = {}
        for mz_2, intensity_2 in spectrum_data[molecule_id]:
            for mz in list(mz_intensities.keys()):
                if abs(mz - mz_2) <= tolerance:
                    mz_intensities[mz] += intensity_2
                    current_mz_intensities[mz] = mz_intensities[mz]
        mz_intensities = current_mz_intensities

        current_nl_intensities = {}
        for nl_2, nl_intensity_2 in nl_data[molecule_id].items():
            for nl in list(nl_intensities.keys()):
                if abs(nl - nl_2) <= tolerance:
                    nl_intensities[nl] += nl_intensity_2
                    current_nl_intensities[nl] = nl_intensities[nl]
        nl_intensities = current_nl_intensities

    common_mz = {mz: intensity for mz, intensity in mz_intensities.items()}
    common_nl = {nl: intensity for nl, intensity in nl_intensities.items()}

    return common_mz, common_nl
