import re
import os
from decimal import Decimal, ROUND_DOWN
from django.http import JsonResponse, HttpResponse

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

def generate_neutral_loss(request):
    if request.method == 'POST':
        use_first_step_output = request.POST.get('useFirstStepOutput', False)
        top_n = int(request.POST.get('topN'))
        min_neutral_loss = float(request.POST.get('minNeutralLoss'))
        request.session['minNeutralLoss'] = min_neutral_loss
        request.session['topN'] = top_n

        if use_first_step_output:
            user_directory = request.session.get('user_directory')
            if not user_directory:
                return JsonResponse({'status': 'error', 'message': 'User directory not found'}, status=404)

            output_file = os.path.join(user_directory, "output.log")

            if not os.path.exists(output_file):
                return JsonResponse({'status': 'error', 'message': 'Output log file not found'}, status=404)

            neutral_loss_links = parse_output_log_for_neutral_loss(output_file, top_n)

        else:
            uploaded_file = request.FILES.get('file')
            if not uploaded_file:
                return JsonResponse({'status': 'error', 'message': 'No file uploaded'}, status=400)

            file_content = uploaded_file.read().decode('utf-8')
            request.session['mgf_file_content'] = file_content

            neutral_loss_links = parse_mgf_for_neutral_loss(file_content, top_n)

        return JsonResponse({
            'status': 'success',
            'message': 'Neutral loss spectrum generated successfully.',
            'data': neutral_loss_links
        })

    return JsonResponse({'status': 'error', 'message': 'Invalid request method.'}, status=405)


def parse_output_log_for_neutral_loss(file_path, top_n):
    with open(file_path, 'r') as file:
        content = file.read()

    molecules = []
    current_molecule = None
    current_energy = None

    for line in content.splitlines():
        line = line.strip()
        if line.startswith("#ID="):
            if current_molecule:
                molecules.append(current_molecule)
            molecule_id = line.split("=")[-1]
            current_molecule = {
                "id": molecule_id,
                "formula": "",
                "pepmass": "",
                "spectra": {
                    "10eV": f"<a href='#' onclick=\"showNLSpectrum('{molecule_id}', '0', {top_n})\">Spectrum</a>",
                    "20eV": f"<a href='#' onclick=\"showNLSpectrum('{molecule_id}', '1', {top_n})\">Spectrum</a>",
                    "40eV": f"<a href='#' onclick=\"showNLSpectrum('{molecule_id}', '2', {top_n})\">Spectrum</a>"
                }
            }
        elif line.startswith("#Formula="):
            current_molecule["formula"] = line.split("=")[-1]
        elif line.startswith("#PMass="):
            current_molecule["pepmass"] = line.split("=")[-1]

    if current_molecule:
        molecules.append(current_molecule)

    return molecules


def generate_neutral_losses(top_ions, min_neutral_loss, charge=1):
    neutral_losses = {}
    for i in range(len(top_ions)):
        mz1, intensity1 = top_ions[i]
        for j in range(len(top_ions)):
            if i != j:
                mz2, intensity2 = top_ions[j]
                neutral_loss = (Decimal(mz1) - Decimal(mz2)) * Decimal(charge)

                if neutral_loss >= Decimal(min_neutral_loss):
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


def calculate_neutral_loss_percentages(neutral_losses):
    if not neutral_losses or all(len(nl_list) == 0 for nl_list in neutral_losses.values()):
        return {}

    nl_percentages = {}

    total_intensities = {nl: sum(intensity for _, _, intensity in nl_list) for nl, nl_list in neutral_losses.items()}

    max_intensity = max(total_intensities.values())

    for nl, total_intensity in total_intensities.items():
        nl_percentages[float(nl)] = (total_intensity / max_intensity) * 100

    return nl_percentages


def show_nl_spectrum(request):
    molecule_id = request.GET.get('molecule_id')
    energy_level = request.GET.get('energy_level', None)
    user_directory = request.session.get('user_directory')
    mgf_file_content = request.session.get('mgf_file_content')
    top_n = request.session.get('topN')
    min_neutral_loss = request.session.get('minNeutralLoss')

    if not user_directory and not mgf_file_content:
        return JsonResponse({'status': 'error', 'message': 'File not found'}, status=404)

    session_data_key = None

    if molecule_id.startswith("SCANS-"):
        scan_number = molecule_id.split("-")[-1]

        spectrum_data = []
        collecting_data = False
        charge = 1
        pepmass = None

        lines = mgf_file_content.splitlines()
        for line in lines:
            line = line.strip()

            if line.startswith("PEPMASS=") and pepmass is None:
                pepmass = float(line.split("=")[-1].split()[0])

            if line.startswith(f"SCANS={scan_number}"):
                collecting_data = True

            elif collecting_data:
                if line.startswith("CHARGE="):
                    charge_str = line.split("=")[-1].replace('+', '').strip()
                    charge = int(charge_str) if charge_str.isdigit() else 1
                elif line and line[0].isdigit():
                    mz, intensity = map(float, line.split())
                    spectrum_data.append([mz, intensity])
                elif line == "END IONS":
                    break

        if not spectrum_data or pepmass is None:
            return JsonResponse({'status': 'error', 'message': 'Spectrum data not found'}, status=404)

        max_intensity = max([intensity for _, intensity in spectrum_data])
        spectrum_data = [[mz, (intensity / max_intensity) * 100] for mz, intensity in spectrum_data]

        top_ions = sorted(spectrum_data, key=lambda x: x[1], reverse=True)[:top_n]
        neutral_losses = generate_neutral_losses(top_ions, min_neutral_loss, charge)
        nl_percentages = calculate_neutral_loss_percentages(neutral_losses)

    else:
        energy_mapping = {
            '0': '10eV',
            '1': '20eV',
            '2': '40eV'
        }

        output_file = os.path.join(user_directory, "output.log")
        if not os.path.exists(output_file):
            return JsonResponse({'status': 'error', 'message': 'Output log file not found'}, status=404)

        with open(output_file, 'r') as file:
            content = file.read()

        spectrum_data = []
        current_energy = None
        collecting_data = False
        pepmass = None
        target_energy = energy_mapping.get(energy_level)
        for line in content.splitlines():
            line = line.strip()

            if line.startswith(f"#ID={molecule_id}"):
                collecting_data = True
                continue

            if collecting_data:
                if line.startswith("energy0"):
                    current_energy = "10eV"
                elif line.startswith("energy1"):
                    current_energy = "20eV"
                elif line.startswith("energy2"):
                    current_energy = "40eV"

                elif line.startswith("#PMass="):
                    pepmass = float(line.split("=")[-1])

                elif current_energy == target_energy and collecting_data:
                    if line:
                        mz, intensity = map(float, line.split())
                        spectrum_data.append([mz, intensity])
                    elif line == '':
                        break
                elif line == "":
                    collecting_data = False

        if not spectrum_data or pepmass is None:
            return JsonResponse({'status': 'error', 'message': 'Spectrum data not found'}, status=404)

        max_intensity = max([intensity for _, intensity in spectrum_data])
        spectrum_data = [[mz, (intensity / max_intensity) * 100] for mz, intensity in spectrum_data]

        top_ions = sorted(spectrum_data, key=lambda x: x[1], reverse=True)[:top_n]
        neutral_losses = generate_neutral_losses(top_ions, min_neutral_loss)
        nl_percentages = calculate_neutral_loss_percentages(neutral_losses)

    return JsonResponse({
        'status': 'success',
        'spectrum_data': spectrum_data,
        'nl_data': nl_percentages,
        'pepmass': pepmass,
        'session_key': session_data_key
    })


def download_neutral_loss_log(request):
    user_directory = request.session.get('user_directory')
    top_n = request.session.get('topN')
    min_neutral_loss = request.session.get('minNeutralLoss')

    if not user_directory or not top_n or not min_neutral_loss:
        return HttpResponse("Missing session data", status=400)

    output_log_file = os.path.join(user_directory, 'output.log')
    if not os.path.exists(output_log_file):
        return HttpResponse("Output log file not found", status=404)

    with open(output_log_file, 'r') as file:
        lines = file.readlines()

    output_lines = []
    current_molecule_header = []
    spectrum_data = {}
    current_energy = None

    def process_current_molecule(header, spectrum_data):
        result = []
        result.extend(header)

        result.append(f"#neutral losses data calculated by VMN")
        result.append(f"#top N = {top_n} Min Neutral Loss = {min_neutral_loss}")

        for energy, data in spectrum_data.items():
            if not data:
                continue

            result.append(energy)
            max_intensity = max([intensity for _, intensity in data])
            normalized_spectrum = [[mz, (intensity / max_intensity) * 100] for mz, intensity in data]
            top_ions = sorted(normalized_spectrum, key=lambda x: x[1], reverse=True)[:top_n]

            neutral_losses = generate_neutral_losses(top_ions, min_neutral_loss)
            nl_percentages = calculate_neutral_loss_percentages(neutral_losses)

            sorted_nl_percentages = sorted(nl_percentages.items(), key=lambda x: x[1], reverse=True)

            for nl, intensity in sorted_nl_percentages:
                result.append(f"{nl:.5f} {intensity:.5f}")

        result.append("")
        return result

    for line in lines:
        stripped_line = line.strip()

        if stripped_line.startswith("#ID="):
            if current_molecule_header:
                output_lines.extend(process_current_molecule(current_molecule_header, spectrum_data))

            current_molecule_header = []
            spectrum_data = {"energy0": [], "energy1": [], "energy2": []}
            current_energy = None
            current_molecule_header.append(stripped_line)

        elif stripped_line.startswith(("#SMILES", "#InChiKey", "#Formula", "#PMass")):
            current_molecule_header.append(stripped_line)

        elif stripped_line.startswith("energy"):
            current_energy = stripped_line
            spectrum_data[current_energy] = []

        elif re.match(r'^\d', stripped_line):
            if current_energy:
                mz, intensity = map(float, stripped_line.split())
                spectrum_data[current_energy].append([mz, intensity])

    if current_molecule_header:
        output_lines.extend(process_current_molecule(current_molecule_header, spectrum_data))

    response = HttpResponse("\n".join(output_lines), content_type='text/plain')
    response['Content-Disposition'] = 'attachment; filename="neutral_loss.log"'
    return response


def download_neutral_loss_mgf(request):
    mgf_file_content = request.session.get('mgf_file_content')
    top_n = request.session.get('topN')
    min_neutral_loss = request.session.get('minNeutralLoss')

    if not mgf_file_content or top_n is None or min_neutral_loss is None:
        return HttpResponse("Missing required data in session", status=400)

    output_lines = []
    molecule_count = 0
    spectrum_data = []

    lines = mgf_file_content.splitlines()
    for i, line in enumerate(lines):
        line = line.strip()

        if line.startswith("BEGIN IONS"):
            molecule_count += 1
            spectrum_data = []
            output_lines.append(line)

        elif line.startswith(("TITLE", "PEPMASS", "SCANS", "RTINSECONDS", "CHARGE", "MSLEVEL", "MERGED_STATS")):
            output_lines.append(line)

        elif re.match(r'^\d+\.\d+\s+\d+', line):
            try:
                mz, intensity = map(float, line.split())
                spectrum_data.append([mz, intensity])
            except ValueError:
                print(f"Skipping invalid line: {line}")

        elif line == "END IONS":
            if not spectrum_data:
                print(f"No spectrum data found for molecule {molecule_count}. Skipping.")
                continue

            try:
                max_intensity = max([intensity for _, intensity in spectrum_data])
                normalized_spectrum = [[mz, (intensity / max_intensity) * 100] for mz, intensity in spectrum_data]
                top_ions = sorted(normalized_spectrum, key=lambda x: x[1], reverse=True)[:top_n]

                neutral_losses = generate_neutral_losses(top_ions, min_neutral_loss)
                nl_percentages = calculate_neutral_loss_percentages(neutral_losses)

                sorted_nl_percentages = sorted(nl_percentages.items(), key=lambda x: x[1], reverse=True)

                output_lines.append("#neutral losses data calculated by VMN")
                output_lines.append(f"#top N = {top_n} Min Neutral Loss = {min_neutral_loss}")
                for nl, intensity in sorted_nl_percentages:
                    output_lines.append(f"{nl:.5f} {intensity:.5f}")
                output_lines.append(line)
                output_lines.append("")
            except ValueError as e:
                print(f"Error processing molecule {molecule_count}: {e}")
                continue

    response = HttpResponse("\n".join(output_lines), content_type='text/plain')
    response['Content-Disposition'] = 'attachment; filename="neutral_loss.mgf"'
    return response