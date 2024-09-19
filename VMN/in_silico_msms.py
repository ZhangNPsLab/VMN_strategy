import re
from django.conf import settings
import uuid
from django.http import JsonResponse, HttpResponse, FileResponse
import os
from VMN.tasks import run_simulation_task

def simulate_data(request):
    if request.method == 'POST':
        smiles_list = []
        uploaded_file = request.FILES.get('file')
        if uploaded_file:
            for line in uploaded_file:
                decoded_line = line.decode('utf-8').strip()
                if decoded_line:
                    smiles_list.append(decoded_line)

        user_input = request.POST.get('smiles', '').strip()
        if user_input:
            smiles_list.extend(re.split(r'\s+', user_input))

        if not smiles_list:
            return JsonResponse({'status': 'error', 'message': 'No SMILES provided.'}, status=400)

        user_directory = os.path.join(settings.MEDIA_ROOT, 'temp', str(uuid.uuid4()))
        os.makedirs(user_directory, exist_ok=True)

        request.session['user_directory'] = user_directory

        molecule_file_path = os.path.join(user_directory, "molecule.txt")
        with open(molecule_file_path, 'w') as f:
            for idx, smiles in enumerate(smiles_list, start=1):
                f.write(f"Molecule{idx} {smiles}\n")

        # 启动异步任务
        task = run_simulation_task.delay(molecule_file_path, user_directory)

        # 返回任务ID
        return JsonResponse({'status': 'success', 'task_id': task.id})


def check_task_status(request, task_id):
    print(f"Checking task status for task_id: {task_id}")
    task_result = run_simulation_task.AsyncResult(task_id)
    print(f"Task state: {task_result.state}")

    if task_result.state == 'PENDING':
        response = {'state': task_result.state, 'status': 'Pending...'}
    elif task_result.state == 'SUCCESS':
        result = task_result.result
        if result['status'] == 'SUCCESS':
            file_path = result['result']
            print(f"Task completed. Output file path: {file_path}")
            parsed_result = parse_output_log(file_path)
            print(f"Parsed result: {parsed_result}")
            response = {
                'state': task_result.state,
                'molecules': parsed_result
            }
        else:
            response = {
                'state': 'FAILURE',
                'status': result.get('error', 'Unknown error')
            }
    elif task_result.state == 'FAILURE':
        response = {
            'state': task_result.state,
            'status': str(task_result.info),
        }
    else:
        response = {'state': task_result.state, 'status': 'Processing...'}

    return JsonResponse(response)

def parse_output_log(file_path):
    with open(file_path, 'r') as file:
        content = file.read()

    molecules = []
    current_molecule = None

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
                    "10eV": f"<a href='#' onclick=\"showSpectrum('{molecule_id}', '0')\">Spectrum</a>",
                    "20eV": f"<a href='#' onclick=\"showSpectrum('{molecule_id}', '1')\">Spectrum</a>",
                    "40eV": f"<a href='#' onclick=\"showSpectrum('{molecule_id}', '2')\">Spectrum</a>"
                }
            }
        elif line.startswith("#Formula="):
            current_molecule["formula"] = line.split("=")[-1]
        elif line.startswith("#PMass="):
            current_molecule["pepmass"] = line.split("=")[-1]

    if current_molecule:
        molecules.append(current_molecule)

    return molecules

def show_spectrum(request):
    molecule_id = request.GET.get('molecule_id')
    energy_level = request.GET.get('energy_level')

    energy_mapping = {
        '0': '10eV',
        '1': '20eV',
        '2': '40eV'
    }

    user_directory = request.session.get('user_directory')
    if not user_directory:
        return JsonResponse({'status': 'error', 'message': 'User directory not found'}, status=404)

    output_file = os.path.join(user_directory, "output.log")
    if not os.path.exists(output_file):
        return JsonResponse({'status': 'error', 'message': 'Output log file not found'}, status=404)

    # 解析 output.log 文件，提取相应的质谱数据
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

    if not spectrum_data or not pepmass:
        return JsonResponse({'status': 'error', 'message': 'Spectrum data not found'}, status=404)

    return JsonResponse({
        'status': 'success',
        'spectrum_data': spectrum_data,
        'pepmass': pepmass,
        'energy_level': target_energy
    })


def download_output_log(request):
    user_directory = request.session.get('user_directory')
    if not user_directory:
        return HttpResponse('User directory not found', status=404)

    output_file = os.path.join(user_directory, "output.log")
    if not os.path.exists(output_file):
        return HttpResponse('File not found', status=404)

    response = FileResponse(open(output_file, 'rb'), content_type='text/plain')
    response['Content-Disposition'] = 'attachment; filename="output.log"'
    return response
