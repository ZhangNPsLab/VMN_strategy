from django.shortcuts import render
from django.conf import settings
import uuid
import csv
from pyteomics import mgf
import numpy as np
from django.http import JsonResponse, HttpResponse
from SMMN.utils import module4net
import os
import pandas as pd
import shutil
import zipfile

def show_filter(request):
    if request.method == 'POST':
        ion_match_count = int(request.POST.get('ionMatchCount', 0))
        nl_match_count = int(request.POST.get('nlMatchCount', 0))
        min_normalized_intensity = float(request.POST.get('minNormalizedIntensity', 0.02))
        tolerance = float(request.POST.get('tolerance', 0.02))
        cosine_score = float(request.POST.get('cosineScore', 0.7))
        filter_model = request.POST.get('filterModel', 'MN')
        andOrvalue = int(request.POST.get('andOrValue', 0))
        characteristic_ions = request.POST.get('characteristicIon', '')
        characteristic_nl = request.POST.get('characteristicNL', '')

        common_ions = [float(ion) for ion in characteristic_ions.split() if
                       float(ion) > 0] if characteristic_ions else []
        common_neutral_losses = [float(nl) for nl in characteristic_nl.split() if
                                 float(nl) > 0] if characteristic_nl else []

        mgf_file = request.FILES.get('mgfFile', None)
        csv_file = request.FILES.get('csvFile', None) if filter_model == 'FBMN' else None

        user_directory = request.session.get('user_directory')
        if not os.path.exists(user_directory):
            user_directory = os.path.join(settings.MEDIA_ROOT, str(uuid.uuid4()))
            os.makedirs(user_directory, exist_ok=True)

        if mgf_file:
            input_mgf = os.path.join(user_directory, mgf_file.name)
            with open(input_mgf, 'wb+') as f:
                for chunk in mgf_file.chunks():
                    f.write(chunk)

            replace_feature_id_with_title(input_mgf)

            with mgf.read(input_mgf) as spectra:
                spectra_data = list(spectra)

            filtered_spectra, metadata = process_mgf_file(
                spectra_data,
                ion_match_count,
                nl_match_count,
                common_ions,
                common_neutral_losses,
                tolerance,
                min_normalized_intensity,
                andOrvalue
            )

            output_mgf = os.path.join(user_directory, 'filtered_spectra.mgf')
            output_metadata_csv = os.path.join(user_directory, 'metadata.csv')
            write_filtered_spectra(filtered_spectra, output_mgf)
            write_metadata(metadata, output_metadata_csv)

            if filter_model == 'FBMN' and csv_file:
                input_csv = os.path.join(user_directory, csv_file.name)
                with open(input_csv, 'wb+') as f:
                    for chunk in csv_file.chunks():
                        f.write(chunk)

                titles_to_keep = [spectrum['params']['title'] for spectrum in filtered_spectra if
                                  'title' in spectrum['params']]
                output_filtered_csv = os.path.join(user_directory, 'filtered_data.csv')

                filter_csv(input_csv, titles_to_keep, output_filtered_csv)

            spectra_collection = module4net.load_mgf_file(output_mgf)
            cosine_score_threshold = cosine_score
            peak_tolerance = tolerance
            top_k = 10

            all_matches = module4net.generate_all_matches(spectra_collection, peak_tolerance, cosine_score_threshold,
                                                          top_k)

            csv_filename = module4net.match_to_csv(all_matches)

            output_csv = os.path.join(user_directory, 'match.csv')
            shutil.copy(csv_filename, output_csv)
            G = module4net.draw_network(csv_filename, cosine_score, component_size=5, peak_matching_rate=0.0,
                                        structure_mz=0)
            network_html = module4net.draw_interactive_network_with_communities(G, k=10)

            return render(request, 'network.html', {'network_html': network_html})

        else:
            return JsonResponse({'status': 'error', 'message': 'MGF file is required for MN or FBMN filter.'},
                                status=400)

    return JsonResponse({'status': 'error', 'message': 'Invalid request method.'}, status=405)

def filter_csv(input_csv, titles_to_keep, output_csv):
    df = pd.read_csv(input_csv)
    titles_to_keep = list(map(str, titles_to_keep))
    df['row ID'] = df['row ID'].astype(str)
    filtered_df = df[df['row ID'].isin(titles_to_keep)]
    filtered_df.to_csv(output_csv, index=False)


def replace_feature_id_with_title(input_file):
    temp_file = input_file + ".temp"
    with open(input_file, 'r') as infile, open(temp_file, 'w') as outfile:
        for line in infile:
            outfile.write(line.replace("FEATURE_ID=", "TITLE="))
    os.replace(temp_file, input_file)


def process_mgf_file(spectra_data, ion_threshold, neutral_loss_threshold, common_ions, common_neutral_losses, tolerance,
                     min_normalized_intensity, andOrvalue):
    filtered_spectra = []
    metadata = []

    for spectrum in spectra_data:
        ion_score = calculate_ion_score(spectrum, common_ions, tolerance, min_normalized_intensity)
        neutral_loss_score = calculate_neutral_loss_score(spectrum, common_neutral_losses, tolerance,
                                                          min_normalized_intensity)
        spectrum['ion_score'] = ion_score
        spectrum['neutral_loss_score'] = neutral_loss_score

        classification = 0

        if andOrvalue == 0:
            if ion_score >= ion_threshold and neutral_loss_score >= neutral_loss_threshold:
                classification = 3
        elif andOrvalue == 1:
            if ion_score >= ion_threshold and neutral_loss_score >= neutral_loss_threshold:
                classification = 3
            elif ion_score >= ion_threshold:
                classification = 1
            elif neutral_loss_score >= neutral_loss_threshold:
                classification = 2

        if classification > 0:
            filtered_spectra.append(spectrum)
            metadata.append({
                'row ID': spectrum['params'].get('title', ''),
                'row m/z': spectrum['params'].get('pepmass', [None])[0],
                'row retention time': spectrum['params'].get('rtinseconds', 0) / 60,
                'classification': classification
            })

    return filtered_spectra, metadata


def calculate_ion_score(spectrum, common_ions, tolerance, min_normalized_intensity):
    ion_score = 0
    intensity_array = spectrum['intensity array']
    mz_array = spectrum['m/z array']
    params = spectrum.get('params', {})
    pepmass = params.get('pepmass', [0])[0]

    if len(intensity_array) == 0:
        return ion_score

    valid_intensities = [intensity for intensity, mz in zip(intensity_array, mz_array) if abs(mz - pepmass) >= 0]
    max_intensity = max(valid_intensities) if valid_intensities else 0
    min_intensity = np.min(intensity_array)

    if max_intensity == min_intensity:
        normalized_intensity = np.zeros_like(intensity_array)
    else:
        normalized_intensity = (intensity_array - min_intensity) / (max_intensity - min_intensity)

    mz_list = [fragment_ion for y, fragment_ion in enumerate(spectrum['m/z array']) if
               normalized_intensity[y] > min_normalized_intensity]

    matched_ions = set()

    for mz in mz_list:
        for ion in common_ions:
            if abs(mz - ion) < tolerance and ion not in matched_ions:
                ion_score += 1
                matched_ions.add(ion)

    return ion_score


def calculate_neutral_loss_score(spectrum, common_neutral_losses, tolerance, min_normalized_intensity):
    neutral_loss_score = 0
    intensity_array = spectrum['intensity array']
    mz_array = spectrum['m/z array']
    params = spectrum.get('params', {})
    pepmass = params.get('pepmass', [0])[0]

    if len(intensity_array) == 0:
        return neutral_loss_score

    valid_intensities = [intensity for intensity, mz in zip(intensity_array, mz_array) if abs(mz - pepmass) >= 0]
    sorted_intensities = sorted(valid_intensities, reverse=True)
    max_intensity = sorted_intensities[0] if sorted_intensities else 0
    min_intensity = np.min(intensity_array)

    if max_intensity == min_intensity:
        normalized_intensity = np.zeros_like(intensity_array)
    else:
        normalized_intensity = (intensity_array - min_intensity) / (max_intensity - min_intensity)

    mz_list = [fragment_ion for y, fragment_ion in enumerate(spectrum['m/z array']) if
               normalized_intensity[y] > min_normalized_intensity]

    neutral_losses = [(mz1 - mz2, mz1, mz2) for i, mz1 in enumerate(mz_list) for j, mz2 in enumerate(mz_list) if
                      i != j]

    matched_losses = set()

    for nl, mz1, mz2 in neutral_losses:
        for common_nl in common_neutral_losses:
            if abs(nl - common_nl) < tolerance and common_nl not in matched_losses:
                neutral_loss_score += 1
                matched_losses.add(common_nl)
    return neutral_loss_score


def write_filtered_spectra(filtered_spectra, output_mgf):
    with open(output_mgf, 'w') as f:
        mgf.write(filtered_spectra, output=f)


def write_metadata(metadata, output_csv):
    keys = ['row ID', 'row m/z', 'row retention time', 'classification']
    with open(output_csv, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=keys)
        writer.writeheader()
        writer.writerows(metadata)


def download_filter_data(request):
    user_directory = request.session.get('user_directory')
    if not user_directory:
        return HttpResponse("User directory not found", status=404)

    filtered_data_csv = os.path.join(user_directory, 'filtered_data.csv')
    filtered_spectra_mgf = os.path.join(user_directory, 'filtered_spectra.mgf')
    metadata_csv = os.path.join(user_directory, 'metadata.csv')

    zip_filename = os.path.join(user_directory, 'filtered_data.zip')

    with zipfile.ZipFile(zip_filename, 'w') as zip_file:
        if os.path.exists(filtered_data_csv):
            zip_file.write(filtered_data_csv, 'filtered_data.csv')
        if os.path.exists(filtered_spectra_mgf):
            zip_file.write(filtered_spectra_mgf, 'filtered_spectra.mgf')
        if os.path.exists(metadata_csv):
            zip_file.write(metadata_csv, 'metadata.csv')

    if os.path.exists(zip_filename):
        with open(zip_filename, 'rb') as f:
            response = HttpResponse(f.read(), content_type='application/zip')
            response['Content-Disposition'] = f'attachment; filename="filtered_data.zip"'
            return response
    else:
        return HttpResponse("Error creating ZIP file", status=500)