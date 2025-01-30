import pandas as pd
import math
import networkx as nx
from pyvis.network import Network
import bisect
import re
from collections import namedtuple
import os


class Spectrum:
    def __init__(self, filename, rt, scan, unknown_param, peaks, mz, charge, some_other_param):
        self.filename = filename
        self.rt = rt
        self.scan = scan
        self.unknown_param = unknown_param
        self.peaks = peaks
        self.mz = mz
        self.charge = charge
        self.some_other_param = some_other_param


Match = namedtuple('Match', ['peak1', 'peak2', 'score'])
Peak = namedtuple('Peak', ['mz', 'intensity'])
Alignment = namedtuple('Alignment', ['peak1', 'peak2'])


def generate_spectrum_network(mgf_file, cosine_score, peak_tolerance, top_k=10, component_size=5,
                              peak_matching_rate=0, structure_mz=0, k=10.0):
    spectra_collection = load_mgf_file(mgf_file)

    all_matches = generate_all_matches(spectra_collection, peak_tolerance, cosine_score, top_k)

    csv_filename = match_to_csv(all_matches)

    G = draw_network(csv_filename, component_size, peak_matching_rate, structure_mz)

    network_html = draw_interactive_network_with_communities(G, k)

    return network_html




def load_mgf_file(filename):
    spectra = []
    current_spectrum = None

    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()

            if line == 'BEGIN IONS':
                current_spectrum = {
                    'peaks': [],
                    'mz': 0,
                    'charge': 0,
                    'rt': 0,
                    'scan': 0,
                    'filename': filename  # 保存当前文件名
                }

            elif line == 'END IONS':
                if current_spectrum and len(current_spectrum['peaks']) > 0:
                    spectra.append(current_spectrum)

            elif line.startswith('PEPMASS='):
                current_spectrum['mz'] = float(line.split('=')[1])

            elif line.startswith('CHARGE='):
                charge_str = line.split('=')[1].replace('+', '').replace('-', '')
                current_spectrum['charge'] = int(charge_str)

            elif line.startswith('RTINSECONDS='):
                current_spectrum['rt'] = float(line.split('=')[1])

            elif line.startswith('SCANS='):
                current_spectrum['scan'] = int(line.split('=')[1])

            elif re.match(r'\d+\.\d+ \d+', line):
                mz, intensity = map(float, line.split())
                current_spectrum['peaks'].append((mz, intensity))

    return spectra

def generate_all_matches(spectra_collection, peak_tolerance, cosine_score_threshold, top_k):
    all_matches = []

    for i, base_spectrum in enumerate(spectra_collection):
        print('base_spectrum', base_spectrum)
        match_list = []

        for j, spectrum in enumerate(spectra_collection):
            if i != j:
                cosine_score, matched_peaks = score_alignment(base_spectrum, spectrum, peak_tolerance)
                if cosine_score >= cosine_score_threshold:
                    match_obj = {}
                    match_obj["filename"] = base_spectrum['filename']
                    match_obj["scan"] = base_spectrum['scan']
                    match_obj["queryfilename"] = spectrum['filename']
                    match_obj["queryscan"] = spectrum['scan']
                    match_obj["mz1"] = spectrum['mz']
                    match_obj["rt1"] = spectrum['rt']
                    match_obj["mz2"] = base_spectrum['mz']
                    match_obj["rt2"] = base_spectrum['rt']
                    match_obj["cosine"] = cosine_score
                    match_obj["matchedpeaks"] = matched_peaks
                    match_obj["mzerror"] = abs(base_spectrum['mz'] - spectrum['mz'])
                    match_obj['intensity'] = sum(peak[1] for peak in spectrum['peaks'])
                    match_obj['source'] = "classical_molecular"
                    match_obj["Peak-matching Rate"] = None
                    match_list.append(match_obj)

        if len(match_list) == 0:
            default_match = {
                "filename": None,
                "scan": None,
                "queryfilename": base_spectrum['filename'],
                "queryscan": base_spectrum['scan'],
                "mz1": base_spectrum['mz'],
                "rt1": base_spectrum['rt'],
                "mz2": None,
                "rt2": None,
                "cosine": 0,
                "matchedpeaks": 0,
                "mzerror": None,
                "intensity": sum(peak[1] for peak in base_spectrum['peaks']),
                "source": "classical_molecular",
                "Peak-matching Rate": None
            }
            match_list.append(default_match)

        match_list = sorted(match_list, key=lambda x: x['cosine'], reverse=True)[:top_k]
        all_matches.extend(match_list)

    return all_matches

def score_alignment(spectrum1, spectrum2, tolerance):
    spec1_peaks = convert_to_peaks(spectrum1['peaks'])
    spec2_peaks = convert_to_peaks(spectrum2['peaks'])
    pm1 = spectrum1['mz']
    pm2 = spectrum2['mz']

    return calculate_alignment(spec1_peaks, spec2_peaks, pm1, pm2, tolerance)

def calculate_alignment(spec1, spec2, pm1, pm2, tolerance, max_charge_consideration=1):
    if len(spec1) == 0 or len(spec2) == 0:
        return 0.0, []

    spec1_n = sqrt_normalize_spectrum(spec1)
    spec2_n = sqrt_normalize_spectrum(spec2)

    shift = pm1 - pm2
    zero_shift_alignments = find_match_peaks_efficient(spec1_n, spec2_n, 0, tolerance)
    real_shift_alignments = find_match_peaks_efficient(spec1_n, spec2_n, shift, tolerance) if abs(
        shift) > tolerance else []

    if max_charge_consideration > 1:
        for charge_considered in range(2, max_charge_consideration + 1):
            real_shift_alignments += find_match_peaks_efficient(spec1_n, spec2_n, shift / charge_considered, tolerance)

    real_shift_alignments = list(set(real_shift_alignments))

    zero_shift_match = [alignment_to_match(spec1_n, spec2_n, alignment) for alignment in zero_shift_alignments]
    real_shift_match = [alignment_to_match(spec1_n, spec2_n, alignment) for alignment in real_shift_alignments]

    all_possible_match_scores = zero_shift_match + real_shift_match
    all_possible_match_scores.sort(key=lambda x: x.score, reverse=True)

    total_score = 0.0
    reported_alignments = []
    spec1_peak_used = set()
    spec2_peak_used = set()

    for match in all_possible_match_scores:
        if match.peak1 not in spec1_peak_used and match.peak2 not in spec2_peak_used:
            spec1_peak_used.add(match.peak1)
            spec2_peak_used.add(match.peak2)
            reported_alignments.append(match)
            total_score += match.score

    return total_score, reported_alignments


def sqrt_normalize_spectrum(spectrum):
    acc_norm = sum([s.intensity for s in spectrum])
    normed_value = math.sqrt(acc_norm)
    return [Peak(s.mz, math.sqrt(s.intensity) / normed_value) for s in spectrum]


def find_match_peaks_efficient(spec1, spec2, shift, tolerance):
    adj_tolerance = tolerance + 0.000001
    spec2_mass_list = [peak.mz for peak in spec2]
    alignment_mapping = []

    for i, peak in enumerate(spec1):
        left_bound = peak.mz - shift - adj_tolerance
        right_bound = peak.mz - shift + adj_tolerance
        left_index = bisect.bisect_left(spec2_mass_list, left_bound)
        right_index = bisect.bisect_right(spec2_mass_list, right_bound)

        for j in range(left_index, right_index):
            alignment_mapping.append(Alignment(i, j))

    return alignment_mapping


def alignment_to_match(spec1_n, spec2_n, alignment):
    s1_peak = spec1_n[alignment.peak1].intensity
    s2_peak = spec2_n[alignment.peak2].intensity
    match_score = s1_peak * s2_peak
    return Match(alignment.peak1, alignment.peak2, match_score)

def convert_to_peaks(peak_tuples):
    return [Peak(*p) for p in peak_tuples]

def match_to_csv(all_matches):
    csv_filename = 'match.csv'
    df = pd.DataFrame(all_matches)
    df.columns = ['Filename', 'CLUSTERID2', 'Query Filename', 'CLUSTERID1', "mz1", "rt1", "mz2", "rt2", 'Cosine',
                  'Matched Peaks',
                  'DeltaMZ',
                  'explained_intensity', 'source', "Peak-matching Rate"]
    df_selected = df[
        ['CLUSTERID1', 'CLUSTERID2', 'Cosine', "mz1", "rt1", "mz2", "rt2", 'DeltaMZ', 'explained_intensity', 'source',
         "Peak-matching Rate"]]
    df_selected.to_csv(csv_filename, index=False)
    return csv_filename


# 绘制网络图
def draw_network(csv_filename, cosine_threshold, component_size=5, peak_matching_rate=0.0, structure_mz=0):
    df = pd.read_csv(csv_filename)

    G = nx.MultiGraph()

    for i in range(len(df)):
        node1 = standardize_node(df["CLUSTERID1"].iloc[i])
        node2 = standardize_node(df["CLUSTERID2"].iloc[i])
        source = df["source"].iloc[i]

        # 根据来源设置颜色
        color = "#87CEFA" if source == "classical_molecular" else "#CD5C5C"
        edge_color = "#87CEFA" if source == "classical_molecular" else "#CD5C5C"
        edge_style = "solid" if source == "classical_molecular" else "dashed"
        EdgeType = "classical_molecular" if source == "classical_molecular" else "csmn"

        valid_node1 = is_valid_node(node1)
        valid_node2 = is_valid_node(node2)

        mass_difference = df["DeltaMZ"].iloc[i]
        retention_time1 = df["rt1"].iloc[i]
        mz1 = df["mz1"].iloc[i]
        retention_time2 = df["rt2"].iloc[i]
        mz2 = df["mz2"].iloc[i]

        cosine_score = float(df["Cosine"].iloc[i] or 0)

        # 添加节点和边，根据来源设置属性
        if valid_node1:
            if not G.has_node(node1):
                color = "#CD5C5C" if float(mz1) > structure_mz else color
                G.add_node(node1, retention_time=retention_time1, mz=mz1, color=color, shape="dot")

        if valid_node2:
            if not G.has_node(node2):
                color = "#CD5C5C" if float(mz2) > structure_mz else color
                G.add_node(node2, retention_time=retention_time2, mz=mz2, color=color, shape="dot")

        if valid_node1 and valid_node2:
            if cosine_score >= cosine_threshold:
                if not G.has_edge(node1, node2) and not G.has_edge(node2, node1):
                    G.add_edge(node1, node2, mass_difference=mass_difference, cosine_score=cosine_score,
                               component=-1, EdgeType=EdgeType, EdgeScore=cosine_score, style=edge_style, color=edge_color)

    connected_components = list(nx.connected_components(G))
    for component in connected_components:
        if len(component) > component_size:
            edges_to_remove = [(u, v) for u, v, d in G.edges(component, data=True) if 'cosine_score' in d]
            edges_to_remove.sort(key=lambda x: x[2]['cosine_score'])
            while len(component) > component_size and edges_to_remove:
                edge_with_lowest_score = edges_to_remove.pop(0)
                G.remove_edge(edge_with_lowest_score[0], edge_with_lowest_score[1])
                component = max(nx.connected_components(G), key=len)

    # 保存图形
    nx.write_graphml(G, "ClassicalNetwork.graphml")
    return G



def is_valid_node(node):
    if pd.isna(node) or node == "":
        return False
    try:
        float(node)
        return True
    except ValueError:
        return False



def standardize_node(node):
    try:
        return str(int(float(node)))
    except ValueError:
        return node


# 绘制交互式网络
def draw_interactive_network_with_communities(G, k):
    net = Network(height="680px", width="100%", notebook=True)
    for node, attrs in G.nodes(data=True):
        net.add_node(str(node), label=str(node), title=f"MZ: {attrs['mz']} RT: {attrs['retention_time']}s", color=attrs['color'])

    for edge in G.edges(data=True):
        net.add_edge(str(edge[0]), str(edge[1]), title=f"Cosine Score: {edge[2]['cosine_score']}",
                     color=edge[2]['color'])

    return net.generate_html()


def parse_table_with_headers(filename, skip_incomplete_lines=False, debug=False, delimiter="\t"):
    input_file = None
    try:
        input_file = open(filename, "r", encoding='ascii', errors='ignore')
        print(f"Successfully opened file: {filename}")
    except FileNotFoundError:
        print(f"File not found: {filename}")
        project_root = os.path.abspath(os.path.dirname(__file__))
        full_filename = os.path.join(project_root, filename)
        print(f"Trying with full path: {full_filename}")
        try:
            input_file = open(full_filename, "r", encoding='ascii', errors='ignore')
        except FileNotFoundError:
            print(f"File still not found: {full_filename}")
            return -1, None


def parse_table_without_headers(filename):
    input_file = open(filename, "r")

    line_count = 0
    column_values = {}
    for line in input_file:
        line_splits = line.rstrip().split("\t")

        line_count += 1
        if line_count == 1:
            for i in range(len(line_splits)):
                column_values[i] = []

        column_count = 0
        for line_split in line_splits:
            column_values[column_count].append(line_split)
            column_count += 1

    return (line_count - 1, column_values)
