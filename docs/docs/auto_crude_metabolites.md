# Auto Crude Metabolites Filter

The **Auto Crude Metabolites Filter** module helps researchers filter crude metabolite data using **characteristic ion** and **characteristic NL** data to build the target metabolite network.

---

## Parameter Inputs

To use this module, users need to input the following parameters:

### 1. **Characteristic Ion and Characteristic NL**
- Users can input either **characteristic ion** or **characteristic NL** (both cannot be empty).
- Multiple entries should be separated by spaces.
- The module will filter the crude metabolites based on the provided characteristic data.

### 2. **Filter Model**
This module supports two filter models:
- **MN Model**: Only a `.mgf` file of the crude metabolites is required.
- **FBMN Model**: Both a `.mgf` file and a corresponding `.csv` file (processed using tools like MZmine) need to be uploaded.

### 3. **Match Count**
- Define the required number of matches with the characteristic data.
- Molecules whose characteristic data matches count higher than the user-defined threshold will be filtered.
  
Two match modes are available:
- **AND Mode**: Both the **characteristic ion** and **characteristic NL** must match the threshold.
- **OR Mode**: Either the **characteristic ion** or **characteristic NL** must match the threshold.

> **Note:**  
> If users only want to filter based on **characteristic ion**, they should still select the **AND** mode. Selecting **OR** mode would result in all molecules being filtered as every molecule satisfies `characteristic NL >= 0`. The same logic applies for **characteristic NL** filtering.

### 4. **Minimum Normalized Intensity**
- Set a minimum normalized intensity threshold.
- Ions with intensities below this value will be filtered out to eliminate low-intensity noise ions.

### 5. **Tolerance**
- Define a tolerance range within which ions or neutral losses will be considered as matching with **characteristic ion** or **characteristic NL**.
- The default value is `0.02`.

### 6. **Cosine Score**
- Set a **cosine similarity threshold**.
- Compounds with cosine similarity higher than the threshold will be connected, and a molecular network preview will be generated for the filtered compounds.

---

## Viewing and Downloading the Results

Filtered target compounds will be displayed as a **molecular network** in the "Preview of Filter MN" section on the right. Users can:

- Preview the filtered compounds' **m/z**, **retention time (rt)**, **cosine score**, and other information.
- Download the filtered `.mgf` and `.csv` files for further analysis or upload to the **GNPS platform** to generate the final **MZmol**.

---

## Metadata File and Classification

Upon downloading, a `metadata.csv` file will be generated containing a "classification" column. The numbers in this column indicate the specific classification of the molecule:

- **Number 3**: Both the **characteristic NL** and **characteristic ion** match counts exceed the threshold.
- **Number 2**: Only the **characteristic NL** match count exceeds the threshold.
- **Number 1**: Only the **characteristic ion** match count exceeds the threshold.

---

### Quick Tips:
- Use the **AND Mode** to filter based on both characteristic ion and neutral loss.
- Set appropriate **cosine score thresholds** to ensure the best network generation results.
- Preview the molecular network to quickly identify key connections and data relationships.

This module allows for advanced filtering of crude metabolites, providing researchers with a streamlined approach to building targeted metabolite networks.
