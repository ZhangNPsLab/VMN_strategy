# VMN Data Modules

MZmol's main functions encompass four key aspects:

1. **In-silico MS/MS data**
2. **Automated generation of neutral loss spectra**
3. **Automated capture of characteristic mass spectrometry data**
4. **Automated filtering of target metabolite networks**

Each module can be used independently or in combination to aid researchers in generating and analyzing mass spectrometry data efficiently.

---

## 1. In silico MS/MS Data of Virtual Molecules

- Users can **input the SMILES** of virtual molecules directly, with multiple SMILES separated by spaces.
- Alternatively, users can **upload a .txt file** containing SMILES for batch generation of in silico MS/MS data.

This module allows for the simulation of MS/MS spectra for molecules that cannot be easily analyzed experimentally.

---

## 2. Auto Neutral Losses Spectrum

- Users can generate **neutral losses data** using the in-silico MS/MS data produced by the first module.
- Alternatively, users can **import experimental MS/MS data** in the form of a `.mgf` file to generate corresponding neutral losses spectra.

This helps users quickly generate and analyze neutral losses from their compounds.

---

## 3. Auto Characteristic Data

- No data import is required in this module.
- Users can select the **neutral losses data** generated in the second module to automatically produce **characteristic mass spectrometry data**.

This module captures common fragmentation patterns for detailed analysis.

---

## 4. Auto Crude Metabolites Filter

The data required for this module depends on the selected filtering mode:

- **MN Mode**: Only a `.mgf` file of crude metabolites is required.
- **FBMN Mode**: Users need both a `.mgf` file of crude metabolites and the corresponding **CSV file** (which can be processed using software like MZmine).

This module is designed to filter and isolate targeted metabolite networks from complex datasets.

---

**Note:**  
Each module is flexible, allowing for independent usage or integration with other modules to provide comprehensive analysis capabilities in a virtual molecular network approach.
