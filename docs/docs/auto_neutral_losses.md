# Auto Neutral Losses Spectrum

The **Auto Neutral Losses Spectrum** module is designed to help users obtain neutral losses data, allowing better utilization of the **neutral loss characteristics** of compounds.

---

## Data Sources

When generating a neutral losses spectrum, users have two options:

- Use the **in-silico MS/MS data** generated in the first module.
- Import **experimentally collected compound MS/MS data** in the form of a `.mgf` file.

---

## Parameter Inputs

To generate the neutral losses spectrum, users need to specify the following two parameters:

### 1. **Top N (positive integer)**
- This parameter selects the **Top N fragment ions** with the highest intensity to generate the neutral loss spectrum.
- This helps to **reduce noise interference** in the data.
- Typically, Top N can be set to `30`.

### 2. **Min Neutral Loss (positive number)**
- This parameter sets the **minimum neutral loss mass** that will be accepted in the spectrum.
- A typical value for Min Neutral Loss is `50`, as lower-mass neutral losses (such as O, H₂O, CH₂, CH₃, and C₂H₄) are less characteristic and frequently found in natural products.

> **Note:**  
> Selecting appropriate values for these parameters ensures better quality results and less noise in the neutral loss spectrum.

---

## Understanding the Generated Data

In the generated **neutral losses data**:
- **Neutral loss** represents the **mass difference** between precursor fragment ions and product fragment ions.
- **Intensity** is calculated as the **average intensity** of the precursor and product ions.

If multiple precursor-product ion pairs generate the same neutral loss, the intensity is the **sum** of those ion pair intensities. The final intensity is **normalized to 100**, with the maximum intensity set at 100.

---

## Viewing and Downloading the Data

Once the neutral losses spectrum is generated, the following options are available:

- The **generated spectrum** is displayed as **hyperlinks** in the results table on the right side of the screen.
- By clicking on the **Spectrum** link in the table, users can:
  - View the corresponding **neutral losses spectrum**.
  - Explore the **MS/MS spectrum** for more detailed interaction with the data.

Additionally, users can download the generated neutral losses data by clicking on the **Download button** on the right side of the table.

---

### Quick Tips:
- Use higher values for **Top N** to include more fragments and gain a broader overview of the spectrum.
- Ensure **Min Neutral Loss** is set appropriately to exclude common low-mass neutral losses, which may not provide useful insight.

This module makes it easier for users to analyze and interpret neutral losses, allowing for enhanced data filtering and exploration.
