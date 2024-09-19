# Auto Characteristic Data

The **Auto Characteristic Data** module is designed to help users automatically retrieve **characteristic ions** and **neutral losses** that are commonly present in selected molecules. This data provides valuable insights into the **core structure** of compounds and their **fragmentation patterns**.

---

## Steps to Use the Module

Follow these steps to retrieve characteristic data for the selected molecules:

### 1. **Select Molecules for Analysis**
- Users can select **"all"** molecules or **"select"** specific molecules by inputting their IDs.
  
### 2. **Choose Energy Level**
- Select the energy level for the analysis:
  - **10 eV**
  - **20 eV**
  - **40 eV**
- For **experimental data**, users should select `"expt"` in the energy level option.

### 3. **Define Tolerance Range**
- This parameter defines the **mass difference** within which fragment ions or neutral losses are considered identical.
  - For **in-silico MS/MS data**, the tolerance can typically be set around **10E-5**.
  - For **experimental MS/MS data**, the tolerance should be set to **0.02** or lower for higher precision.

---

## Result Table Overview

The results will be displayed in a table with the following columns:

- **Characteristic Ion**:  
  Lists the fragment ions commonly found in the selected molecules.
  
- **Characteristic NL**:  
  Displays the neutral losses that are shared among the selected molecules.
  
- **Average Intensity**:  
  Shows the **sum of the intensities** of the characteristic ions or neutral losses across the selected molecules.

---

## Interacting with the Data

Once the characteristic data is generated:

- Users can interact with the **result table** to view and analyze the data.
- **Download the Data**:  
  By clicking the **Download button**, users can download the characteristic data in a convenient format for further analysis.

---

### Quick Tips:
- Use a **lower tolerance range** for experimental data to ensure higher precision.
- Select the appropriate energy level to match the experimental conditions or simulation parameters.

The **Auto Characteristic Data** module makes it easy to automatically identify key fragmentation patterns and neutral losses, offering deeper insights into molecular structures.
