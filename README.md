# Zooplankton Abundance and Diversity in Seagrass Habitats – Data and R Code for Journal Submission

## 📂 Repository Overview
This repository contains the data and R scripts used for the analysis presented in the manuscript:
**"Zooplankton Abundance and Diversity in Native and Invasive Seagrass Habitats and the Implications for Juvenile Snapper Diet."**

## 📊 Data Files
- **`zooplankton_data.xlsx`**: The raw data collected from passive zooplankton traps deployed in native (*Syringodium filiforme*, *Thalassia testudinum*) and invasive (*Halophila stipulacea*) seagrass beds.
- **`metadata.csv`**: Descriptions of columns and variables in the dataset.

## 🖥️ R Scripts
- **`zooplankton_analysis.R`**: Primary script for data cleaning, statistical analysis, and visualization.
- **`zooplankton_simper.R`**: SIMPER analysis script comparing zooplankton communities across seagrass types.
- **`zooplankton_nmds.R`**: NMDS ordination analysis script.

## 📝 Usage
1. **Clone Repository:**  
   ```bash
   git clone https://github.com/GpONeal/Zooplankton_Seagrass_Study.git
   ```
2. **Open RStudio:** Load `zooplankton_analysis.R` and run line-by-line or source the script.
3. **Install Dependencies:** The scripts require packages such as `tidyverse`, `vegan`, and `ggplot2`. Install them via:
   ```r
   install.packages(c("tidyverse", "vegan", "ggplot2"))
   ```

## 📈 Results Summary
- **Zooplankton Abundance:** Higher copepod abundance in *H. stipulacea*, reduced amphipods and decapods.
- **Diversity:** No significant differences across seagrass types.
- **SIMPER Analysis:** Harpacticoid copepods were primary contributors to dissimilarities.

## 🧾 Citation
If you use this repository, please cite:
> O'Neill, G.P., Costa, S.V., & Nemeth, R.S. (2024). Zooplankton abundance and diversity in native and invasive seagrass habitats. *Journal of Marine Ecology*. DOI: [pending]

## 📜 License
This project is licensed under the MIT License.

## 🤝 Acknowledgments
- Data collection and fieldwork by Sophia V. Costa and Dr. Richard Nemeth.
- R script contributions from Dr. Kayla Blincow.
- Funding from NOAA Coral Reef Conservation Program (NA19NOS4820118) and NSF Virgin Islands EPSCoR (NSF 1946412).

---
### ✅ Ready for Journal Submission
This repository is prepared for open-access sharing as per journal data availability policies.
