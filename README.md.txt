# 🧬 Lung Cancer Gene Expression Analysis

This project focuses on identifying differentially expressed genes (DEGs) in **lung cancer** using bioinformatics techniques.  
The analysis includes data processing, visualization (heatmap and volcano plot), and protein–protein interaction (PPI) network analysis.

---

## 📊 Project Workflow

1. **Data Collection**  
   - Gene expression dataset was downloaded from the **GEO (Gene Expression Omnibus)** database.

2. **Differential Expression Analysis**  
   - DEGs were identified using statistical thresholds:  
     - |log2FC| ≥ 1  
     - adjusted p-value ≤ 0.05  
   - Results were visualized using:
     - **Volcano Plot** — shows upregulated (red) and downregulated (blue) genes.
     - **Heatmap** — displays the top 30 most significant DEGs.

3. **Protein–Protein Interaction (PPI) Network**  
   - Significant DEGs were imported into **STRING** and visualized using **Cytoscape**.  
   - Hub genes were identified based on network connectivity.

---

## 🔬 Key Findings

- **Upregulated genes (red)** are highly expressed in tumor samples, suggesting potential roles in cancer progression.  
- **Downregulated genes (blue)** show reduced expression and might act as tumor suppressors.  
- Hub genes such as **SPP1**, **COL10A1**, **CTHRC1**, and **GOLM1** may serve as:
  - 🎯 **Therapeutic targets**
  - 🧫 **Potential biomarkers** for lung cancer detection or prognosis

---

## 🧰 Tools and Packages

- **R (Bioconductor)**
  - `limma`
  - `ggplot2`
  - `pheatmap`
- **STRING database** – for protein–protein interactions  
- **Cytoscape** – for network visualization  

---

## 🖼️ Visual Results

### Volcano Plot
![Volcano Plot](volcano plot for genes affect on lung cancer)

### Heatmap (Top 30 Genes)
![Heatmap](Heatmap 30 genes in lung cancer)

### PPI Network
![PPI Network](generatetaskspecificdownloadfile)

---

## 📁 Files Included

| File | Description |
|------|--------------|
| `analysis_code.R` | Full R script used for DEG analysis and plotting |
| `DEGs.csv` | List of significant differentially expressed genes |
| `Heatmap.png` | Top 30 DEGs visualized as a heatmap |
| `VolcanoPlot.png` | Volcano plot showing gene regulation |
| `PPI.png` | Protein–Protein Interaction network from STRING |

---

## 👨‍🔬 Author
**Sameh Mohamed**  
Biotechnology Student — Zagazig University  
Passionate about Bioinformatics and Molecular Biology  

📎 [LinkedIn](https://www.linkedin.com/in/sameh-mohamed-82437629a/) | [GitHub](foudas400-blip)

---

## 🧠 Future Work
- Perform **functional enrichment analysis (GO/KEGG)** for DEGs  
- Validate hub genes experimentally or via additional datasets  
- Integrate survival analysis to confirm prognostic biomarkers

---

⭐ *If you found this project useful, please star the repository!*
