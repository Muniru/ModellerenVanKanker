# Plan van Aanpak: Medicijnhergebruik bij Kanker (Python Focus)

## **Doelstelling**
Het identificeren van bestaande medicijnen die hergebruikt kunnen worden voor de behandeling van specifieke kankertypes door middel van bioinformatica-analyse van genexpressie- en medicijndata, uitgevoerd met Python.

---

## **Werkwijze**
Het project wordt in vijf fasen uitgevoerd: dataverzameling, datavoorbereiding, analyse, interpretatie en rapportage.

---

## **Fasen en Taken**

### **Week 1: Voorbereiding en Dataverzameling**
- **Stap 1.1: Literatuuronderzoek**
  - Onderzoek het concept van medicijnhergebruik, en ontdek hoe genexpressiedata en interactienetwerken kunnen worden gebruikt.
  - Bronnen: PubMed, Google Scholar, en reviewartikelen.

- **Stap 1.2: Datasetselectie**
  - **Doel:** Verzamel RNA-seq data en medicijndata.
  - **Acties:**
    1. Zoek naar RNA-seq datasets via [GEO](https://www.ncbi.nlm.nih.gov/geo/) of [TCGA](https://portal.gdc.cancer.gov/).
    2. Download medicijndata via [LINCS L1000](https://lincsproject.org/LINCS/) of [DrugBank](https://go.drugbank.com/).
  - **Output:** Genormaliseerde RNA-seq dataset en medicijnkenmerken.

- **Stap 1.3: Toolselectie**
  - Python-pakketten:
    - **Voor data-analyse:** `pandas`, `numpy`, `matplotlib`, `seaborn`, `scipy`.
    - **Voor RNA-seq:** `scanpy`, `sklearn`, `statsmodels`.
    - **Voor netwerkvisualisatie:** `networkx`, `pyvis`.
  - Installeer en test alle benodigde tools.

---

### **Week 2: Data-analyse**
- **Stap 2.1: Data preprocessing**
  - **RNA-seq data:**
    - Lees en controleer de dataset (bijv. met `pandas`).
    - Normaliseer de count data naar TPM of log-transformaties.
    - Gebruik `scanpy` voor filtering en kwaliteitscontrole.
  - **Medicijndata:**
    - Controleer op duplicaten, ontbrekende waarden en standaardiseer eigenschappen.

- **Stap 2.2: Differentiële genexpressieanalyse**
  - Identificeer significant verschillend tot expressie gebrachte genen (DEGs) tussen tumor- en normaalweefsel.
  - Gebruik een statistische test zoals een t-test of een log-likelihood test (`scipy.stats` of `statsmodels`).

- **Stap 2.3: Medicijnmatching**
  - Gebruik medicijndata (bijv. van LINCS of DrugBank) om te zoeken naar medicijnen die DEGs beïnvloeden.
  - Correlaties en netwerkanalyse uitvoeren tussen medicijnen en genen.

---

### **Week 3: Validatie en Rapportage**
- **Stap 3.1: Gen-Medicijn Netwerken**
  - Maak een interactienetwerk van genen en hun beïnvloedende medicijnen met behulp van `networkx` of `pyvis`.

- **Stap 3.2: Hypothesevorming**
  - Formuleer hypothesen over de potentiële mechanismen van hergebruikte medicijnen.
  - Controleer in de literatuur of de medicijnen al gerelateerd zijn aan kankeronderzoek.

- **Stap 3.3: Rapportage**
  - Schrijf een rapport en maak visualisaties:
    - **Heatmaps:** Genexpressie (gebruik `seaborn.heatmap`).
    - **Netwerken:** Medicijn-gen interacties (gebruik `networkx`).
  - Structuur:
    - **Introductie**
    - **Methodologie**
    - **Resultaten**
    - **Discussie**
    - **Conclusie**

---

## **Tijdsplanning**

| Week      | Activiteit                              | Output                                     |
|-----------|----------------------------------------|-------------------------------------------|
| Week 1    | Literatuur, data en tools verzamelen   | Dataset, tools en kennisbasis             |
| Week 2    | Data preprocessing en analyse          | Lijst met DEGs en kandidaat-medicijnen    |
| Week 3    | Validatie, netwerkconstructie, rapport | Rapport en visualisaties                  |

---

## **Verwachte Resultaten**
- Een lijst van bestaande medicijnen die mogelijk hergebruikt kunnen worden voor de behandeling van kanker.
- Analyse met Python-script en visualisaties voor validatie.



---
## Datasets

- Count data: [https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE62944](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE62944)
- 