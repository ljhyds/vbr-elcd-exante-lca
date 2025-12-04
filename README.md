# VBR Ex-Ante LCA with ELCD 3.2 using openLCA 2.5 and Brightway

![Python](https://img.shields.io/badge/Python-3.11-blue.svg)
![openLCA](https://img.shields.io/badge/openLCA-2.5.x-2ea44f.svg)
![Brightway](https://img.shields.io/badge/Brightway-4.x-orange.svg)

Reproducible free-toolchain workflow for the ex-ante LCA of Vitrified Bauxite Residue (VBR).

---

## 1. Overview

This repository provides a fully reproducible Python workflow to perform an ex-ante LCA of **Vitrified Bauxite Residue (VBR)** using only free tools:

* **openLCA 2.5** (via IPC for accessing ELCD)
* **ELCD 3.2 background database** (obtained through openLCA Nexus)
* **Brightway** (for foreground modelling and LCIA computation)

The workflow consists of two scripts:

* `import_elcd_lcia.py` — imports ELCD background processes and LCIA methods into Brightway using openLCA IPC
* `run_vbr_lca.py` — builds the VBR foreground system and performs LCIA

This enables an open and reproducible reconstruction of the core modelling logic used in recent ex-ante LCA studies on vitrification-based emerging materials.

---

## 2. Requirements

### Software Requirements

* **openLCA 2.5.x**
* **Python 3.11** (recommended)
* **Brightway2 / Brightway3 stack**
* **openLCA IPC client libraries (olca-ipc, olca-schema)**

Install Python dependencies:

```bash
pip install -r requirements.txt
```

---

## 3. Getting the ELCD 3.2 Database

The ELCD 3.2 database (originally released by the European Commission / JRC) was obtained via the **openLCA Nexus** platform.
The database is **not redistributed** in this repository due to licensing restrictions.

Users must download the ELCD 3.2 database in **openLCA database format** (as provided on the openLCA Nexus platform) and import it manually into **openLCA 2.5** before running the scripts.

---

## 4. openLCA 2.5 Setup (IPC)

1. Launch **openLCA 2.5**
2. Import the **ELCD 3.2** ILCD database
3. Activate the database
4. Start the IPC Server:

   * `Tools → Developer Tools → IPC Server`
   * Set port = **8080**

The importer script expects openLCA to be reachable at:

```
localhost:8080
```

---

## 5. Step A — Import ELCD Background & LCIA into Brightway

Run:

```bash
python import_elcd_lcia.py
```

This script will:

* Construct a **UUID‑preserving biosphere** from ELCD elementary flows
* Import ELCD **background processes** via openLCA IPC
* Import LCIA methods (e.g., ReCiPe 2016 Midpoint (H), Global warming) with **UUID‑preserved LCIA–flow matching** to avoid mismatches
* Preserve **quantitative reference amounts** (not forced to 1.0)

After completion, the Brightway project will contain:

* **an extended biosphere3**, populated with UUID‑preserved elementary flows imported from openLCA

* `elcd_bw_ipc` (background database)

* LCIA methods imported from `ELCD`

## 6. Step B — Build VBR Foreground & Run LCIA

Run:

```bash
python run_vbr_lca.py
```

This script:

* Builds foreground unit processes for:

  * mixing
  * pelletising
  * vitrification furnace
  * cooling
  * milling
* Links these processes to ELCD electricity and materials
* Adds direct CO₂ emissions for the vitrification furnace
* Runs LCIA using the imported ELCD-based ReCiPe method
* Produces impact results for the VBR functional unit (1 kg)

---

## 7. Reproducibility Notes

Tested versions:

* **openLCA:** 2.5.x
* **Python:** 3.11.x
* **Brightway:**

  * bw2data ≥ 4.0
  * bw2calc ≥ 2.0
  * bw2io ≥ 0.9
  * bw_processing ≥ 0.9
  * matrix_utils ≥ 0.3
* **IPC Libraries:**

  * olca-ipc ≥ 2.0
  * olca-schema ≥ 2.0

If version conflicts occur, pin exact versions from `requirements.txt`.

---

## 8. License & Data Policy

All Python code in this repository is released under the **MIT License**.

The ELCD database and related LCIA content remain under their original **European Commission / JRC licenses**.
They **must be downloaded from the official openLCA Nexus platform** and are **not included** in this repository.

---

## 9. Citation

If you use or extend this workflow, please cite:

* Georgiades, A., Myers, R. J., et al. (2025). *Life cycle assessment of vitrified bauxite residue across laboratory, pilot, and industrial scales.* The International Journal of Life Cycle Assessment. 【Supplementary Information Sections S2–S5 were used for furnace modelling, Q_min calculation, grinding energy assumptions, raw material formulation, proxy technologies, and electricity scenarios.】

* ELCD 3.2 dataset (European Commission / Joint Research Centre). Obtained via openLCA Nexus.

* openLCA 2.5 — GreenDelta GmbH.

* Mutel, C. (2017). *Brightway: An open-source LCA framework for Python.* Journal of Open Source Software.

* Bond, F. C. (1952) and Charles, R. J. (1957) — Classical milling energy equations used for grinding energy calculations..

---

## 10. Parameter Sources (Physical & Thermodynamic Constants)

The physical constants used in the foreground model fall into two groups:

* **Standard chemical/thermodynamic constants** (stoichiometric CO₂ factors, molecular weights, specific heat of water, latent heat of vaporization)

  * Taken from general reference data (IUPAC Atomic Weights, NIST Chemistry WebBook).

* **Model-specific ranges for vitrified bauxite residue** (specific heat capacity of solid oxides/slag precursors, vitrification furnace temperature, etc.)

  * Based on Georgiades et al. (2025) main text and Supplementary Information (e.g. Sections S2–S5).

* **Grinding energy relationship**

  * Based on the classical Bond (1952) and Charles (1957) milling energy equations.

* **Industrial furnace efficiency ranges**

  * General **industrial furnace efficiency ranges** (≈75–82%) are informed by Best Available Techniques (BAT) reference documents for high-temperature industrial furnace operations.

These references correspond directly to the parameters implemented in `run_vbr_lca.py`, ensuring that all physical and thermodynamic constants used in the model are clearly traceable to their underlying sources.

---
