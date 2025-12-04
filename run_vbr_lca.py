from collections import defaultdict
from bw2data import projects, Database, Method, databases, methods, get_activity
from bw2calc import LCA
import matplotlib.pyplot as plt

# -----------------------------
# 1. Project & background setup
# -----------------------------

PROJECT_NAME = "vbr_project"
BG_DB_NAME = "elcd_bw_ipc"
FG_DB_NAME = "vbr_fg_processes"

def set_project():
    projects.set_current(PROJECT_NAME)

def get_db(name):
    if name not in databases:
        raise RuntimeError(
            f"Database '{name}' not found. Check that it is imported in project '{PROJECT_NAME}'."
        )
    return Database(name)

# Georgiades et al. 2025 industrial-scale upscaling scenarios (2020 baseline)
# Source: Georgiades et al. 2025, SI Table S21 (climate change, kg CO2-eq per kg VBR)
MYERS_SCENARIOS_2020 = {
    "The reproduced model": None,
    "M3.1.1 (process modelling)": {
        "Other upstream": 0.05,
        "Furnace electricity": 0.21,
        "Furnace heat (natural gas)": 0.0,
        "Furnace heat (other fuels)": 0.0,
        "Furnace direct": 0.21,
    },
    "M3.2.1 (linear regression)": {
        "Other upstream": 0.05,
        "Furnace electricity": 0.19,
        "Furnace heat (natural gas)": 0.0,
        "Furnace heat (other fuels)": 0.0,
        "Furnace direct": 0.21,
    },
    "M3.2.2 (power regression)": {
        "Other upstream": 0.05,
        "Furnace electricity": 0.20,
        "Furnace heat (natural gas)": 0.0,
        "Furnace heat (other fuels)": 0.0,
        "Furnace direct": 0.21,
    },
    "M3.3.1 (clinker proxy – Centre et al.)": {
        "Other upstream": 0.07,
        "Furnace electricity": 0.03,
        "Furnace heat (natural gas)": 0.0,
        "Furnace heat (other fuels)": 0.18,
        "Furnace direct": 0.21,
    },
    "M3.3.2 (clinker proxy – ecoinvent)": {
        "Other upstream": 0.07,
        "Furnace electricity": 0.05,
        "Furnace heat (natural gas)": 0.0,
        "Furnace heat (other fuels)": 0.38,
        "Furnace direct": 0.21,
    },
    "M3.3.3 (ceramic proxy – BAT)": {
        "Other upstream": 0.10,
        "Furnace electricity": 0.03,
        "Furnace heat (natural gas)": 0.23,
        "Furnace heat (other fuels)": 0.0,
        "Furnace direct": 0.21,
    },
    "M3.3.4 (ceramic proxy – ecoinvent)": {
        "Other upstream": 0.13,
        "Furnace electricity": 0.12,
        "Furnace heat (natural gas)": 0.33,
        "Furnace heat (other fuels)": 0.0,
        "Furnace direct": 0.21,
    },
}

# -----------------------------------
# 2. Find background & biosphere data
# -----------------------------------

def find_bg_activity(db, name_substring, unit=None, location=None):
    """
    Return first "real" background activity whose name contains substring
    (and optional unit/location filters).

    - Skips any activity whose name contains 'dummy' (e.g. 'Dummy_cullet float glass').
    - Skips activities that have no self-production exchange (to avoid broken processes).
    """
    name_substring = name_substring.lower()
    for act in db:
        name = act['name'].lower()
        if "dummy" in name:
            continue
        if name_substring not in name:
            continue
        if unit and act['unit'] != unit:
            continue
        if location and act.get('location') != location:
            continue
        prods = [
            exc for exc in act.exchanges()
            if exc['type'] == 'production' and exc['input'] == act.key
        ]
        if len(prods) == 0:
            continue
        return act

    raise RuntimeError(
        f"Real activity containing '{name_substring}' not found in DB '{db.name}'. "
        "Checked: skipped 'dummy' names and activities with no production exchanges."
    )

# -----------------------------
# 3. Parameters & furnace model
# -----------------------------

def compute_Qmin_and_Efurn(p):
    # 1) Moisture content and wet BR mass
    m_BR_wet = 1.0 / (1.0 - p["w_moist"])    # kg wet BR per kg dry BR
    m_water = m_BR_wet * p["w_moist"]        # kg water per kg dry BR solids

    # === Fe2O3–C reduction stoichiometry ===
    # Amount of Fe2O3 contained in 1 kg of dry BR solids
    m_Fe2O3 = p["m_solid"] * p["w_Fe2O3"]            # kg Fe2O3 per batch
    n_Fe2O3 = m_Fe2O3 / p["M_Fe2O3"]                 # mol Fe2O3

    # Moles of carbon available from coke
    n_C = p["m_coke"] / p["M_C"]                     # mol C

    # Reaction: Fe2O3 + C → 2 FeO + CO
    # Stoichiometry requires 1 mol C per 1 mol Fe2O3
    n_red = min(n_Fe2O3, n_C)                        # mol of Fe2O3 that is reduced (C-limited)

    # Endothermic enthalpy for Fe2O3 → FeO reduction (+, heat required)
    Q_red = n_red * p["dH_red_Fe2O3_to_FeO"]         # kJ

    # Mass loss due to oxygen removal (per mol Fe2O3 reduced)
    m_O_loss_solid = n_red * p["m_loss_per_mol_Fe2O3"]

    # 2) Solid mass to be preheated (BR + CaCO3 + SiO2)
    m_pre = p["m_solid"] + p["m_CaCO3"] + p["m_silica"]

    # === Final VBR solid mass after Fe2O3 reduction ===
    m_solid_reduced = p["m_solid"] - m_O_loss_solid
    m_vbr = m_solid_reduced + p["m_silica"] + 0.56 * p["m_CaCO3"]

    # 3) Latent heat of evaporation and sensible heat
    Q_latent = m_water * p["h_vap"]

    T_start = p.get("T_in_preheat", p["T_amb"])
    Q_sensible = m_pre * p["cp_solid"] * (p["T_furn"] - T_start)

    # Heat required for CaCO3 decomposition
    n_CaCO3 = p["m_CaCO3"] / p["M_CaCO3"]
    Q_CaCO3 = n_CaCO3 * p["dH_CaCO3"]

    # 4) Fusion heat (optional)
    if p.get("ignore_fusion", False):
        Q_fusion = 0.0
    else:
        Q_fusion = p["L_fus"] * m_vbr

    # 5) Enthalpy of coke combustion: C → CO2
    #    All coke carbon is assumed to oxidise to CO2 eventually.
    n_C_to_CO2 = n_C                                 # mol C → CO2
    Q_C_comb = n_C_to_CO2 * p["dH_C_to_CO2"]         # kJ (negative, exothermic)

    # 6) Minimum furnace energy demand (enthalpy-based)
    #    Process heat = sensible + latent + CaCO3 + fusion + Fe-oxide reduction + C combustion
    Q_min_batch = Q_sensible + Q_latent + Q_CaCO3 + Q_fusion + Q_red + Q_C_comb
    Q_min_batch_kWh = Q_min_batch / 3600.0

    # Convert to per-kg-VBR basis (functional unit)
    Q_min_kWh_per_kg = Q_min_batch_kWh / m_vbr
    E_furn_ind_per_kg = Q_min_kWh_per_kg / p["eta_furn_ind"]

    return {
        "m_BR_wet": m_BR_wet,
        "m_water": m_water,
        "m_vbr": m_vbr,
        "Q_min_kWh": Q_min_kWh_per_kg,
        "E_furn_ind": E_furn_ind_per_kg,
    }


def get_default_params():
    params = {
        "w_moist": 0.09,           # Initial moisture content of BR (mass fraction)
        "cp_solid": 1.0,           # Specific heat capacity of solids (kJ/kg/K)
        "T_amb": 25.0,             # Ambient temperature (°C)
        "T_furn": 1200.0,          # Target furnace temperature (°C)

        # Inlet temperature of preheated solids
        "T_in_preheat": 300.0,     # °C

        # Batch basis: 1 kg dry BR solids
        "m_solid": 1.1,            # kg dry BR solids per batch

        # === Fe2O3–C reduction parameters ===
        "w_Fe2O3": 0.48,           # Fe2O3 weight fraction in BR solids
        "M_Fe2O3": 0.1597,         # Molar mass of Fe2O3 (kg/mol)
        "M_C": 0.0120,             # Molar mass of carbon (kg/mol)
        "dH_red_Fe2O3_to_FeO": 170.0,   # Enthalpy for Fe2O3 → FeO (kJ/mol)
        "m_loss_per_mol_Fe2O3": 0.016,  # Solid mass lost per mol Fe2O3 reduced (kg)

        "h_vap": 2257.0,           # Latent heat of water evaporation (kJ/kg)

        "m_CaCO3": 0.22,           # CaCO3 addition (kg per kg BR solids)
        "M_CaCO3": 0.100,          # Molar mass of CaCO3 (kg/mol)
        "dH_CaCO3": 180.0,         # Decomposition enthalpy of CaCO3 (kJ/mol)

        "m_silica": 0.15,          # Silica addition (kg/kg)

        # Coke input – used in stoichiometry and C → CO2 enthalpy
        "m_coke": 0.0204,          # Coke input (kg/kg)
        "LHV_coke": 28200.0,       # Kept for reference; not used in Q_min
        "dH_C_to_CO2": -393.0,     # Enthalpy of C(s) + O2 → CO2 (kJ/mol, exothermic)

        "L_fus": 500.0,            # Heat of fusion of slag (kJ/kg)

        "ignore_fusion": True,     # Ignore fusion heat (set Q_fusion = 0)
        "eta_furn_ind": 0.75,      # Induction furnace efficiency

        # Electricity consumption for upstream steps
        "E_mix": 0.025,            # Mixing (kWh/kg)
        "E_pellet": 0.05,          # Pelletisation (kWh/kg)
        "E_preheat_extra": 0.00,   # Additional preheating electricity (kWh/kg)
        "E_cool": 0.00,            # Cooling electricity (kWh/kg)
        "E_mill": 0.05,            # Final milling (kWh/kg)

        "transport_distance_km": 50.0,   # Transport distance for BR (km)
    }

    res = compute_Qmin_and_Efurn(params)
    params["m_vbr"] = res["m_vbr"]
    params["Q_min_kWh"] = res["Q_min_kWh"]
    params["E_furn"] = res["E_furn_ind"]
    return params




# -------------------------------
# 4. Create foreground DB
# -------------------------------

PRIMARY_LCIA_METHOD = ('ELCD_LCIA', 'ReCiPe 2016 Midpoint (H)', 'Global warming')
CO2_ONLY_METHOD = ('VBR_INTERNAL', 'CO2-only', 'climate change')

def ensure_co2_only_method():
    """Register CO2-only method: CF=1 for all CO2 flows in biosphere3."""
    if CO2_ONLY_METHOD in methods:
        return
    bio = get_db("biosphere3")
    cfs = []
    for flow in bio:
        if "carbon dioxide" in flow['name'].lower() and flow['unit'] == "kilogram":
            cfs.append((flow.key, 1.0))
    if not cfs:
        raise RuntimeError("No CO2 flows found for CO2-only method.")
    m = Method(CO2_ONLY_METHOD)
    m.register()
    m.write(cfs)

def pick_co2_for_direct(el_act, method_key):
    bio = get_db("biosphere3")
    m = Method(method_key)
    try:
        cf_dict = dict(m.load())
    except Exception:
        cf_dict = {}

    candidates = []
    for exc in el_act.exchanges():
        if exc["type"] != "biosphere":
            continue
        flow = exc.input
        name = flow.get("name", "").lower()
        if "carbon dioxide" not in name:
            continue
        if flow.get("unit") != "kilogram":
            continue

        fk = flow.key
        cf = cf_dict.get(fk, None)
        if cf is None:
            continue
        try:
            cf = float(cf)
        except Exception:
            continue
        if cf <= 0:
            continue

        candidates.append((abs(cf - 1.0), cf, flow))

    if candidates:
        candidates.sort(key=lambda x: (x[0], -x[1]))
        _, cf, flow = candidates[0]
        print("[CO2-PICK] Using electricity-mix CO2 with positive CF:")
        print(f"  name={flow['name']} | key={flow.key} | CF={cf}")
        return flow

    for flow in bio:
        if "carbon dioxide" not in flow["name"].lower():
            continue
        if flow.get("unit") != "kilogram":
            continue
        fk = flow.key
        cf = cf_dict.get(fk, None)
        if cf is None:
            continue
        try:
            cf = float(cf)
        except Exception:
            continue
        if cf <= 0:
            continue
        print("[CO2-PICK] Fallback to biosphere3 CO2 with positive CF:")
        print(f"  name={flow['name']} | key={flow.key} | CF={cf}")
        return flow

    raise RuntimeError("No CO2 flow with positive CF found in method.")

def build_foreground_db(params):
    bg = get_db(BG_DB_NAME)

    el_eu = find_bg_activity(
        bg,
        "Electricity Mix, consumption mix, at consumer, AC, < 1kV",
        location="EU-27",
    )
    limestone = find_bg_activity(bg, "Calcium carbonate > 63 microns, at plant")
    silica = find_bg_activity(bg, "Very fine milled silica sand d50 = 20 micrometer")
    coal = find_bg_activity(bg, "Hard Coal, consumption mix, at power plant")
    lorry = find_bg_activity(bg, "Lorry transport, Euro 0, 1, 2, 3, 4 mix, 22 t total weight")

    bio_co2 = pick_co2_for_direct(el_eu, PRIMARY_LCIA_METHOD)

    def kWh_to_MJ(kwh):
        # IMPORTANT: ELCD electricity behaves as kWh-equivalent → no 3.6 scaling
        return kwh * 3.6

    ton_km_per_kg = params["transport_distance_km"] / 1000.0
    total_ton_km = ton_km_per_kg * (params["m_CaCO3"] + params["m_silica"] + params["m_coke"])

    m_vbr = params["m_vbr"]

    co2_per_kg_coke = 44.0 / 12.0
    direct_co2_coke = (params["m_coke"] * co2_per_kg_coke) / m_vbr

    co2_per_kg_caco3 = 44.0 / 100.0
    direct_co2_caco3 = (params["m_CaCO3"] * co2_per_kg_caco3) / m_vbr

    fg_data = {}

    def add_fg_activity(code, name, exchanges):
        fg_data[(FG_DB_NAME, code)] = {
            "name": name,
            "unit": "kilogram",
            "location": "GLO",
            "reference product": name,
            "exchanges": [
                {"input": (FG_DB_NAME, code), "amount": 1.0, "type": "production"},
            ] + exchanges,
        }

    add_fg_activity("mixing", "Mixing for VBR, 1 kg", [
        {"input": el_eu.key, "amount": kWh_to_MJ(params["E_mix"]), "unit": "MJ", "type": "technosphere"},
    ])

    add_fg_activity("pelletising", "Pelletising for VBR, 1 kg", [
        {"input": el_eu.key, "amount": kWh_to_MJ(params["E_pellet"]), "unit": "MJ", "type": "technosphere"},
    ])

    furnace_exchanges = [
        {"input": el_eu.key, "amount": kWh_to_MJ(params["E_furn"] + params["E_preheat_extra"]), "unit": "kWh", "type": "technosphere"},
        {"input": limestone.key, "amount": params["m_CaCO3"], "unit": "kilogram", "type": "technosphere"},
        {"input": silica.key, "amount": params["m_silica"], "unit": "kilogram", "type": "technosphere"},
        {"input": coal.key, "amount": params["m_coke"], "unit": "kilogram", "type": "technosphere"},
        {"input": bio_co2.key, "amount": direct_co2_caco3, "unit": "kilogram", "type": "biosphere"},
        {"input": bio_co2.key, "amount": direct_co2_coke, "unit": "kilogram", "type": "biosphere"},
    ]
    add_fg_activity("furnace", "Furnace & high-T processing for VBR, 1 kg", furnace_exchanges)

    add_fg_activity("cooling", "Cooling for VBR, 1 kg", [
        {"input": el_eu.key, "amount": kWh_to_MJ(params["E_cool"]), "unit": "MJ", "type": "technosphere"},
    ])

    add_fg_activity("milling", "Milling for VBR, 1 kg", [
        {"input": el_eu.key, "amount": kWh_to_MJ(params["E_mill"]), "unit": "MJ", "type": "technosphere"},
    ])

    add_fg_activity("transport_raw_materials",
        f"Transport of raw materials ({params['transport_distance_km']} km) for VBR, 1 kg", [
        {"input": lorry.key, "amount": total_ton_km, "unit": "ton kilometer", "type": "technosphere"},
    ])

    fg_db = Database(FG_DB_NAME)
    fg_db.write(fg_data)

    print("\n[INFO] Selected background activities:")
    print(f"  Electricity: {el_eu['name']} | unit={el_eu['unit']} | loc={el_eu.get('location')}")
    print(f"  Limestone:   {limestone['name']} | unit={limestone['unit']}")
    print(f"  Silica:      {silica['name']} | unit={silica['unit']}")
    print(f"  Coke/Coal:   {coal['name']} | unit={coal['unit']}")
    print(f"  Lorry:       {lorry['name']} | unit={lorry['unit']}")
    print("")
    print(f"[INFO] Approx. solid VBR mass per batch (for info): {params['m_vbr']:.3f} kg")

    return fg_db


def choose_lcia_method():
    if PRIMARY_LCIA_METHOD in methods:
        print(f"[INFO] Using imported LCIA method: {PRIMARY_LCIA_METHOD}")
        return PRIMARY_LCIA_METHOD
    raise RuntimeError(f"Primary LCIA method not found: {PRIMARY_LCIA_METHOD}")


# --------------------
# 6. Plotting (unchanged)
# --------------------

def plot_process_co2(process_co2_dict,
                     title="VBR process-wise climate impacts",
                     outfile="vbr_process_co2.png"):
    processes, values = zip(*sorted(
        process_co2_dict.items(),
        key=lambda x: x[1],
        reverse=True
    ))

    plt.figure(figsize=(10, 6))
    bars = plt.bar(range(len(processes)), values)
    plt.xticks(range(len(processes)), processes, rotation=45, ha="right")
    plt.ylabel("kg CO$_2$-eq per kg VBR")
    plt.title(title)

    for bar, val in zip(bars, values):
        height = bar.get_height()
        plt.text(
            bar.get_x() + bar.get_width() / 2,
            height,
            f"{val:.3f}",
            ha="center",
            va="bottom",
            fontsize=8
        )

    plt.tight_layout()
    plt.savefig(outfile, dpi=300)
    plt.close()
    print(f"[INFO] Saved bar chart to {outfile}")


def plot_rep_vs_myers_all(process_co2_dict,
                           rep_furnace_direct,
                           title="The reproduced model vs Myers scenarios (2020 baseline)",
                           outfile="vbr_vs_myers_all_scenarios.png"):

    import numpy as np

    energy_label = "Furnace energy (elec/heat)"
    furnace_total = process_co2_dict.get("Furnace", 0.0)
    total = sum(process_co2_dict.values())
    furnace_energy = max(furnace_total - rep_furnace_direct, 0.0)
    others = total - furnace_energy - rep_furnace_direct

    scenarios_collapsed = {}
    scenarios_collapsed["The reproduced model"] = {
        "Others": others,
        energy_label: furnace_energy,
        "Furnace direct": rep_furnace_direct,
    }

    for name, comp in MYERS_SCENARIOS_2020.items():
        if name == "The reproduced model":
            continue
        furnace_energy_myers = (
            comp.get("Furnace electricity", 0.0)
            + comp.get("Furnace heat (natural gas)", 0.0)
            + comp.get("Furnace heat (other fuels)", 0.0)
        )
        scenarios_collapsed[name] = {
            "Others": comp.get("Other upstream", 0.0),
            energy_label: furnace_energy_myers,
            "Furnace direct": comp.get("Furnace direct", 0.0),
        }

    scenario_names = ["The reproduced model"] + [
        k for k in MYERS_SCENARIOS_2020.keys() if k != "The reproduced model"
    ]
    category_order = ["Others", energy_label, "Furnace direct"]

    n_scen = len(scenario_names)
    y = np.arange(n_scen)

    fig, ax = plt.subplots(figsize=(12, 7))
    left = np.zeros(n_scen)
    for cat in category_order:
        widths = np.array(
            [scenarios_collapsed[scen].get(cat, 0.0) for scen in scenario_names]
        )
        ax.barh(y, widths, left=left, label=cat)
        left += widths

    xlabels = []
    for scen in scenario_names:
        if scen == "The reproduced model":
            xlabels.append("The reproduced model")
        else:
            if scen.startswith("M3."):
                parts = scen.split(" ", 1)
                short = parts[0].replace("M", "")
                if len(parts) > 1:
                    full = parts[1].strip("()")
                    xlabels.append(f"{short}: {full}")
                else:
                    xlabels.append(short)
            else:
                xlabels.append(scen)

    ax.set_yticks(y)
    ax.set_yticklabels(xlabels)
    ax.set_xlabel("kg CO$_2$-eq per kg VBR")
    ax.set_title(title)

    ax.legend(
        loc="upper left",
        bbox_to_anchor=(1.02, 1.0),
        borderaxespad=0.0,
        fontsize=8,
    )

    fig.tight_layout()
    fig.savefig(outfile, dpi=300)
    plt.close(fig)
    print(f"[INFO] Saved HORIZONTAL multi-scenario comparison chart to {outfile}")


def debug_co2_coverage(activity_key, method):
    lca = LCA({activity_key: 1.0}, method)
    lca.lci()
    lca.lcia()

    # contributions by biosphere flow
    inv = lca.characterized_inventory
    bio = get_db("biosphere3")

    co2_keys = [f.key for f in bio if "carbon dioxide" in f['name'].lower() and f['unit']=="kilogram"]

    co2_contrib = 0.0
    total = float(lca.score)

    # inv is sparse matrix; easiest: sum columns for CO2 keys
    for k in co2_keys:
        try:
            idx = lca.biosphere_dict[k]
            co2_contrib += inv[idx].sum()
        except KeyError:
            pass

    print(f"[DEBUG] Method={method}")
    print(f"  total score = {total:.6f}")
    print(f"  CO2 contrib = {co2_contrib:.6f}")
    print(f"  CO2 share   = {co2_contrib/total if total else 0:.3f}")


# --------------------
# 8. Main
# --------------------

def main():
    set_project()
    _ = get_db(BG_DB_NAME)
    _ = get_db("biosphere3")

    params = get_default_params()

    print(f"[INFO] Q_min (kWh per kg VBR): {params['Q_min_kWh']:.4f}")
    print(f"[INFO] Furnace electricity (kWh per kg VBR): {params['E_furn']:.4f}")
    print(f"[INFO] Approx. VBR solid mass per batch: {params['m_vbr']:.4f} kg")

    build_foreground_db(params)

    LCIA_METHOD = choose_lcia_method()
    ensure_co2_only_method()

    process_keys = {
        "Mixing": (FG_DB_NAME, "mixing"),
        "Pelletising": (FG_DB_NAME, "pelletising"),
        "Furnace": (FG_DB_NAME, "furnace"),
        "Cooling": (FG_DB_NAME, "cooling"),
        "Milling": (FG_DB_NAME, "milling"),
        "Transport (raw materials)": (FG_DB_NAME, "transport_raw_materials"),
    }

    # process-wise scores (unchanged)
    process_co2 = {}
    for label, key in process_keys.items():
        lca = LCA({key: 1.0}, LCIA_METHOD)
        lca.lci()
        lca.lcia()
        process_co2[label] = float(lca.score)
        print(f"[RESULT] {label:25s}: {lca.score:.6f} kg CO2-eq / kg VBR")

    total = sum(process_co2.values())
    print(f"[RESULT] TOTAL: {total:.6f} kg CO2-eq / kg VBR")

    # furnace direct (unchanged)
    m_vbr = params["m_vbr"]
    co2_per_kg_coke = 44.0 / 12.0
    co2_per_kg_caco3 = 44.0 / 100.0
    furnace_direct_batch = (
        params["m_coke"] * co2_per_kg_coke
        + params["m_CaCO3"] * co2_per_kg_caco3
    )
    furnace_direct = furnace_direct_batch / m_vbr
    print(f"[RESULT] Furnace direct emissions (stoichiometric, foreground only): "
          f"{furnace_direct:.6f} kg CO2 / kg VBR")

    plot_process_co2(process_co2)
    plot_rep_vs_myers_all(process_co2, furnace_direct)


    # ------------------------------------------
    # Run all 18 ReCiPe 2016 Midpoint (H) + compare to Myers 2020 baseline (S19)
    # ------------------------------------------
    import pandas as pd

    def get_all_recipe_midpoint_methods():
        """Return all ReCiPe 2016 Midpoint (H) methods from ELCD_LCIA."""
        out = []
        for m in methods:
            if (
                m[0] == "ELCD_LCIA"
                and "ReCiPe 2016 Midpoint (H)" in m[1]
            ):
                out.append(m)
        return sorted(out)

    recipe_methods = get_all_recipe_midpoint_methods()
    print("\n[INFO] Found ReCiPe Midpoint(H) method count:", len(recipe_methods))

    # ----  18-category results
    rep_results_18 = {}
    for m in recipe_methods:
        cat = m[2]
        total_cat = 0.0
        for label, key in process_keys.items():
            lca = LCA({key: 1.0}, m)
            lca.lci()
            lca.lcia()
            total_cat += float(lca.score)
        rep_results_18[cat] = total_cat

    # ---- Myers 2020 baseline (industrial, by upscaling approach)
    # Source: Georgiades et al. 2025 SI Table S19.  (units exactly as in S19)
    MYERS_S19_2020 = {
        "M3.1.1": {  # process modelling
            "Terrestrial acidification": 8.20e-4,
            "Global warming": 0.472,
            "Freshwater ecotoxicity": 0.009,
            "Marine ecotoxicity": 0.012,
            "Terrestrial ecotoxicity": 0.403,
            "Fossil resource scarcity": 0.076,
            "Freshwater eutrophication": 2.00e-4,
            "Marine eutrophication": 1.61e-5,
            "Human toxicity, carcinogenic": 0.015,
            "Human toxicity, non-carcinogenic": 0.273,
            "Ionising radiation": 0.129,
            "Agricultural land occupation": 0.007,
            "Mineral resource scarcity": 0.002,
            "Ozone depletion": 1.06e-7,
            "Particulate matter formation": 3.30e-4,
            "Photochemical ozone formation, human health": 4.60e-4,
            "Photochemical ozone formation, terrestrial ecosystems": 4.80e-4,
            "Water consumption": 0.003,
        },
        "M3.2.1": {  # linear regression
            "Terrestrial acidification": 7.50e-4,
            "Global warming": 0.448,
            "Freshwater ecotoxicity": 0.008,
            "Marine ecotoxicity": 0.011,
            "Terrestrial ecotoxicity": 0.370,
            "Fossil resource scarcity": 0.070,
            "Freshwater eutrophication": 1.80e-4,
            "Marine eutrophication": 1.46e-5,
            "Human toxicity, carcinogenic": 0.014,
            "Human toxicity, non-carcinogenic": 0.247,
            "Ionising radiation": 0.117,
            "Agricultural land occupation": 0.007,
            "Mineral resource scarcity": 0.002,
            "Ozone depletion": 9.61e-8,
            "Particulate matter formation": 3.00e-4,
            "Photochemical ozone formation, human health": 4.20e-4,
            "Photochemical ozone formation, terrestrial ecosystems": 4.40e-4,
            "Water consumption": 0.002,
        },
        "M3.2.2": {  # power regression
            "Terrestrial acidification": 7.80e-4,
            "Global warming": 0.458,
            "Freshwater ecotoxicity": 0.009,
            "Marine ecotoxicity": 0.012,
            "Terrestrial ecotoxicity": 0.383,
            "Fossil resource scarcity": 0.073,
            "Freshwater eutrophication": 1.90e-4,
            "Marine eutrophication": 1.52e-5,
            "Human toxicity, carcinogenic": 0.014,
            "Human toxicity, non-carcinogenic": 0.258,
            "Ionising radiation": 0.122,
            "Agricultural land occupation": 0.007,
            "Mineral resource scarcity": 0.002,
            "Ozone depletion": 1.00e-7,
            "Particulate matter formation": 3.10e-4,
            "Photochemical ozone formation, human health": 4.40e-4,
            "Photochemical ozone formation, terrestrial ecosystems": 4.60e-4,
            "Water consumption": 0.003,
        },
        "M3.3.1": {  # clinker proxy – Centre et al.
            "Terrestrial acidification": 1.38e-3,
            "Global warming": 0.480,
            "Freshwater ecotoxicity": 0.005,
            "Marine ecotoxicity": 0.006,
            "Terrestrial ecotoxicity": 0.620,
            "Fossil resource scarcity": 0.073,
            "Freshwater eutrophication": 1.00e-3,
            "Marine eutrophication": 8.97e-6,
            "Human toxicity, carcinogenic": 0.009,
            "Human toxicity, non-carcinogenic": 0.172,
            "Ionising radiation": 0.015,
            "Agricultural land occupation": 0.003,
            "Mineral resource scarcity": 0.001,
            "Ozone depletion": 5.00e-8,
            "Particulate matter formation": 5.10e-4,
            "Photochemical ozone formation, human health": 6.00e-4,
            "Photochemical ozone formation, terrestrial ecosystems": 6.40e-4,
            "Water consumption": 0.001,
        },
        "M3.3.2": {  # clinker proxy – ecoinvent
            "Terrestrial acidification": 2.79e-3,
            "Global warming": 0.709,
            "Freshwater ecotoxicity": 0.010,
            "Marine ecotoxicity": 0.014,
            "Terrestrial ecotoxicity": 0.442,
            "Fossil resource scarcity": 0.110,
            "Freshwater eutrophication": 2.70e-3,
            "Marine eutrophication": 1.76e-5,
            "Human toxicity, carcinogenic": 0.021,
            "Human toxicity, non-carcinogenic": 0.427,
            "Ionising radiation": 0.027,
            "Agricultural land occupation": 0.005,
            "Mineral resource scarcity": 0.001,
            "Ozone depletion": 9.14e-8,
            "Particulate matter formation": 9.80e-4,
            "Photochemical ozone formation, human health": 1.00e-3,
            "Photochemical ozone formation, terrestrial ecosystems": 1.00e-3,
            "Water consumption": 8.60e-4,
        },
        "M3.3.3": {  # ceramic proxy – BAT
            "Terrestrial acidification": 3.30e-3,
            "Global warming": 0.563,
            "Freshwater ecotoxicity": 0.003,
            "Marine ecotoxicity": 0.004,
            "Terrestrial ecotoxicity": 0.286,
            "Fossil resource scarcity": 0.125,
            "Freshwater eutrophication": 5.00e-4,
            "Marine eutrophication": 3.97e-6,
            "Human toxicity, carcinogenic": 0.008,
            "Human toxicity, non-carcinogenic": 0.084,
            "Ionising radiation": 0.016,
            "Agricultural land occupation": 0.002,
            "Mineral resource scarcity": 0.002,
            "Ozone depletion": 4.44e-8,
            "Particulate matter formation": 1.70e-4,
            "Photochemical ozone formation, human health": 3.90e-4,
            "Photochemical ozone formation, terrestrial ecosystems": 4.40e-4,
            "Water consumption": 6.60e-4,
        },
        "M3.3.4": {  # ceramic proxy – ecoinvent
            "Terrestrial acidification": 6.60e-4,
            "Global warming": 0.785,
            "Freshwater ecotoxicity": 0.007,
            "Marine ecotoxicity": 0.010,
            "Terrestrial ecotoxicity": 0.393,
            "Fossil resource scarcity": 0.192,
            "Freshwater eutrophication": 1.30e-4,
            "Marine eutrophication": 1.01e-5,
            "Human toxicity, carcinogenic": 0.015,
            "Human toxicity, non-carcinogenic": 0.195,
            "Ionising radiation": 0.069,
            "Agricultural land occupation": 0.005,
            "Mineral resource scarcity": 0.004,
            "Ozone depletion": 9.53e-8,
            "Particulate matter formation": 3.00e-4,
            "Photochemical ozone formation, human health": 6.40e-4,
            "Photochemical ozone formation, terrestrial ecosystems": 7.10e-4,
            "Water consumption": 0.002,
        },
    }

    # ---- build comparison table
    df = pd.DataFrame({"The reproduced model": rep_results_18})


    # ----------------------------------------------------------
    # 9. Fix category-name mismatch between ELCD ReCiPe and Myers
    # ----------------------------------------------------------

    name_map = {
        # ELCD ReCiPe name ------------------------------- Myers Table S19 name
        "Fine particulate matter formation": "Particulate matter formation",
        "Human carcinogenic toxicity": "Human toxicity, carcinogenic",
        "Human non-carcinogenic toxicity": "Human toxicity, non-carcinogenic",
        "Ozone formation, Human health": "Photochemical ozone formation, human health",
        "Ozone formation, Terrestrial ecosystems": "Photochemical ozone formation, terrestrial ecosystems",
        "Stratospheric ozone depletion": "Ozone depletion",
        "Land use": "Agricultural land occupation",
        "Ionizing radiation": "Ionising radiation",
    }

    # rename model categories using the above map
    df.rename(index=name_map, inplace=True)

    # Map ELCD ReCiPe category names to Myers short labels if needed.
    # If names already match, this just works.
    for scen, vals in MYERS_S19_2020.items():
        df[scen] = pd.Series(vals)

    # quick delta vs M3.1.1 (process modelling baseline)
    if "Global warming" in df.index:
        df["Δ vs M3.1.1"] = df["The reproduced model"] - df["M3.1.1"]

    df = df.sort_index()

    pd.set_option("display.float_format", lambda x: f"{x:.2e}")

    print("\n[COMPARISON TABLE]")
    print(df)

    df.to_excel("vbr_recipe18_vs_myers2020.xlsx")
    print("[INFO] Saved comparison to vbr_recipe18_vs_myers2020.xlsx")


    debug_co2_coverage((FG_DB_NAME, "furnace"), PRIMARY_LCIA_METHOD)
    debug_co2_coverage((FG_DB_NAME, "furnace"), CO2_ONLY_METHOD)


    # 10) Top LCIA contributors
    from bw2data import get_activity

    def print_top_contributors(fu_key, method, top_n=30, title=None):
        """
        Print top characterized biosphere contributors for a given FU activity.
        Robust across BW versions and key shapes.
        Shows extra identifiers: name, categories, unit.
        """
        lca_tc = LCA({fu_key: 1.0}, method)
        lca_tc.lci()
        lca_tc.lcia()

        biosphere_map = getattr(lca_tc, "biosphere_dict", None)
        if biosphere_map is None:
            biosphere_map = lca_tc.dicts.biosphere  # older BW

        cf = lca_tc.characterized_inventory

        data = []
        inv = {v: k for k, v in biosphere_map.items()}

        try:
            cf_coo = cf.tocoo()
            for idx, val in zip(cf_coo.row, cf_coo.data):
                if val == 0:
                    continue
                fk = inv.get(idx)
                if fk is None:
                    continue
                # normalize weird key shapes
                if isinstance(fk, tuple) and len(fk) == 1 and isinstance(fk[0], tuple):
                    fk = fk[0]
                data.append((fk, float(val)))
        except Exception:
            for fk, idx in biosphere_map.items():
                try:
                    val = float(cf[idx])
                except Exception:
                    continue
                if val == 0:
                    continue
                if isinstance(fk, tuple) and len(fk) == 1 and isinstance(fk[0], tuple):
                    fk = fk[0]
                data.append((fk, val))

        data.sort(key=lambda x: abs(x[1]), reverse=True)

        print("\n========== TOP {} CONTRIBUTORS {} ==========".format(
            top_n, f"({title})" if title else ""
        ))

        for fk, amt in data[:top_n]:
            flow = None
            try:
                flow = get_activity(fk)  # works for tuple key or int id
            except Exception:
                flow = None

            if flow:
                name = flow.get("name", "NA")
                cats = flow.get("categories", None)
                cats_str = "/".join(cats) if isinstance(cats, (list, tuple)) else str(cats)
                unit = flow.get("unit", "NA")
                db = fk[0] if isinstance(fk, tuple) and len(fk) > 0 else "NA"
                code = fk[1] if isinstance(fk, tuple) and len(fk) > 1 else str(fk)
                print(f"  - {name} | cats={cats_str} | unit={unit} | {code} | {db} | {amt:.6f}")
            else:
                # fallback if flow lookup fails
                if isinstance(fk, tuple) and len(fk) == 2:
                    print(f"  - {fk[1]} | {fk[0]} | {amt:.6f}")
                else:
                    print(f"  - {fk} | {amt:.6f}")

        print("================================================\n")

    print_top_contributors((FG_DB_NAME, "furnace"), PRIMARY_LCIA_METHOD, top_n=10,
                           title="PRIMARY_LCIA_METHOD")
    print_top_contributors((FG_DB_NAME, "furnace"), CO2_ONLY_METHOD, top_n=10,
                           title="CO2_ONLY_METHOD")


if __name__ == "__main__":
    main()
