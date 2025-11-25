# 22_import_elcd_lcia_via_ipc_ver3_uuid_preserve.py
# ------------------------------------------------------------
# openLCA IPC → Brightway background DB + LCIA importer (ver3 UUID-PRESERVE)
#
# Core idea:
#  1) Import ELCD processes as before.
#  2) EXTEND biosphere3 by adding ALL ELCD elementary flows:
#       - bw flow code = ELCD flow UUID (string)
#       - name / categories / CAS / formula copied if available
#  3) When importing processes, biosphere exchanges point to these new flows directly.
#  4) When importing LCIA, use LCIA factor UUID → biosphere3(code=same UUID).
#  5) Handle CF duplication by averaging identical biosphere keys within a method.
#
# Result: true 1:1:1 chain, no mapping collapse, no CF inflation.
# ------------------------------------------------------------

import sys, traceback, json, os, re
from collections import defaultdict

import olca_ipc as olca
import olca_schema as o
from bw2data import projects, Database, databases, Method, methods


# ================= USER SETTINGS =================
IPC_PORT = 8080

# Brightway
PROJECT_NAME = "vbr_project"
BW_BG_DB_NAME = "elcd_bw_ipc"
BW_BIOSPHERE_DB = "biosphere3"

# openLCA DB that contains ELCD + LCIA
SOURCE_OLCA_DB = "ELCD_LCIA"

# Optional filters
SKIP_DUMMY = True
REQUIRE_PRODUCTION = True
MAX_ACTIVITIES = None

# LCIA import switches
IMPORT_LCIA = True
LCIA_METHOD_NAME_FILTER = None
DELETE_EXISTING_LCIA_PREFIX = ("VBR",)

# JSON outputs
ELCD_FLOW_MAP_JSON = "elcd_uuid_to_biosphere3_key.json"
LCIA_FLOW_MAP_JSON = "lcia_uuid_to_biosphere3_key.json"
# ================================================


# ==================================================
# IPC helpers
# ==================================================
def connect_olca():
    return olca.Client(IPC_PORT)

def set_bw_project():
    projects.set_current(PROJECT_NAME)

def set_olca_db(client):
    if hasattr(client, "set_database"):
        client.set_database(SOURCE_OLCA_DB)

def get_all_process_descs(client):
    return client.get_descriptors(o.Process)

def load_full_process(client, desc):
    return client.get(o.Process, desc.id)

def get_all_flow_descs(client):
    return client.get_descriptors(o.Flow)

def load_full_flow(client, flow_uuid):
    return client.get(o.Flow, flow_uuid)

def get_all_impact_method_descs(client):
    return client.get_descriptors(o.ImpactMethod)

def load_full_impact_method(client, desc):
    return client.get(o.ImpactMethod, desc.id)

def get_full_impact_category(client, cat_stub):
    cid = getattr(cat_stub, "id", None)
    if not cid:
        return cat_stub
    try:
        return client.get(o.ImpactCategory, cid)
    except Exception:
        return cat_stub


# ==================================================
# JSON IO
# ==================================================
def save_map_json(flow_map, path):
    out = {k: [v[0], v[1]] for k, v in flow_map.items()}
    with open(path, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2, ensure_ascii=False)
    print(f"[DONE] Saved map: {len(out)} entries → {path}")

def load_map_json(path):
    if not os.path.exists(path):
        print(f"[WARN] No mapping JSON found: {path}")
        return {}
    with open(path, "r", encoding="utf-8") as f:
        raw = json.load(f)
    rebuilt = {k: (v[0], v[1]) for k, v in raw.items()}
    print(f"[INFO] Loaded map: {len(rebuilt)} entries ← {path}")
    return rebuilt


# ==================================================
# Category helpers (best-effort copy)
# ==================================================
def get_full_category_path(cat_obj):
    """Return tuple like ('emissions to air', 'unspecified') if possible."""
    if not cat_obj:
        return ()
    if isinstance(cat_obj, str):
        parts = [p.strip() for p in re.split(r"[\\/>\|]+", cat_obj) if p.strip()]
        return tuple(parts)
    parts = []
    cur = cat_obj
    safety = 0
    while cur and safety < 20:
        nm = getattr(cur, "name", None)
        if nm:
            parts.append(str(nm).strip())
        cur = getattr(cur, "category", None) or getattr(cur, "parent", None)
        safety += 1
    return tuple(reversed(parts))


# ==================================================
# 1) EXTEND biosphere3 with ELCD elementary flows
# ==================================================
def extend_biosphere3_with_elcd_flows(client):
    bio = Database(BW_BIOSPHERE_DB)

    flow_descs = get_all_flow_descs(client)
    print(f"[INFO] Flow descriptors found: {len(flow_descs)}")

    # Existing biosphere3 codes for fast lookup
    existing_codes = set()
    for f in bio:
        existing_codes.add(f.key[1])

    elcd_flow_map = {}
    created = 0
    skipped = 0
    failed = 0

    for fd in flow_descs:
        fl = load_full_flow(client, fd.id)
        if not fl or not fl.name:
            continue
        if fl.flow_type != o.FlowType.ELEMENTARY_FLOW:
            continue

        code = str(fl.id)
        bw_key = (BW_BIOSPHERE_DB, code)
        elcd_flow_map[code] = bw_key  # UUID -> biosphere3 key (same UUID)

        if code in existing_codes:
            skipped += 1
            continue

        cat_path = get_full_category_path(getattr(fl, "category", None))
        cats = tuple(cat_path) if cat_path else ()

        data = {
            "name": (fl.name or "").strip(),
            "unit": fl.reference_unit.name if getattr(fl, "reference_unit", None) else None,
            "categories": cats,
            "type": "emission",  # safe default; BW doesn't require strict enum here
        }

        # copy optional metadata if present
        cas = getattr(fl, "cas_number", None) or getattr(fl, "casNumber", None)
        if cas:
            data["CAS number"] = str(cas).strip()

        formula = getattr(fl, "formula", None)
        if formula:
            data["formula"] = str(formula).strip()

        syns = getattr(fl, "synonyms", None)
        if syns:
            data["synonyms"] = syns if isinstance(syns, list) else [str(syns)]

        try:
            bio.new_activity(code=code, **data).save()
            created += 1
        except Exception:
            failed += 1
            print(f"[DEBUG] Failed to create biosphere3 flow: {fl.name} | UUID={code}")

    print(f"[INFO] ELCD biosphere extension created: {created}")
    print(f"[INFO] ELCD biosphere extension skipped (already existed): {skipped}")
    print(f"[INFO] ELCD biosphere extension failures: {failed}")

    save_map_json(elcd_flow_map, ELCD_FLOW_MAP_JSON)
    return elcd_flow_map


# ==================================================
# Exchange conversion (now UUID-direct)
# ==================================================
def convert_exchange(exc, proc_key_map, uuid_to_bio_key):
    # TECHNOSPHERE
    if exc.flow and exc.flow.flow_type == o.FlowType.PRODUCT_FLOW:
        provider_id = getattr(exc, "default_provider", None)
        if provider_id and provider_id.id in proc_key_map:
            return {
                "input": proc_key_map[provider_id.id],
                "amount": float(exc.amount or 0),
                "type": "technosphere",
                "unit": exc.unit.name if exc.unit else None,
            }
        return None

    # BIOSPHERE: use UUID-direct biosphere3 key
    if exc.flow and exc.flow.flow_type in (o.FlowType.ELEMENTARY_FLOW, o.FlowType.WASTE_FLOW):
        fid = str(exc.flow.id)
        if fid not in uuid_to_bio_key:
            return None
        return {
            "input": uuid_to_bio_key[fid],
            "amount": float(exc.amount or 0),
            "type": "biosphere",
            "unit": exc.unit.name if exc.unit else None,
        }

    return None


# ==================================================
# LCIA helpers
# ==================================================
def maybe_delete_existing_lcia(prefix):
    if not prefix:
        return
    to_delete = [m for m in methods if tuple(m)[:len(prefix)] == tuple(prefix)]
    for m in to_delete:
        try:
            del methods[m]
        except:
            pass
    if to_delete:
        print(f"[WARN] Deleted old methods with prefix {prefix}: {len(to_delete)}")

def collect_lcia_factor_flow_ids(client):
    method_descs = get_all_impact_method_descs(client)
    lcia_flow_ids = set()

    for md in method_descs:
        im = load_full_impact_method(client, md)
        if not im:
            continue

        for cat_stub in getattr(im, "impact_categories", None) or []:
            cat = get_full_impact_category(client, cat_stub)

            for fac in getattr(cat, "impact_factors", None) or []:
                fl = getattr(fac, "flow", None)
                if fl and getattr(fl, "id", None):
                    lcia_flow_ids.add(str(fl.id))

    return lcia_flow_ids


def average_duplicate_cfs(cf_list):
    """
    cf_list: list of (bw_key, value)
    Returns list where duplicated keys are averaged.
    """
    bucket = defaultdict(list)
    for k, v in cf_list:
        bucket[k].append(v)

    out = []
    for k, vs in bucket.items():
        if len(vs) == 1:
            out.append((k, vs[0]))
        else:
            out.append((k, sum(vs) / len(vs)))
    return out


def import_lcia_methods_uuid_direct(client, uuid_to_bio_key):
    method_descs = get_all_impact_method_descs(client)
    print(f"[INFO] LCIA descriptors: {len(method_descs)}")

    imported = 0
    skipped = 0

    for md in method_descs:
        im = load_full_impact_method(client, md)
        if not im or not im.name:
            continue

        if LCIA_METHOD_NAME_FILTER and LCIA_METHOD_NAME_FILTER.lower() not in im.name.lower():
            continue

        for cat_stub in getattr(im, "impact_categories", None) or []:
            cat = get_full_impact_category(client, cat_stub)
            cat_name = (getattr(cat, "name", None) or "").strip()
            if not cat_name:
                continue

            bw_method_key = (SOURCE_OLCA_DB, im.name, cat_name)
            raw_cf_list = []

            for fac in getattr(cat, "impact_factors", None) or []:
                fl = getattr(fac, "flow", None)
                if not fl:
                    continue
                fid = str(fl.id)
                if fid not in uuid_to_bio_key:
                    continue

                val = getattr(fac, "value", None)
                if val is None:
                    continue

                raw_cf_list.append((uuid_to_bio_key[fid], float(val)))

            if not raw_cf_list:
                skipped += 1
                continue

            # IMPORTANT: average duplicated keys to avoid inflation
            cf_list = average_duplicate_cfs(raw_cf_list)

            m = Method(bw_method_key)
            m.register()
            m.write(cf_list)
            imported += 1

    print(f"[DONE] LCIA import: {imported} imported, {skipped} skipped")


# ==================================================
# MAIN
# ==================================================
def main():
    set_bw_project()
    print(f"[INFO] BW project       : {PROJECT_NAME}")
    print(f"[INFO] BW biosphere DB : {BW_BIOSPHERE_DB}")
    print(f"[INFO] BW BG DB name   : {BW_BG_DB_NAME}")

    client = connect_olca()
    set_olca_db(client)
    print(f"[INFO] IPC port        : {IPC_PORT}")
    print(f"[INFO] openLCA DB      : {SOURCE_OLCA_DB}\n")

    # ----------------------------------------------------------
    # A) Extend biosphere3 with ELCD elementary flows (UUID-preserve)
    # ----------------------------------------------------------
    uuid_to_bio_key = extend_biosphere3_with_elcd_flows(client)

    # ----------------------------------------------------------
    # B) Import ELCD processes → BW background DB
    # ----------------------------------------------------------
    proc_descs = get_all_process_descs(client)
    if MAX_ACTIVITIES:
        proc_descs = proc_descs[:MAX_ACTIVITIES]
    print(f"[INFO] Process descriptors found: {len(proc_descs)}")

    if BW_BG_DB_NAME in databases:
        print(f"[WARN] BW DB '{BW_BG_DB_NAME}' exists → deleting and recreating.")
        del databases[BW_BG_DB_NAME]
    bw_db = Database(BW_BG_DB_NAME)

    proc_key_map = {}
    bw_data = {}

    # shells
    for d in proc_descs:
        p = load_full_process(client, d)
        if not p:
            continue

        name = (p.name or "").strip()
        if SKIP_DUMMY and "dummy" in name.lower():
            continue

        bw_key = (BW_BG_DB_NAME, p.id)
        proc_key_map[p.id] = bw_key

        ref_unit = "kilogram"
        ref_prod = None
        qr = getattr(p, "quantitative_reference", None)

        if qr:
            if getattr(qr, "unit", None):
                ref_unit = qr.unit.name
            if getattr(qr, "flow", None):
                ref_prod = qr.flow.name

        if not ref_prod:
            for exc in p.exchanges or []:
                if exc.flow and not exc.is_input:
                    ref_prod = exc.flow.name
                    if exc.unit:
                        ref_unit = exc.unit.name
                    break

        bw_data[bw_key] = {
            "name": name,
            "unit": ref_unit,
            "location": p.location.name if p.location else None,
            "reference product": ref_prod,
            "exchanges": []
        }

    print(f"[INFO] BW shells created: {len(bw_data)}")

    # ----------------------------------------------------------
    # C) Fill exchanges and write BW BG DB (biosphere = UUID-direct)
    # ----------------------------------------------------------
    kept, skipped_prod = 0, 0

    for d in proc_descs:
        p = load_full_process(client, d)
        if not p or p.id not in proc_key_map:
            continue

        bw_key = proc_key_map[p.id]
        ex_list = []

        qr = getattr(p, "quantitative_reference", None)
        prod_added = False

        if qr and getattr(qr, "flow", None):
            # A-style: preserve openLCA quantitative reference amount
            prod_amt = float(getattr(qr, "amount", None) or 1.0)
            if prod_amt == 0:
                prod_amt = 1.0

            ex_list.append({
                "input": bw_key,
                "amount": prod_amt,
                "type": "production",
                "unit": qr.unit.name if qr.unit else None,
            })
            prod_added = True
    
        if not prod_added:
            for exc in p.exchanges or []:
                if exc.flow and not exc.is_input:
                    prod_amt = float(getattr(exc, "amount", None) or 1.0)
                    if prod_amt == 0:
                        prod_amt = 1.0

                    ex_list.append({
                        "input": bw_key,
                        "amount": prod_amt,
                        "type": "production",
                        "unit": exc.unit.name if exc.unit else None,
                    })
                    prod_added = True
                    break

        for exc in p.exchanges or []:
            bw_exc = convert_exchange(exc, proc_key_map, uuid_to_bio_key)
            if bw_exc:
                ex_list.append(bw_exc)

        if REQUIRE_PRODUCTION and not any(e["type"] == "production" for e in ex_list):
            skipped_prod += 1
            continue

        bw_data[bw_key]["exchanges"] = ex_list
        kept += 1

    print(f"[INFO] Activities kept: {kept}")
    print(f"[INFO] Skipped (no production): {skipped_prod}")

    bw_db.write(bw_data)
    print(f"[DONE] BW background DB created: '{BW_BG_DB_NAME}'")

    # ----------------------------------------------------------
    # D) LCIA import (UUID-direct + average duplicates)
    # ----------------------------------------------------------
    if IMPORT_LCIA:
        if DELETE_EXISTING_LCIA_PREFIX:
            maybe_delete_existing_lcia(DELETE_EXISTING_LCIA_PREFIX)

        lcia_flow_ids = collect_lcia_factor_flow_ids(client)
        print(f"[INFO] Unique LCIA factor flows: {len(lcia_flow_ids)}")

        # Map is trivial now: LCIA UUID -> biosphere3(UUID)
        lcia_flow_map = {}
        missing = 0
        for fid in lcia_flow_ids:
            if fid in uuid_to_bio_key:
                lcia_flow_map[fid] = uuid_to_bio_key[fid]
            else:
                missing += 1

        print(f"[INFO] LCIA mapped flows: {len(lcia_flow_map)}")
        print(f"[INFO] LCIA missing (no biosphere UUID): {missing}")

        save_map_json(lcia_flow_map, LCIA_FLOW_MAP_JSON)

        import_lcia_methods_uuid_direct(client, lcia_flow_map)


if __name__ == "__main__":
    try:
        main()
    except Exception:
        print("[FATAL] IPC→BW import failed.")
        traceback.print_exc()
        sys.exit(1)
