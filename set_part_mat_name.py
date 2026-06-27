#saiftherehman
"""
Set Part, Material, and Property Name by Free-Mesh Proximity (ANSA v22.1.5)

Three independent operations, each gated by a flag at the top of the file:

  SET_PART            reassign FE elements to the matched CAD ANSAPART
  SET_MATERIALS       write CAD material MID onto LSDYNA properties
                        (passes int ID — ANSA rejects cross-deck MAT1 entity)
                        then SynchronizeMaterials propagates to NASTRAN /
                        ABAQUS / PAM-CRASH with delete_released=True
  SET_PROPERTY_NAMES  copy the ANSAPART's Name onto each property card
                        (NASTRAN + LSDYNA), after running it through
                        _sanitize_name (collapses non-alnum -> underscore,
                        strips duplicated prefix and trailing CAD format
                        hints, merges short revision-code runs)

The cluster→part voting runs in every configuration; the flags only gate
the write steps and their downstream effects.

Deck usage:
  FE_DECK  = NASTRAN — FE elements live here
  GEO_DECK = NASTRAN — CAD geometry (FACE, ANSAPART), temp mesh shells
  MAT_DECK = LSDYNA  — primary deck for material assignment

Run inside ANSA: Scripts → Run Script → select this file.
"""

from ansa import base, constants, mesh
import math
import re
from collections import defaultdict


# ─── Configuration ────────────────────────────────────────────────────────────

# Bare-identifier aliases so any of these work for flag values:
#   T / F      Yes / No      YES / NO      Y / N      True / False
T = Y = Yes = YES = True
F = N = No = NO = False

# What to do for each matched cluster.  Cluster→part voting still runs in all
# cases; these flags only gate the write steps and their downstream effects.
# Each flag accepts the aliases above OR strings:  "T"/"F", "yes"/"no", "1"/"0".
SET_PART            = T   # reassign FE elements to the matched CAD ANSAPART
SET_MATERIALS       = T   # write CAD material MID onto LSDYNA properties
SET_PROPERTY_NAMES  = T   # copy ANSAPART name onto NASTRAN+LSDYNA property cards


def _truthy(v):
    if isinstance(v, str):
        return v.strip().lower() in ("t", "true", "yes", "y", "1", "on")
    return bool(v)


SET_PART           = _truthy(SET_PART)
SET_MATERIALS      = _truthy(SET_MATERIALS)
SET_PROPERTY_NAMES = _truthy(SET_PROPERTY_NAMES)

FE_DECK  = constants.NASTRAN  # deck that owns the FE elements
GEO_DECK = constants.NASTRAN  # deck used to access CAD geometry
MAT_DECK = constants.LSDYNA   # deck for material assignment (LS-DYNA)

# Temporary mesh element size for CAD reference (model units, e.g. mm).
TEMP_MESH_SIZE = 10.0

# FE element types to reassign
FE_ELEMENT_TYPES = ["SHELL", "SOLID"]

# NASTRAN node connectivity fields (G1–G20 covers all shell/solid types)
NODE_FIELDS = ["G1", "G2", "G3", "G4", "G5", "G6", "G7", "G8",
               "G9", "G10", "G11", "G12", "G13", "G14", "G15", "G16",
               "G17", "G18", "G19", "G20"]

# NASTRAN node connectivity fields — used only for geometry-side temp mesh shells
_GEO_NODE_FIELDS = ["G1", "G2", "G3", "G4", "G5", "G6", "G7", "G8"]

# Warn if nearest reference point is farther than this (model units)
DISTANCE_WARNING = 50.0

# True = preview; SetEntityPart and material assignment are NOT called
DRY_RUN = False


# ─── Spatial index ────────────────────────────────────────────────────────────

class SpatialGrid:
    """O(1) average nearest-point lookup via grid hashing."""

    def __init__(self, cell_size):
        self._cs = cell_size
        self._cells = defaultdict(list)

    def _key(self, pos):
        cs = self._cs
        return (int(math.floor(pos[0] / cs)),
                int(math.floor(pos[1] / cs)),
                int(math.floor(pos[2] / cs)))

    def insert(self, pos, cad_part):
        self._cells[self._key(pos)].append((pos, cad_part))

    def nearest(self, qpos):
        ix, iy, iz = self._key(qpos)
        best_part, best_sq = None, float("inf")
        for r in range(30):
            for dx in range(-r, r + 1):
                for dy in range(-r, r + 1):
                    for dz in range(-r, r + 1):
                        if max(abs(dx), abs(dy), abs(dz)) != r:
                            continue
                        for pos, part in self._cells.get((ix+dx, iy+dy, iz+dz), []):
                            sq = ((pos[0]-qpos[0])**2 +
                                  (pos[1]-qpos[1])**2 +
                                  (pos[2]-qpos[2])**2)
                            if sq < best_sq:
                                best_sq, best_part = sq, part
            if best_part is not None and r * self._cs > math.sqrt(best_sq):
                break
        return best_part, math.sqrt(best_sq)


# ─── Helpers ──────────────────────────────────────────────────────────────────

def _node_entities(elem):
    """Node entities of an FE element (LS-DYNA deck, N1–N8)."""
    vals = elem.get_entity_values(FE_DECK, NODE_FIELDS)
    return [v for v in vals.values() if v is not None and hasattr(v, "_id")]


def _geo_node_entities(elem):
    """Node entities of a geometry temp-mesh shell (NASTRAN deck, G1–G8)."""
    vals = elem.get_entity_values(GEO_DECK, _GEO_NODE_FIELDS)
    return [v for v in vals.values() if v is not None and hasattr(v, "_id")]


# ─── Step 1: Snapshot existing FE elements ────────────────────────────────────

def collect_all_fe_elements():
    elems = []
    for etype in FE_ELEMENT_TYPES:
        elems.extend(base.CollectEntities(FE_DECK, None, etype) or [])
    return elems


# ─── Step 2: FE cluster detection via Neighb ─────────────────────────────────

def find_fe_clusters(all_fe_elems):
    """
    Find truly connected FE components using ANSA Neighb.
    Or([elem]) isolates one element; Neighb("ALL") flood-fills to all
    mesh-connected neighbours; CollectEntities(filter_visible=True) harvests it.
    """
    remaining = {e._id: e for e in all_fe_elems}
    clusters = []

    while remaining:
        _, elem = next(iter(remaining.items()))

        base.Or(entities=[elem])
        base.Neighb("ALL")

        cluster = []
        for etype in FE_ELEMENT_TYPES:
            cluster.extend(
                base.CollectEntities(FE_DECK, None, etype,
                                     filter_visible=True) or []
            )

        for c in cluster:
            remaining.pop(c._id, None)

        clusters.append(cluster)

    return clusters


# ─── Step 3: Temp-mesh reference clouds per CAD part ─────────────────────────

def collect_cad_parts():
    """
    CAD parts = ANSAPART entities that own FACE geometry.
    Returns: {part_entity: {'name': str, 'faces': [entity, ...]}}
    """
    all_parts = base.CollectEntities(GEO_DECK, None, "ANSAPART") or []
    cad_parts = {}
    for part in all_parts:
        faces = base.CollectEntities(GEO_DECK, part, "FACE",
                                     recursive=True) or []
        if not faces:
            continue
        vals = part.get_entity_values(GEO_DECK, ["Name"])
        name = vals.get("Name") or f"Part_{part._id}"
        cad_parts[part] = {"name": name, "faces": faces}
    return cad_parts


def build_reference_clouds(cad_parts, existing_fe_ids):
    """
    For each CAD part:
      base.Or(faces) → mesh.CreateFreeMesh() → collect node positions
      → mesh.EraseMesh()
      → capture material via CollectEntities(LSDYNA, part, mat_from_entities=True)
        (reads from existing FE already in the part, not from temp shells)
    Temp shells are in GEO_DECK (NASTRAN); node fields are G1–G8.
    Returns: (clouds, cad_materials)
      clouds        = {part_entity: [(x,y,z), ...]}
      cad_materials = {part_entity: material_entity_or_None}
    """
    mesh.SetMeshParamTargetLength("absolute", TEMP_MESH_SIZE)
    clouds = {}
    cad_materials = {}

    for part, data in cad_parts.items():
        name  = data["name"]
        faces = data["faces"]

        base.Or(entities=faces)
        mesh.CreateFreeMesh()

        positions  = []
        seen_nodes = set()

        for face in faces:
            temp_shells = base.CollectEntities(GEO_DECK, face, "SHELL") or []
            for shell in temp_shells:
                if shell._id in existing_fe_ids:
                    continue
                for node in _geo_node_entities(shell):
                    if node._id not in seen_nodes:
                        seen_nodes.add(node._id)
                        pos = node.position
                        if pos is not None:
                            positions.append((float(pos[0]),
                                              float(pos[1]),
                                              float(pos[2])))

        mesh.EraseMesh()

        # Capture the CAD part's material from its EXISTING FE elements (LSDYNA
        # deck). This must happen before SetEntityPart shuffles elements around.
        mats = base.CollectEntities(MAT_DECK, part, "__MATERIALS__",
                                     mat_from_entities=True) or []
        mat_found = mats[0] if mats else None
        cad_materials[part] = mat_found

        if positions:
            clouds[part] = positions
            mat_info = f"  mat_id={mat_found._id}" if mat_found else "  [no mat]"
            print(f"    '{name}':  {len(faces)} faces "
                  f"→ {len(positions)} ref pts  [free mesh]{mat_info}")
        else:
            print(f"  [WARN] '{name}': no reference points — skipped.")

    return clouds, cad_materials


def build_spatial_index(clouds):
    grid = SpatialGrid(TEMP_MESH_SIZE)
    total = 0
    for part, positions in clouds.items():
        for pos in positions:
            grid.insert(pos, part)
        total += len(positions)
    print(f"  Spatial index: {total:,} reference points across {len(clouds)} parts.")
    return grid


# ─── Name sanitization ────────────────────────────────────────────────────────
#
# Verified against the examples in naming_example.txt:
#   A1925406408_9_41_0009.0041_A1925406408_9_41_0009.0041_HALTER BATTERIE
#     -> A1925406408_9_41_0009_0041_HALTER_BATTERIE
#   9M1.803.443--01, FEDERBEINAUFNAHME [prd-p1-01971834]
#     -> 9M1_803_443_01_FEDERBEINAUFNAHME_prd_p1_01971834
#   A5304012000_DICASTAL_03_VA_7X23_PRT_20260409_stp.igs
#     -> A5304012000_DICASTAL_03_VA_7X23_PRT_20260409
#   1FA_857_003____DMU_TM__007_____Z_INSTRUMENTENTAFE______VZM7639_VZM7578
#     -> 1FA_857_003_DMU_TM_007_Z_INSTRUMENTENTAFE_VZM7639_VZM7578
#   P5FHFL5 OPTIMIERUNGSERGEBNIS  A 14 B NONFRG 5P  0
#     -> P5FHFL5_A14B_OPTIMIERUNGSERGEBNIS_NONFRG_5P_0

_FORMAT_HINTS = {"stp", "step", "igs", "iges", "prt", "sldprt", "ipt",
                  "stl", "obj", "dxf", "x_t", "x_b"}
_REV_PATTERN = re.compile(r'^[A-Za-z]\d+[A-Za-z]$')


def _merge_short_tokens(tokens):
    """Merge single-alpha + short-tokens runs into one (A_14_B -> A14B)."""
    def _short(t):
        return len(t) == 1 or (t.isdigit() and len(t) <= 2)
    out, i = [], 0
    while i < len(tokens):
        t = tokens[i]
        if len(t) == 1 and t.isalpha():
            group, j = [t], i + 1
            while j < len(tokens) and _short(tokens[j]):
                group.append(tokens[j]); j += 1
            out.append("".join(group) if len(group) > 1 else t)
            i = j
        else:
            out.append(t); i += 1
    return out


def _remove_duplicate_prefix(tokens):
    """Strip a repeated leading prefix: [A,B,C,A,B,C,X] -> [A,B,C,X]."""
    n = len(tokens)
    for k in range(n // 2, 0, -1):
        if tokens[:k] == tokens[k:2 * k]:
            return tokens[k:]
    return tokens


def _sanitize_name(name):
    """Convert a raw ANSA part name to a clean property name (see examples
    in naming_example.txt)."""
    if not name:
        return ""
    cleaned = re.sub(r'[^A-Za-z0-9]+', '_', name).strip('_')
    if not cleaned:
        return ""
    tokens = cleaned.split('_')
    tokens = _remove_duplicate_prefix(tokens)
    tokens = _merge_short_tokens(tokens)
    while tokens and tokens[-1].lower() in _FORMAT_HINTS:
        tokens.pop()
    # Drop any later occurrence of the first token (handles part numbers
    # that appear 3+ times after a 2x prefix-dedup leaves one trailing copy).
    if len(tokens) > 1:
        first = tokens[0]
        tokens = [first] + [t for t in tokens[1:] if t != first]
    if len(tokens) > 2:
        for i, t in enumerate(tokens):
            if i > 1 and _REV_PATTERN.match(t):
                tokens.insert(1, tokens.pop(i)); break
    return "_".join(tokens)[:80]


# ─── Material application ─────────────────────────────────────────────────────

def _part_all_pids(cad_part):
    """List every PID on an ANSAPART (typically both CAD PID and FE PID)."""
    try:
        raw = cad_part.get_entity_values(GEO_DECK, ["PID"]).get("PID")
    except Exception:
        return []
    if raw is None:
        return []
    if isinstance(raw, int):
        return [raw]
    pids = []
    for token in str(raw).replace(";", ",").split(","):
        try:
            pid = int(token.strip())
            if pid > 0:
                pids.append(pid)
        except ValueError:
            pass
    return pids


def _rename_prop(pid, name):
    """Set Name on this PID's property card in both NASTRAN and LSDYNA decks."""
    for deck in (FE_DECK, MAT_DECK):
        prop = base.GetEntity(deck, "__PROPERTIES__", pid)
        if prop is not None:
            base.SetEntityCardValues(deck, prop, {"Name": name})


def apply_material(cad_part, elems, geo_mat=None, part_name=None):
    """Apply MID and/or Name to the LSDYNA properties this part owns,
    gated by the SET_MATERIALS / SET_PROPERTY_NAMES flags at the top.

    Returns (ok, message, modified_props) where modified_props are LSDYNA
    property entities, fed to SynchronizeMaterials at the end."""
    if not SET_MATERIALS and not SET_PROPERTY_NAMES:
        return True, "(material + naming both disabled)", []
    if SET_MATERIALS and geo_mat is None:
        return False, "no material on CAD part", []
    mat_id = geo_mat._id if (SET_MATERIALS and geo_mat is not None) else None
    name   = (_sanitize_name(part_name)
              if (SET_PROPERTY_NAMES and part_name) else None)

    modified, seen = [], set()

    def _touch(prop, pid):
        if SET_MATERIALS and mat_id is not None:
            base.SetEntityCardValues(MAT_DECK, prop, {"MID": mat_id})
        if name:
            _rename_prop(pid, name)
        modified.append(prop)

    for elem in elems:
        prop = base.GetEntityCardValues(MAT_DECK, elem, ["__prop__"]).get("__prop__")
        if prop is None:
            continue
        if not hasattr(prop, "_id"):
            try:
                prop = base.GetEntity(MAT_DECK, "__PROPERTIES__", int(prop))
            except (TypeError, ValueError):
                prop = None
        if prop is None or prop._id in seen:
            continue
        seen.add(prop._id)
        _touch(prop, prop._id)

    for pid in _part_all_pids(cad_part):
        if pid in seen:
            continue
        prop = base.GetEntity(MAT_DECK, "__PROPERTIES__", pid)
        if prop is None:
            continue
        seen.add(pid)
        _touch(prop, pid)

    if not modified:
        return False, "no LSDYNA properties found", []
    bits = []
    if SET_MATERIALS:      bits.append(f"mat ID={mat_id}")
    if SET_PROPERTY_NAMES: bits.append(f"name='{name}'")
    bits.append(f"{len(modified)} prop(s)")
    return True, "  ".join(bits), modified


# ─── Step 5: Vote-based matching, assignment, and material copy ───────────────

def _vote(elements, grid):
    votes    = defaultdict(int)
    min_dist = float("inf")
    seen     = set()
    for elem in elements:
        for node in _node_entities(elem):
            if node._id not in seen:
                seen.add(node._id)
                pos = node.position
                if pos is None:
                    continue
                cad_part, dist = grid.nearest(pos)
                if cad_part is not None:
                    votes[cad_part] += 1
                    if dist < min_dist:
                        min_dist = dist
    return votes, min_dist


def assign_clusters(clusters, grid, part_names, cad_materials):
    results = []
    for i, elems in enumerate(clusters):
        votes, min_dist = _vote(elems, grid)
        if not votes:
            results.append((i, elems, None, 0.0, min_dist))
            continue
        best       = max(votes, key=votes.get)
        confidence = votes[best] / max(sum(votes.values()), 1)
        results.append((i, elems, best, confidence, min_dist))

    # Detect part conflicts
    cad_usage = defaultdict(list)
    for i, elems, cad_part, *_ in results:
        if cad_part is not None:
            cad_usage[cad_part].append(i)
    for cad_part, idxs in cad_usage.items():
        if len(idxs) > 1:
            name  = part_names[cad_part]
            sizes = [len(results[j][1]) for j in idxs]
            print(f"  [WARN] {len(idxs)} clusters → '{name}' "
                  f"(element counts: {sizes})")

    total, summary, all_modified_props = 0, {}, []
    for i, elems, cad_part, confidence, min_dist in results:
        if cad_part is None:
            print(f"  [WARN] Cluster {i+1:>3}: no votes — skipped.")
            continue

        name = part_names[cad_part]
        dist_flag = (f"  [WARN dist={min_dist:.1f}]"
                     if min_dist > DISTANCE_WARNING else "")
        print(f"  Cluster {i+1:>3} ({len(elems):>7,} elems)"
              f" → '{name}'  conf={confidence*100:.0f}%{dist_flag}")

        if DRY_RUN:
            summary[name] = summary.get(name, 0) + len(elems)
            continue

        part_ok = bool(base.SetEntityPart(elems, cad_part))
        if not part_ok:
            print(f"  [ERROR] SetEntityPart failed for cluster {i+1}.")

        if part_ok:
            total += len(elems)
            summary[name] = summary.get(name, 0) + len(elems)

            geo_mat = cad_materials.get(cad_part)
            ok, msg, modified_props = apply_material(cad_part, elems, geo_mat,
                                                      part_name=name)
            if ok:
                print(f"           {msg}")
                all_modified_props.extend(modified_props)
            elif SET_MATERIALS or SET_PROPERTY_NAMES:
                print(f"           [WARN] {msg}")

    return total, summary, all_modified_props


# ─── Main ─────────────────────────────────────────────────────────────────────

def main():
    sep = "─" * 62
    tag = "  [DRY RUN]" if DRY_RUN else ""

    if not (SET_PART or SET_MATERIALS or SET_PROPERTY_NAMES):
        print(sep)
        print("  All three flags disabled "
              "(SET_PART / SET_MATERIALS / SET_PROPERTY_NAMES = F).")
        print("  Nothing to do — exiting.")
        print(sep)
        return

    # NOTE: SetCurrentDeck(LSDYNA) is called at the END, not here.
    # Calling it early causes CreateFreeMesh to put temp shells in the LSDYNA
    # deck, making CollectEntities(NASTRAN, face, "SHELL") return nothing.
    # All API calls in this script specify deck explicitly — the working deck
    # only needs to be LSDYNA for the user after the script completes.

    flags = (f"part={'Y' if SET_PART else 'N'}  "
             f"materials={'Y' if SET_MATERIALS else 'N'}  "
             f"names={'Y' if SET_PROPERTY_NAMES else 'N'}")
    print(sep)
    print(f"  Set Part + Material + Name  (LS-DYNA working deck){tag}")
    print(f"  {flags}")
    print(sep)

    # ── Fast path: SET_PART=F ──────────────────────────────────────────────
    # No need to collect FE elements, cluster them, or build reference
    # clouds — those exist solely to figure out which cluster maps to which
    # CAD part for SetEntityPart.  Material and name writes only need the
    # CAD part itself + its PIDs.
    if not SET_PART:
        print(f"\n[1/2] Identifying CAD parts …")
        cad_parts = collect_cad_parts()
        if not cad_parts:
            print("  ERROR: no CAD parts with faces found."); return
        part_names = {p: d["name"] for p, d in cad_parts.items()}
        print(f"  {len(cad_parts)} CAD part(s):")
        for data in cad_parts.values():
            print(f"    '{data['name']}':  {len(data['faces'])} faces")

        # We need cad_materials (each part -> LSDYNA material) for the MID
        # writes.  build_reference_clouds collects this but also does free-
        # meshing we don't need; just read mat_from_entities directly.
        cad_materials = {}
        for part in cad_parts:
            mats = base.CollectEntities(MAT_DECK, part, "__MATERIALS__",
                                         mat_from_entities=True) or []
            cad_materials[part] = mats[0] if mats else None

        print(f"\n[2/2] Writing properties …")
        total, summary, modified_props = 0, {}, []
        for part, data in cad_parts.items():
            name = part_names[part]
            ok, msg, mods = apply_material(part, [], cad_materials.get(part),
                                            part_name=name)
            print(f"  '{name}': {msg}")
            if ok:
                summary[name] = len(mods)
                modified_props.extend(mods)

        if modified_props and not DRY_RUN and SET_MATERIALS:
            print(f"\n  SynchronizeMaterials (LSDYNA → all decks): "
                  f"{len(modified_props)} prop(s) …")
            base.SynchronizeMaterials(
                modified_props, MAT_DECK, 0,
                True, False, False, True,
            )

        base.SetCurrentDeck(constants.LSDYNA)
        print(f"\n{sep}")
        print("  Results (no SetEntityPart):")
        for name, count in sorted(summary.items()):
            print(f"    {name:<46}  {count:>4} prop(s)")
        print(sep)
        print("  ANSA working deck switched to LS-DYNA.")
        return

    # ── Full pipeline (SET_PART=T) ─────────────────────────────────────────
    print(f"\n[1/5] Collecting FE elements …")
    all_fe_elems = collect_all_fe_elements()
    if not all_fe_elems:
        print("  ERROR: no FE elements found."); return
    existing_fe_ids = {e._id for e in all_fe_elems}
    print(f"  {len(all_fe_elems):,} element(s)  |  snapshot: {len(existing_fe_ids):,} IDs saved.")

    print(f"\n[2/5] Finding connected FE clusters via Neighb …")
    clusters = find_fe_clusters(all_fe_elems)
    print(f"  {len(clusters)} connected cluster(s):")
    for i, elems in enumerate(clusters):
        print(f"    Cluster {i+1}: {len(elems):,} elements")

    print(f"\n[3/5] Identifying CAD parts …")
    cad_parts = collect_cad_parts()
    if not cad_parts:
        print("  ERROR: no CAD parts with faces found."); return
    part_names = {p: d["name"] for p, d in cad_parts.items()}
    print(f"  {len(cad_parts)} CAD part(s):")
    for data in cad_parts.values():
        print(f"    '{data['name']}':  {len(data['faces'])} faces")

    if len(clusters) != len(cad_parts):
        print(f"\n  [NOTE] {len(clusters)} FE cluster(s) vs "
              f"{len(cad_parts)} CAD part(s) — not 1:1.")

    print(f"\n[4/5] Building reference clouds via temp free mesh "
          f"(size={TEMP_MESH_SIZE}) …")
    clouds, cad_materials = build_reference_clouds(cad_parts, existing_fe_ids)
    if not clouds:
        print("  ERROR: no reference points generated."); return
    grid = build_spatial_index(clouds)

    print(f"\n[5/5] Assigning {len(clusters)} cluster(s) + copying materials …")
    total, summary, modified_props = assign_clusters(
        clusters, grid, part_names, cad_materials)

    if all_fe_elems:
        base.Or(entities=all_fe_elems)

    # One call to ANSA's built-in cross-deck sync.  Source is LSDYNA (where
    # apply_material just wrote integer MIDs); ANSA materializes the deck-
    # native materials and wires the NASTRAN / ABAQUS / PAM-CRASH property
    # cards.  delete_released=True drops orphans left over from the rewrite.
    if modified_props and not DRY_RUN and SET_MATERIALS:
        print(f"\n  SynchronizeMaterials (LSDYNA → all decks): "
              f"{len(modified_props)} prop(s) …")
        base.SynchronizeMaterials(
            modified_props, MAT_DECK, 0,
            True, False, False, True,
        )

    # Switch working deck to LS-DYNA NOW — after all visibility/mesh ops are done
    base.SetCurrentDeck(constants.LSDYNA)

    action = "would be assigned" if DRY_RUN else "assigned"
    print(f"\n{sep}")
    print("  Results:")
    for name, count in sorted(summary.items()):
        print(f"    {name:<46}  {count:>8,} elements")
    print(f"\n  Total: {total:,} elements {action} to {len(summary)} part(s).")
    print(sep)
    print("  ANSA working deck switched to LS-DYNA.")


if __name__ == "__main__":
    main()
