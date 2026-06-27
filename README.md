## **Set Part, Material, and Property Name by Free‑Mesh Proximity (ANSA v22.1.5)**

A Python automation script for **ANSA** that uses **free‑mesh proximity voting** to automatically align FE clusters with their corresponding CAD parts. It then applies **part assignment, material IDs, and clean property names** across multiple solver decks.

---

### **What This Script Does**
- **Set Part** → Reassign FE elements to the matched CAD `ANSAPART`  
- **Set Materials** → Write CAD material IDs onto LS‑DYNA properties, then synchronize across NASTRAN / ABAQUS / PAM‑CRASH  
- **Set Property Names** → Copy sanitized CAD part names onto property cards (NASTRAN + LS‑DYNA)

Each operation is gated by a flag (`SET_PART`, `SET_MATERIALS`, `SET_PROPERTY_NAMES`) so you can enable or disable them independently.

---

### **How It Works**
1. **Cluster Detection** → Finds connected FE clusters via `Neighb()` flood‑fill  
2. **CAD Part Identification** → Collects CAD parts (`ANSAPART`) with FACE geometry  
3. **Reference Clouds** → Builds free‑mesh point clouds for proximity matching  
4. **Voting System** → Matches FE clusters to CAD parts based on nearest‑point votes  
5. **Assignment & Sync** → Applies part assignment, material IDs, and sanitized names; synchronizes materials across decks

---

### **Features**
- Works across **NASTRAN** (FE + CAD geometry) and **LS‑DYNA** (materials) decks  
- **SpatialGrid** hashing for fast nearest‑point lookup  
- **Name sanitization** to clean CAD part names (removes duplicates, format hints, merges short tokens)  
- **Dry‑run mode** for previewing changes without modifying the model  
- Clear console reporting with warnings, confidence scores, and summaries  

---

### **Usage**
Inside ANSA:  
`Scripts → Run Script → Select this file`

Flags at the top of the script control behavior:  

```python
SET_PART            = True   # Assign FE elements to CAD parts
SET_MATERIALS       = True   # Apply CAD material IDs
SET_PROPERTY_NAMES  = True   # Copy sanitized names to properties
```

---

### **Example Output**
```
──────────────────────────────────────────────
  Set Part + Material + Name  (LS-DYNA working deck)
  part=Y  materials=Y  names=Y
──────────────────────────────────────────────

[1/5] Collecting FE elements …
  12,345 element(s)  |  snapshot: 12,345 IDs saved.

[2/5] Finding connected FE clusters via Neighb …
  8 connected cluster(s):
    Cluster 1: 1,234 elements
    Cluster 2: 2,345 elements
    ...

[5/5] Assigning clusters + copying materials …
  Cluster 3 (  456 elems) → 'A5304012000_DICASTAL_03_VA_7X23'  conf=92%
           mat ID=101  name='A5304012000_DICASTAL_03_VA_7X23'  2 prop(s)

Results:
    A5304012000_DICASTAL_03_VA_7X23             456 elements
    HALTER_BATTERIE                             789 elements

Total: 12,345 elements assigned to 8 part(s).
──────────────────────────────────────────────
  ANSA working deck switched to LS-DYNA.
```

---

### **Summary**
This script is ideal for:
- Cleaning up imported FE models  
- Ensuring **1:1 mapping** between FE clusters and CAD parts  
- Propagating **materials and names** consistently across decks  
- Automating repetitive preprocessing tasks in ANSA  
