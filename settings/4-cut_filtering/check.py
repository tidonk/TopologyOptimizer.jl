import os
import re

base_dir = os.path.dirname(os.path.abspath(__file__))
folders = [
    "1-cutfiltering-pb-spsel-rel",
    "2-cutfiltering-pb-spsel-abs",
    "3-cutfiltering-abs-spsel-abs",
    "4-cutfiltering-rel-spsel-rel",
    "5-cutfiltering-rel-spsel-abs",
    "6-cutfiltering-abs-spsel-rel"
]

regexes = {
    "1-cutfiltering-pb-spsel-rel": r"^sp([0-9]{3})\.toml$",
    "2-cutfiltering-pb-spsel-abs": r"^sp([0-9]{2})\.toml$",
    "3-cutfiltering-abs-spsel-abs": r"^cf([0-9]{2})_sp([0-9]{2})\.toml$",
    "4-cutfiltering-rel-spsel-rel": r"^cf([0-9]{3})_sp([0-9]{3})\.toml$",
    "5-cutfiltering-rel-spsel-abs": r"^cf([0-9]{3})_sp([0-9]{2})\.toml$",
    "6-cutfiltering-abs-spsel-rel": r"^cf([0-9]{2})_sp([0-9]{3})\.toml$"
}

def check_file(folder, filename, cf_val, sp_val):
    path = os.path.join(base_dir, folder, filename)
    with open(path, 'r') as f:
        content = f.read()
    
    mismatches = []
    
    # SP check
    if folder in ["1-cutfiltering-pb-spsel-rel", "4-cutfiltering-rel-spsel-rel", "6-cutfiltering-abs-spsel-rel"]:
        # rel
        val = float(sp_val) / 100.0
        pattern = rf"STOP_AFTER_N_HITS_REL\s*=\s*{val:.1f}|{val:.2f}|{val:.3f}"
        if not re.search(rf"STOP_AFTER_N_HITS_REL\s*=\s*{val:.2f}", content) and not re.search(rf"STOP_AFTER_N_HITS_REL\s*=\s*{val:.1f}", content):
             # Handle 0.05 as .05 or 0.05
             search_val = f"{val:.2f}".lstrip('0') if val < 1 and val > 0 else f"{val:.2f}"
             if not (re.search(rf"STOP_AFTER_N_HITS_REL\s*=\s*{val:.2f}", content) or re.search(rf"STOP_AFTER_N_HITS_REL\s*=\s*{val:g}", content)):
                mismatches.append(f"SP mismatch: expected STOP_AFTER_N_HITS_REL = {val}")
    else:
        # abs
        val = int(sp_val)
        if not re.search(rf"STOP_AFTER_N_HITS\s*=\s*{val}", content):
            mismatches.append(f"SP mismatch: expected STOP_AFTER_N_HITS = {val}")

    # CF check
    if folder in ["3-cutfiltering-abs-spsel-abs", "6-cutfiltering-abs-spsel-rel"]:
        # abs
        val = int(cf_val)
        if not re.search(r"LIMIT_MODE\s*=\s*:abs_num", content):
            mismatches.append("CF mismatch: expected LIMIT_MODE = :abs_num")
        if not re.search(rf"LIMIT_ABS_NUM\s*=\s*{val}", content):
            mismatches.append(f"CF mismatch: expected LIMIT_ABS_NUM = {val}")
    elif folder in ["4-cutfiltering-rel-spsel-rel", "5-cutfiltering-rel-spsel-abs"]:
        # rel
        val = float(cf_val) / 100.0
        if not re.search(r"LIMIT_MODE\s*=\s*:rel_num", content):
            mismatches.append("CF mismatch: expected LIMIT_MODE = :rel_num")
        if not (re.search(rf"LIMIT_REL_ALPHA\s*=\s*{val:.2f}", content) or re.search(rf"LIMIT_REL_ALPHA\s*=\s*{val:g}", content)):
            mismatches.append(f"CF mismatch: expected LIMIT_REL_ALPHA = {val}")

    return mismatches

for folder in folders:
    folder_path = os.path.join(base_dir, folder)
    if not os.path.exists(folder_path):
        print(f"Folder {folder} does not exist")
        continue
    
    regex = regexes[folder]
    for filename in os.listdir(folder_path):
        if not filename.endswith(".toml"): continue
        match = re.match(regex, filename)
        if not match:
            print(f"Filename mismatch in {folder}: {filename}")
            continue
        
        groups = match.groups()
        cf_val = None
        sp_val = None
        if len(groups) == 1:
            sp_val = groups[0]
        else:
            cf_val, sp_val = groups
        
        mismatches = check_file(folder, filename, cf_val, sp_val)
        for m in mismatches:
            print(f"Content mismatch in {folder}/{filename}: {m}")

readme_path = os.path.join(base_dir, "README.md")
if os.path.exists(readme_path):
    with open(readme_path, 'r') as f:
        for line in f:
            if "cf001_sp05" in line or "cf02_sp005" in line:
                print(f"README: {line.strip()}")
