# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

#!/usr/bin/env python3
import os
import json
import re
from pathlib import Path
from datetime import datetime

# Configuration
SKILLS_ROOT = Path("Skills")
OUTPUT_FILE = Path("skills_catalog.json")
REQUIRED_FIELDS = ["name", "description", "measurable_outcome", "allowed-tools"]

def parse_frontmatter(content):
    """
    Parses YAML-like frontmatter from a markdown file.
    Returns a dictionary of metadata and the remaining markdown content.
    """
    match = re.match(r'^---\s*\n(.*?)\n---\s*\n(.*)', content, re.DOTALL)
    if not match:
        return None, content

    yaml_text = match.group(1)
    body = match.group(2)
    
    metadata = {}
    current_key = None
    
    # Simple line-by-line parser for basic YAML structure
    for line in yaml_text.split('\n'):
        line = line.rstrip()
        if not line or line.startswith('#'):
            continue
            
        # List item (indentation detection is primitive but sufficient for this schema)
        if line.strip().startswith('- '):
            if current_key:
                val = line.strip()[2:].strip()
                if not isinstance(metadata[current_key], list):
                    metadata[current_key] = []
                metadata[current_key].append(val)
            continue

        # Key-value pair
        if ':' in line:
            key, val = line.split(':', 1)
            key = key.strip()
            val = val.strip()
            
            # Remove quotes if present
            if len(val) >= 2 and ((val.startswith('"') and val.endswith('"')) or (val.startswith("'") and val.endswith("'"))):
                val = val[1:-1]
            
            if not val:
                current_key = key
                metadata[key] = [] # Assume list start or empty object
            else:
                current_key = key
                metadata[key] = val
                
    return metadata, body

def scan_skills(root_dir):
    """
    Recursively finds all SKILL.md files and parses them.
    """
    catalog = []
    errors = []
    
    print(f"Scanning {root_dir} for skills...")
    
    for path in Path(root_dir).rglob("SKILL.md"):
        try:
            with open(path, 'r', encoding='utf-8') as f:
                content = f.read()
                
            metadata, body = parse_frontmatter(content)
            
            if not metadata:
                errors.append(f"Missing frontmatter in {path}")
                continue
                
            # Basic Validation
            missing = [field for field in REQUIRED_FIELDS if field not in metadata]
            if missing:
                errors.append(f"Missing fields {missing} in {path}")
                # We add it anyway, but flagged
                metadata['_validation_error'] = f"Missing: {missing}"
            
            # Add file info
            metadata['file_path'] = str(path)
            metadata['last_modified'] = datetime.fromtimestamp(path.stat().st_mtime).isoformat()
            
            # Add to catalog
            catalog.append(metadata)
            
        except Exception as e:
            errors.append(f"Error processing {path}: {str(e)}")

    return catalog, errors

def main():
    if not SKILLS_ROOT.exists():
        print(f"Error: Skills directory '{SKILLS_ROOT}' not found.")
        return

    catalog, errors = scan_skills(SKILLS_ROOT)
    
    # Save to JSON
    with open(OUTPUT_FILE, 'w', encoding='utf-8') as f:
        json.dump({"generated_at": datetime.now().isoformat(), "skills": catalog}, f, indent=2)
        
    print(f"\nCatalog generation complete.")
    print(f"Total Skills Found: {len(catalog)}")
    print(f"Catalog saved to: {OUTPUT_FILE.absolute()}")
    
    if errors:
        print("\nErrors/Warnings:")
        for err in errors:
            print(f"- {err}")
    else:
        print("\nNo errors found. All skills valid.")

if __name__ == "__main__":
    main()

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
