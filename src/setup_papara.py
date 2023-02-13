#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Download papara binary and setup path
"""

import os
from pathlib import Path
from phyloplacement.utils import terminal_execute, extract_tar_file


parent_dir = Path(__file__).parent
papara_url = "https://cme.h-its.org/exelixis/resource/download/software/papara_nt-2.5-static_x86_64.tar.gz"
papara_dir = parent_dir / "vendor"

if not (papara_dir / "papara_static_x86_64").exists():
    terminal_execute(f"wget {papara_url}", work_dir=parent_dir)
    extract_tar_file(
        tar_file=str(parent_dir / "papara_nt-2.5-static_x86_64.tar.gz"),
        dest_dir=papara_dir,
    )
    os.remove(parent_dir / "papara_nt-2.5-static_x86_64.tar.gz")
else:
    print("Papara already downloaded and ready to use")
