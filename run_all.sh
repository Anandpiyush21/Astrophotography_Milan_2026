#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
for s in 01_prep 02_adaptive_smooth 03_compose 04_figures 05_report; do
  echo "$s"
  .venv/bin/python pipeline/$s.py
done
