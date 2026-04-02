#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import pathlib
import shutil
from typing import List


def collect_related_files(csv_path: pathlib.Path) -> List[pathlib.Path]:
    candidates = [
        csv_path,
        csv_path.with_suffix('.h5'),
        csv_path.with_suffix('.vtk'),
        csv_path.with_suffix('.benchmark.txt'),
        csv_path.with_suffix('.benchmark.csv'),
        csv_path.with_suffix('.graph.txt'),
        csv_path.with_suffix('.graph_matrix.csv'),
        csv_path.with_suffix('.graph_matrix.json'),
        csv_path.with_suffix('.field_adapter.txt'),
        csv_path.with_suffix('.field_adapter.json'),
        csv_path.with_suffix('.boundary_mapping.json'),
        csv_path.with_suffix('.field_request.json'),
        csv_path.with_suffix('.field_result.json'),
        csv_path.with_suffix('.field_result_template.json'),
        csv_path.with_suffix('.volume_stub.json'),
        csv_path.with_suffix('.volume_mesh_stub.json'),
        csv_path.with_suffix('.field_bridge_manifest.json'),
        csv_path.with_suffix('.surface_volume_projection.json'),
        csv_path.with_suffix('.volumetric_adapter.json'),
        csv_path.with_suffix('.volume_request.json'),
        csv_path.with_suffix('.volume_result.json'),
        csv_path.with_suffix('.volume_result_template.json'),
        csv_path.with_suffix('.volume_history.json'),
    ]
    return [path for path in candidates if path.exists()]


def main() -> int:
    parser = argparse.ArgumentParser(
        description='Bundle a SCDAT surface-charging case and its sidecar outputs into one directory.'
    )
    parser.add_argument('csv_path')
    parser.add_argument('output_dir')
    parser.add_argument('--copy', action='store_true', help='Copy files instead of moving them.')
    args = parser.parse_args()

    csv_path = pathlib.Path(args.csv_path)
    output_dir = pathlib.Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    files = collect_related_files(csv_path)
    manifest = {
        'source_csv': str(csv_path),
        'mode': 'copy' if args.copy else 'move',
        'files': [],
    }

    for source in files:
        destination = output_dir / source.name
        if args.copy:
            shutil.copy2(source, destination)
        else:
            if destination.exists():
                destination.unlink()
            shutil.move(str(source), str(destination))
        manifest['files'].append({'name': destination.name, 'path': str(destination)})
        print(f'bundled={destination}')

    manifest_path = output_dir / 'surface_case_manifest.json'
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding='utf-8')
    print(f'manifest={manifest_path}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
