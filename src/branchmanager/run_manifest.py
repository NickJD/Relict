"""Reproducible, machine-readable run manifests for BranchManager workflows."""

from __future__ import annotations

import hashlib
import json
import os
import platform
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Optional

from branchmanager import __version__


STAGE_LABELS = {
    'onboarding': 'Onboarding',
    'paper_trail': 'Paper Trail',
    'merge_meeting': 'Merge Meeting',
    'filing_cabinet': 'Filing Cabinet',
    'performance_review': 'Performance Review',
    'hiring_panel': 'Hiring Panel',
    'quarterly_review': 'Quarterly Review',
    'status_meeting': 'Status Meeting',
    'records_update': 'Records Update',
    'annual_report': 'Annual Report',
    'assistant': 'Assistant to the Branch Manager',
}


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace('+00:00', 'Z')


def sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with open(path, 'rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def file_record(path: str | Path, *, role: str = 'input', required: bool = True) -> dict:
    item = Path(path).expanduser()
    resolved = item.resolve()
    record = {
        'role': role,
        'path': str(resolved),
        'required': bool(required),
        'exists': resolved.exists(),
        'kind': 'file' if resolved.is_file() else ('directory' if resolved.is_dir() else 'missing'),
    }
    if resolved.is_file():
        stat = resolved.stat()
        record.update({'bytes': stat.st_size, 'sha256': sha256_file(resolved)})
    elif resolved.is_dir():
        digest = hashlib.sha256()
        count = 0
        total_bytes = 0
        for child in sorted(item for item in resolved.rglob('*') if item.is_file()):
            relative = child.relative_to(resolved).as_posix()
            child_hash = sha256_file(child)
            digest.update(relative.encode('utf-8'))
            digest.update(b'\0')
            digest.update(child_hash.encode('ascii'))
            digest.update(b'\n')
            count += 1
            total_bytes += child.stat().st_size
        record.update({
            'bytes': total_bytes,
            'file_count': count,
            'sha256': digest.hexdigest(),
        })
    return record


def _command_version(candidates: Iterable[str], version_args: Iterable[str] = ('--version',)) -> dict:
    for candidate in candidates:
        binary = shutil.which(candidate)
        if not binary:
            continue
        try:
            result = subprocess.run(
                [binary, *version_args], capture_output=True, text=True, timeout=20, check=False,
            )
            text = (result.stdout or result.stderr or '').strip().splitlines()
            version = text[0] if text else f'exit={result.returncode}'
        except Exception as exc:
            version = f'could_not_query: {exc}'
        return {'available': True, 'path': binary, 'version': version}
    return {'available': False, 'path': '', 'version': ''}


def external_tool_versions() -> dict:
    return {
        'vsearch': _command_version(('vsearch',)),
        'mafft': _command_version(('mafft',), ('--version',)),
        'fasttree': _command_version(('FastTree', 'fasttree'), ('-help',)),
        'iqtree': _command_version(('iqtree3', 'iqtree2', 'iqtree'), ('--version',)),
    }


def git_revision(start: str | Path) -> str:
    try:
        result = subprocess.run(
            ['git', '-C', str(Path(start).resolve()), 'rev-parse', 'HEAD'],
            capture_output=True, text=True, timeout=10, check=False,
        )
        return result.stdout.strip() if result.returncode == 0 else ''
    except Exception:
        return ''


class RunManifest:
    """Incrementally persist a run's state, inputs, tools, stages, and outputs."""

    def __init__(
        self,
        outdir: str | Path,
        workflow: str,
        *,
        argv: Optional[Iterable[str]] = None,
        project_root: Optional[str | Path] = None,
    ):
        self.outdir = Path(outdir)
        self.outdir.mkdir(parents=True, exist_ok=True)
        self.json_path = self.outdir / 'run_manifest.json'
        self.tsv_path = self.outdir / 'run_manifest.tsv'
        self.data = {
            'schema_version': '1.0',
            'workflow': workflow,
            'stage_label': STAGE_LABELS.get(workflow, workflow.replace('_', ' ').title()),
            'status': 'RUNNING',
            'started_at': utc_now(),
            'completed_at': '',
            'branchmanager_version': __version__,
            'git_revision': git_revision(project_root or Path(__file__).parents[2]),
            'command': list(argv if argv is not None else sys.argv),
            'environment': {
                'python': sys.version.replace('\n', ' '),
                'executable': sys.executable,
                'platform': platform.platform(),
                'hostname': platform.node(),
                'cwd': os.getcwd(),
            },
            'external_tools': external_tool_versions(),
            'inputs': [],
            'stages': [],
            'outputs': [],
            'warnings': [],
            'error': '',
        }
        self.write()

    def add_input(self, path: str | Path, *, role: str = 'input', required: bool = True) -> dict:
        record = file_record(path, role=role, required=required)
        self.data['inputs'].append(record)
        self.write()
        return record

    def add_stage(self, name: str, status: str, *, detail: str = '') -> None:
        self.data['stages'].append({
            'name': name,
            'label': STAGE_LABELS.get(name, name.replace('_', ' ').title()),
            'status': status,
            'time': utc_now(),
            'detail': detail,
        })
        self.write()

    def add_output(self, path: str | Path, *, role: str = 'output', required: bool = True) -> dict:
        record = file_record(path, role=role, required=required)
        self.data['outputs'].append(record)
        self.write()
        return record

    def warn(self, message: object) -> None:
        self.data['warnings'].append(str(message))
        self.write()

    def finish(self, status: str, *, error: object = '') -> None:
        self.data['status'] = str(status).upper()
        self.data['completed_at'] = utc_now()
        self.data['error'] = str(error or '')
        self.write()

    def verify_required_outputs(self) -> list[str]:
        missing = [
            record['path'] for record in self.data['outputs']
            if record.get('required') and not Path(record['path']).is_file()
        ]
        return missing

    def write(self) -> None:
        temporary = self.json_path.with_suffix('.json.tmp')
        temporary.write_text(json.dumps(self.data, indent=2, sort_keys=False) + '\n')
        os.replace(temporary, self.json_path)
        rows = [('Field', 'Value')]
        for key in (
            'schema_version', 'workflow', 'stage_label', 'status', 'started_at', 'completed_at',
            'branchmanager_version', 'git_revision', 'error',
        ):
            rows.append((key, self.data.get(key, '')))
        rows.append(('command', ' '.join(str(item) for item in self.data.get('command', []))))
        for stage in self.data.get('stages', []):
            rows.append((f"stage:{stage['name']}", f"{stage['status']} | {stage.get('detail', '')}"))
        for record in self.data.get('inputs', []):
            rows.append((f"input:{record.get('role', 'input')}", f"{record.get('path')} | sha256={record.get('sha256', 'NA')}"))
        for record in self.data.get('outputs', []):
            rows.append((f"output:{record.get('role', 'output')}", f"{record.get('path')} | exists={record.get('exists')} | sha256={record.get('sha256', 'NA')}"))
        with open(self.tsv_path, 'w') as handle:
            for key, value in rows:
                clean = str(value).replace('\t', ' ').replace('\n', ' ')
                handle.write(f'{key}\t{clean}\n')


def load_manifest(path: str | Path) -> dict:
    with open(path) as handle:
        return json.load(handle)
