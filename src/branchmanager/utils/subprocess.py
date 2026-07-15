import logging
import shlex
import subprocess
from pathlib import Path

logger = logging.getLogger(__name__)


def run_cmd(cmd, *, stdout_path=None):
    """Run an external tool without invoking a shell.

    Existing MAFFT/FastTree call sites historically used a single ``>``
    redirect. It is parsed here and opened directly; all other shell operators
    are rejected so filenames and user paths cannot become shell syntax.
    """
    tokens = shlex.split(cmd) if isinstance(cmd, str) else [str(item) for item in cmd]
    forbidden = {'|', '||', '&&', ';', '<', '>>', '2>', '2>>'}
    if any(token in forbidden for token in tokens):
        raise ValueError(f'Unsupported shell operator in command: {cmd}')
    if '>' in tokens:
        if tokens.count('>') != 1 or stdout_path is not None:
            raise ValueError(f'Invalid output redirection in command: {cmd}')
        index = tokens.index('>')
        if index == 0 or index + 2 != len(tokens):
            raise ValueError(f'Output redirection must be the final command element: {cmd}')
        stdout_path = tokens[index + 1]
        tokens = tokens[:index]
    logger.debug("[CMD] %s", shlex.join(tokens))
    output_handle = None
    try:
        if stdout_path is not None:
            destination = Path(stdout_path)
            destination.parent.mkdir(parents=True, exist_ok=True)
            output_handle = open(destination, 'wb')
        result = subprocess.run(
            tokens,
            stdout=output_handle if output_handle is not None else None,
            stderr=subprocess.PIPE,
            check=False,
        )
    except FileNotFoundError as exc:
        raise RuntimeError(f'Executable not found: {tokens[0]}') from exc
    finally:
        if output_handle is not None:
            output_handle.close()
    if result.returncode != 0:
        stderr = (result.stderr or b'').decode('utf-8', errors='replace').strip()
        detail = f': {stderr}' if stderr else ''
        raise RuntimeError(f"Command failed ({result.returncode}): {shlex.join(tokens)}{detail}")
    return result
