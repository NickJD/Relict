import logging
import subprocess

logger = logging.getLogger(__name__)


def run_cmd(cmd):
    logger.debug("[CMD] %s", cmd)
    result = subprocess.run(cmd, shell=True)

    if result.returncode != 0:
        raise RuntimeError(f"Command failed: {cmd}")

