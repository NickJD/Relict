import subprocess


def run_cmd(cmd):
    print(f"[CMD] {cmd}")
    result = subprocess.run(cmd, shell=True)

    if result.returncode != 0:
        raise RuntimeError(f"Command failed: {cmd}")
