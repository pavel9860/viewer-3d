"""Export CATIA parts and one assembly to binary PLY with the project's own reader (cadcli), for the real-data
kernel workloads. Usage: python3 -I tools/export_parts.py SRC_DIR OUT_DIR CADCLI [ASSEMBLY.CATProduct]

Writes OUT_DIR/parts/NNNN.ply + parts/manifest.tsv (index, source), and OUT_DIR/asm.ply when an assembly is given.
Runs 2 exports in parallel with an 8 GB address-space cap per process.
"""
import concurrent.futures as cf
import os
import resource
import subprocess
import sys


def limit():
    resource.setrlimit(resource.RLIMIT_AS, (8 << 30, 8 << 30))


def main():
    src, out, cli = sys.argv[1:4]
    os.makedirs(os.path.join(out, "parts"), exist_ok=True)
    files = sorted(os.path.join(d, f) for d, _, fs in os.walk(src) for f in fs if f.lower().endswith(".catpart"))

    def run(i):
        dst = os.path.join(out, "parts", f"{i:04d}.ply")
        r = subprocess.run([cli, "ply", files[i], dst], capture_output=True, text=True, preexec_fn=limit)
        return i, r.returncode

    ok = []
    with cf.ThreadPoolExecutor(2) as ex:
        for i, rc in ex.map(run, range(len(files))):
            if rc == 0:
                ok.append(f"{i:04d}.ply\t{files[i]}")
    with open(os.path.join(out, "parts", "manifest.tsv"), "w") as f:
        f.write("\n".join(ok) + "\n")
    if len(sys.argv) > 4:
        subprocess.run([cli, "ply", sys.argv[4], os.path.join(out, "asm.ply")], check=True, preexec_fn=limit)
    print(f"{len(files)} CATParts, {len(ok)} exported")


if __name__ == "__main__":
    main()
