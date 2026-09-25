#!/usr/bin/env python3
# Parse approxmc run outputs (out-approxmc-*/*.out_approxmc[.xz] + .timeout_approxmc[.xz]) into data.sqlite3

import argparse
import collections
import glob
import lzma
import math
import multiprocessing
import os
import re
import sqlite3

DB = "data.sqlite3"
ANSI_RE = re.compile(r"\x1b\[[0-9;]*m|\[0m")

COLS = [
    ("solver", "TEXT NOT NULL"),
    ("dirname", "TEXT NOT NULL"),
    ("fname", "TEXT NOT NULL"),
    ("ver", "TEXT"),
    ("call", "TEXT"),
    ("solved", "INT"),
    ("unsat", "INT"),
    ("time", "FLOAT"),       # user time, NULL if unsolved
    ("timeout_t", "FLOAT"),  # user time even when unsolved
    ("wall_t", "FLOAT"),
    ("mem_MB", "FLOAT"),
    ("signal", "INT"),
    ("exit_status", "INT"),
    ("mem_out", "INT"),
    ("errored", "INT"),
    ("mc_log10", "FLOAT"),
    ("orig_sampl_sz", "INT"),
    ("arjun_sampl_sz", "INT"),
    ("sampl_sz", "INT"),
    ("threshold", "INT"),
    ("sparse", "INT"),
    ("rounds", "INT"),
    ("max_hashes", "INT"),
    ("arjun_time", "FLOAT"),
    ("appmc_time", "FLOAT"),
    ("total_time", "FLOAT"),
]


def xopen(fname):
    if fname.endswith(".xz"):
        return lzma.open(fname, "rt", errors="replace")
    return open(fname, "r", errors="replace")


def parse_wall(s):
    secs = 0.0
    for part in s.split(":"):
        secs = secs * 60 + float(part)
    return secs


def parse_timeout(fname):
    r: dict = {}
    with xopen(fname) as f:
        for line in f:
            line = line.strip()
            if line.startswith("Command terminated by signal"):
                r["signal"] = int(line.split()[4])
            elif line.startswith("User time (seconds)"):
                r["timeout_t"] = float(line.split()[3])
            elif line.startswith("Elapsed (wall clock)"):
                r["wall_t"] = parse_wall(line.split()[7])
            elif line.startswith("Maximum resident set size"):
                r["mem_MB"] = int(line.split()[5]) / 1000
            elif line.startswith("Exit status:"):
                r["exit_status"] = int(line.split()[2])
    return r


def parse_out(fname):
    r: dict = {"mem_out": 0, "errored": 0, "solved": 0, "unsat": 0}
    sha = {}
    with xopen(fname) as f:
        for line in f:
            line = ANSI_RE.sub("", line).strip()
            if "bad_alloc" in line or "Cannot allocate memory" in line:
                r["mem_out"] = 1
            if line.startswith("ERROR") or "assertion fail" in line.lower():
                r["errored"] = 1

            if line.startswith("c o [appmc]"):
                if "round:" in line:
                    m = re.search(r"round:\s*(\d+)\s+hashes:\s*(\d+)", line)
                    if m:
                        r["rounds"] = max(r.get("rounds", 0), int(m.group(1)) + 1)
                        r["max_hashes"] = max(r.get("max_hashes", 0), int(m.group(2)))
                elif line.startswith("c o [appmc] Sampling set size:"):
                    r["sampl_sz"] = int(line.split()[6])
                elif line.startswith("c o [appmc] threshold set to"):
                    f_ = line.split()
                    r["threshold"] = int(f_[6])
                    r["sparse"] = int(f_[8])
                elif line.startswith("c o [appmc] ApproxMC T:"):
                    r["appmc_time"] = float(line.split()[5])
            elif line.startswith("c o [arjun] final set size:"):
                m = re.search(r"size:\s*(\d+) percent of original:\s*([\d.]+)", line)
                if m:
                    r["arjun_sampl_sz"] = int(m.group(1))
                    pct = float(m.group(2))
                    if pct > 0:
                        r["orig_sampl_sz"] = round(int(m.group(1)) * 100 / pct)
            elif line.startswith("c o [arjun] Arjun finished. T:"):
                r["arjun_time"] = float(line.split()[6])
            elif line.startswith("c o [appmc+arjun] Total time:"):
                r["total_time"] = float(line.split()[5])
            elif line.startswith("c o executed with command line:"):
                r["call"] = " ".join(line.split()[7:-1])
            elif " SHA1: " in line and line.startswith("c o "):
                name = line.split()[2]
                sha.setdefault(name, line.split()[-1][:8])
            elif line.startswith("s SATISFIABLE"):
                r["solved"] = 1
            elif line.startswith("s UNSATISFIABLE"):
                r["solved"] = 1
                r["unsat"] = 1
            elif line.startswith("s mc "):
                val = line.split()[2]
                if val.isdigit() and val != "0":
                    r["mc_log10"] = len(val) - 1 + math.log10(float(val[0] + "." + val[1:16]))
    r["ver"] = "-".join(sha[k] for k in ("ApproxMC", "Arjun", "CMS") if k in sha) or None
    return r


def key_of(path):
    dirname, base = path.split("/")[-2:]
    return dirname, base.split(".cnf")[0] + ".cnf"


def parse_one(path):
    if ".timeout_approxmc" in path:
        return key_of(path), parse_timeout(path)
    return key_of(path), parse_out(path)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--files", default="out-approxmc*/*cnf*_approxmc*")
    p.add_argument("--reparse", action="store_true", help="Delete the DB and reparse everything")
    p.add_argument("-j", type=int, default=os.cpu_count())
    args = p.parse_args()

    if args.reparse and os.path.exists(DB):
        os.remove(DB)
    con = sqlite3.connect(DB)
    con.execute("CREATE TABLE IF NOT EXISTS data (" + ", ".join(f"{c} {t}" for c, t in COLS) + ")")
    con.execute("CREATE UNIQUE INDEX IF NOT EXISTS data_idx ON data (dirname, fname)")
    done = set(con.execute("SELECT dirname, fname FROM data"))
    if done:
        print(f"Skipping {len(done)} already-parsed instances (use --reparse to redo)")

    todo = [f for f in glob.glob(args.files) if key_of(f) not in done]
    if not todo:
        print("Nothing new to parse")
        return

    runs = collections.defaultdict(dict)
    with multiprocessing.Pool(args.j) as pool:
        for key, r in pool.imap_unordered(parse_one, todo, chunksize=8):
            runs[key].update(r)

    rows = []
    for (dirname, fname), r in runs.items():
        if "solved" not in r or "timeout_t" not in r:
            print(f"WARNING: missing out or timeout file, skipping {dirname}/{fname}")
            continue
        r["solver"] = "approxmc"
        r["dirname"] = dirname
        r["fname"] = fname
        r["time"] = r["timeout_t"] if r["solved"] else None
        rows.append(tuple(r.get(c) for c, _ in COLS))

    con.executemany(f"INSERT INTO data ({','.join(c for c, _ in COLS)}) VALUES "
                    f"({','.join('?' * len(COLS))})", rows)
    con.commit()
    for d, n in sorted(collections.Counter(r[1] for r in rows).items()):
        print(f"  {d}: {n} instances")
    print(f"Inserted {len(rows)} rows into {DB}")


if __name__ == "__main__":
    main()
