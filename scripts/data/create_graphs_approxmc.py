#!/usr/bin/env python3
# Summary tables, gained/lost comparisons and CDF/scatter plots from data.sqlite3

import argparse
import base64
import itertools
import math
import os
import re
import sqlite3

BLUE = "\033[94m"
GREEN = "\033[92m"
RED = "\033[91m"
RESET = "\033[0m"

TMP_DIR = "tmp"
DB = "data.sqlite3"
TIMEOUT = 3600
# two (eps=0.8, delta=0.2) estimates can each be off by 1.8x in opposite directions
MC_TOL = math.log10(1.8 * 1.8)


def short(d):
    return d.replace("out-approxmc-", "")


def print_header(title):
    line = "*" * (len(title) + 8)
    print(f"\n{GREEN}{line}\n*** {title.upper()} ***\n{line}{RESET}")


def print_table(headers, rows):
    rows = [["" if v is None else str(v) for v in r] for r in rows]
    widths = [max(len(h), max((len(r[i]) for r in rows), default=0)) for i, h in enumerate(headers)]
    sep = "+-" + "-+-".join("-" * w for w in widths) + "-+"
    fmt = "| " + " | ".join(f"{{:<{w}}}" for w in widths) + " |"
    print(sep)
    print(fmt.format(*headers))
    print(sep)
    for r in rows:
        print(fmt.format(*r))
    if not rows:
        print(f"{GREEN}EMPTY{RESET}")
    print(sep)


def ft(v, digits=1):
    return None if v is None else f"{v:.{digits}f}"


def display_png(png_file):
    if not os.path.exists(png_file):
        return
    with open(png_file, "rb") as fh:
        data = fh.read()
    w = int.from_bytes(data[16:20], "big")
    h = int.from_bytes(data[20:24], "big")
    print(f"\033]1337;File=inline=1;width={w}px;height={h}px:{base64.b64encode(data).decode()}\a")


def run_gnuplot(gp_file, pdf_file, png_file):
    for f in (pdf_file, png_file):
        if os.path.exists(f):
            os.unlink(f)
    os.system(f"gnuplot {gp_file}")
    print(f"  PDF: {pdf_file}  PNG: {png_file}")
    display_png(png_file)


def gp_str(s):
    return s.replace('"', '\\"').replace("_", "\\\\_")


def get_dirs(con, prefixes):
    dirs = [r[0] for r in con.execute("SELECT DISTINCT dirname FROM data ORDER BY dirname")]
    if not prefixes:
        return dirs
    return [d for d in dirs if any((d + "/").startswith(p) for p in prefixes)]


def dir_label(con, d):
    ver, call = con.execute("SELECT ver, call FROM data WHERE dirname=? AND ver IS NOT NULL "
                            "GROUP BY ver, call ORDER BY count(*) DESC LIMIT 1", (d,)).fetchone() or ("", "")
    return " ".join(x for x in (short(d), (ver or "")[:8], call) if x)


def median(vals):
    vals = sorted(v for v in vals if v is not None)
    if not vals:
        return None
    n = len(vals)
    m = vals[n // 2] if n % 2 else (vals[n // 2 - 1] + vals[n // 2]) / 2
    return int(m) if m == int(m) else m


def print_summary(con, dirs, fl):
    headers = ["dirname", "ver", "call", "solved", "attempted", "sat", "unsat", "PAR2",
               "avg-T-solved", "avg memMB", "memout", "sigALRM", "sigABRT", "errored",
               "med arjT", "med appmcT", "med rounds", "med sampl"]
    rows = []
    for d in dirs:
        r = con.execute(
            f"SELECT max(ver), max(call), sum(solved), count(*), sum(solved AND NOT unsat), sum(unsat),"
            f" sum(coalesce(time, {2 * TIMEOUT}))/count(*), avg(time), avg(mem_MB),"
            f" sum(mem_out), sum(signal=14), sum(signal=6), sum(errored)"
            f" FROM data WHERE dirname=?{fl}", (d,)).fetchone()
        if not r[3]:
            continue
        per = con.execute(f"SELECT arjun_time, appmc_time, rounds, sampl_sz FROM data"
                          f" WHERE dirname=? AND solved{fl}", (d,)).fetchall()
        rows.append((short(d), r[0], r[1], r[2], r[3], r[4], r[5], round(r[6]), ft(r[7]),
                     round(r[8] or 0), r[9], r[10], r[11], r[12],
                     ft(median(p[0] for p in per), 2), ft(median(p[1] for p in per), 2),
                     median(p[2] for p in per), median(p[3] for p in per)))
    rows.sort(key=lambda r: r[7], reverse=True)
    print(f"\n{BLUE}PAR2 = sum(time, or {2 * TIMEOUT}s if unsolved) / attempted. "
          f"Medians over solved instances.{RESET}")
    print_table(headers, rows)
    if len({r[4] for r in rows}) > 1:
        print(f"{RED}WARNING: 'attempted' differs between dirs, PAR2 not comparable{RESET}")


def cdf_plot(con, dirs, fl):
    gp_file, pdf_file, png_file = (f"{TMP_DIR}/cdf.{e}" for e in ("gnuplot", "pdf", "png"))
    series = []
    for d in dirs:
        times = sorted(r[0] for r in con.execute(
            f"SELECT time FROM data WHERE dirname=? AND time IS NOT NULL{fl}", (d,)))
        dat = f"{TMP_DIR}/cdf-{d}.dat"
        with open(dat, "w") as f:
            for i, t in enumerate(times):
                f.write(f"{t}\t{i + 1}\n")
        series.append(f'"{dat}" u 1:2 with lines lw 2 title "{gp_str(dir_label(con, d))} ({len(times)})"')
    with open(gp_file, "w") as f:
        for term, out in [('pdfcairo size 25cm,18cm', pdf_file), ('pngcairo size 800,600', png_file)]:
            f.write(f'set terminal {term}\nset output "{out}"\n')
            f.write('set title "ApproxMC"\nset key bottom right\nset logscale x\nset grid\n')
            f.write('set xlabel "Time (s)"\nset ylabel "Instances counted"\n')
            f.write(f'plot [1:{TIMEOUT}][:] ' + ", \\\n  ".join(series) + "\n\n")
    print(f"\n{BLUE}CDF: instances counted vs. solve time{RESET}")
    run_gnuplot(gp_file, pdf_file, png_file)


def fetch(con, d, fl):
    return {r[0]: r for r in con.execute(
        f"SELECT fname, time, mc_log10, unsat, sampl_sz, orig_sampl_sz, signal, mem_out, errored, timeout_t, arjun_time, appmc_time"
        f" FROM data WHERE dirname=?{fl}", (d,))}


def status(r):
    if r is None:
        return "missing"
    if r[1] is not None:
        return f"{r[1]:.1f}"
    if r[7]:
        return "memout"
    if r[8]:
        return "error"
    if r[6] == 14:
        return "TO"
    return f"sig{r[6]}" if r[6] else "unsolved"


def compare_pair(con, d1, d2, fl, top):
    a, b = fetch(con, d1, fl), fetch(con, d2, fl)
    common = sorted(set(a) & set(b))
    s1, s2 = short(d1), short(d2)
    print(f"\n{BLUE}==== d1 = {dir_label(con, d1)}   vs   d2 = {dir_label(con, d2)} ===={RESET}")
    if len(common) != len(a) or len(common) != len(b):
        print(f"{RED}Only {len(common)} common instances (d1 has {len(a)}, d2 has {len(b)}){RESET}")

    gained = [f for f in common if b[f][1] is not None and a[f][1] is None]
    lost = [f for f in common if a[f][1] is not None and b[f][1] is None]
    both = [f for f in common if a[f][1] is not None and b[f][1] is not None]

    def par2(x):
        return sum(x[f][1] if x[f][1] is not None else 2 * TIMEOUT for f in common) / max(len(common), 1)

    both1 = [f for f in both if max(a[f][1], b[f][1]) >= 1]
    faster = sum(1 for f in both1 if b[f][1] < a[f][1] * 0.9)
    slower = sum(1 for f in both1 if b[f][1] > a[f][1] * 1.1)
    print(f"  solved d1: {len(both) + len(lost)}  d2: {len(both) + len(gained)}   "
          f"{GREEN}gained: {len(gained)}{RESET}  {RED}lost: {len(lost)}{RESET}")
    print(f"  PAR2 d1: {par2(a):.0f}  d2: {par2(b):.0f}")
    def tsum(x, i):
        return sum(x[f][i] or 0 for f in both)

    print(f"  both solved, arjun T d1: {tsum(a, 10):.0f}s  d2: {tsum(b, 10):.0f}s;"
          f" appmc T d1: {tsum(a, 11):.0f}s  d2: {tsum(b, 11):.0f}s")
    print(f"  both solved: {len(both)}, total time d1: {sum(a[f][1] for f in both):.0f}s"
          f"  d2: {sum(b[f][1] for f in both):.0f}s; of those >=1s ({len(both1)}), d2 >10% faster on {faster}, >10% slower on {slower}")

    headers = ["fname", f"T {s1}", f"T {s2}", f"sampl {s1}", f"sampl {s2}", f"orig sampl"]

    def row(f):
        return (f, status(a[f]), status(b[f]), a[f][4], b[f][4], a[f][5] or b[f][5])

    print(f"\n{GREEN}Gained: solved by d2 ({s2}), not by d1 ({s1})  [{len(gained)}]{RESET}")
    print_table(headers, [row(f) for f in sorted(gained, key=lambda f: b[f][1])])
    print(f"\n{RED}Lost: solved by d1 ({s1}), not by d2 ({s2})  [{len(lost)}]{RESET}")
    print_table(headers, [row(f) for f in sorted(lost, key=lambda f: a[f][1])])

    big = [f for f in both if max(a[f][1], b[f][1]) >= 10]
    by_ratio = sorted(big, key=lambda f: b[f][1] / max(a[f][1], 0.01))
    hdr = headers[:3] + ["d2/d1"] + headers[3:]

    def rrow(f):
        r = row(f)
        return r[:3] + (f"{b[f][1] / max(a[f][1], 0.01):.2f}",) + r[3:]

    print(f"\n{GREEN}Top {top} speedups in d2 (both solved, one of them >=10s){RESET}")
    print_table(hdr, [rrow(f) for f in by_ratio[:top] if b[f][1] < a[f][1]])
    print(f"\n{RED}Top {top} slowdowns in d2 (both solved, one of them >=10s){RESET}")
    print_table(hdr, [rrow(f) for f in reversed(by_ratio[-top:]) if b[f][1] > a[f][1]])

    mism = []
    for f in both:
        if a[f][3] != b[f][3]:
            mism.append((f, "UNSAT" if a[f][3] else ft(a[f][2], 3), "UNSAT" if b[f][3] else ft(b[f][2], 3), "SAT/UNSAT"))
        elif a[f][2] is not None and b[f][2] is not None and abs(a[f][2] - b[f][2]) > MC_TOL:
            mism.append((f, ft(a[f][2], 3), ft(b[f][2], 3), f"{abs(a[f][2] - b[f][2]):.3f}"))
    print(f"\n{RED if mism else GREEN}Count mismatches (|log10 diff| > {MC_TOL:.3f}, i.e. >1.8^2 apart)"
          f"  [{len(mism)}]{RESET}")
    print_table(["fname", f"log10 mc {s1}", f"log10 mc {s2}", "diff"], mism)


def scatter_plot(con, d1, d2, fl):
    a, b = fetch(con, d1, fl), fetch(con, d2, fl)
    safe = re.sub(r"[^a-zA-Z0-9_-]", "_", f"{short(d1)}_vs_{short(d2)}")
    dat, gp_file, pdf_file, png_file = (f"{TMP_DIR}/scatter_{safe}.{e}" for e in ("dat", "gnuplot", "pdf", "png"))
    with open(dat, "w") as f:
        for fn in set(a) & set(b):
            f.write(f"{a[fn][1] or TIMEOUT}\t{b[fn][1] or TIMEOUT}\n")
    with open(gp_file, "w") as f:
        for term, out in [('pdfcairo size 15cm,15cm', pdf_file), ('pngcairo size 600,600', png_file)]:
            f.write(f'set terminal {term}\nset output "{out}"\n')
            f.write(f'set title "Solve time (unsolved = {TIMEOUT}s)"\n')
            f.write(f'set xlabel "{gp_str(dir_label(con, d1))} (s)"\nset ylabel "{gp_str(dir_label(con, d2))} (s)"\n')
            f.write(f'set logscale xy\nset xrange [0.1:{TIMEOUT * 1.2}]\nset yrange [0.1:{TIMEOUT * 1.2}]\n')
            f.write('set grid\nset key off\n')
            f.write(f'set arrow 1 from 0.1,0.1 to {TIMEOUT},{TIMEOUT} nohead lc rgb "gray50"\n')
            f.write(f'plot "{dat}" using 1:2 with points pt 7 ps 0.4 lc rgb "blue"\nunset arrow 1\n\n')
    print(f"\n{BLUE}Scatter: {short(d1)} vs {short(d2)}{RESET}")
    run_gnuplot(gp_file, pdf_file, png_file)


def print_failures(con, dirs, fl):
    rows = con.execute(
        f"SELECT dirname, fname, signal, mem_out, errored, exit_status, round(timeout_t,1), round(mem_MB)"
        f" FROM data WHERE dirname IN ({','.join('?' * len(dirs))}) AND NOT solved"
        f" AND (coalesce(signal,0) != 14 OR errored OR mem_out){fl} ORDER BY fname, dirname", dirs).fetchall()
    print(f"\n{BLUE}Unsolved for reasons other than timeout  [{len(rows)}]{RESET}")
    print_table(["dirname", "fname", "signal", "memout", "errored", "exit", "T", "memMB"],
                [(short(r[0]),) + tuple(r[1:]) for r in rows])


def main():
    global TIMEOUT
    p = argparse.ArgumentParser(description="Tables and plots for approxmc benchmark data")
    p.add_argument("--dirs", nargs="+", metavar="PREFIX", default=[], help="Only these dirs (prefix match)")
    p.add_argument("--fname", nargs="+", metavar="PATTERN", default=[],
                   help="Filter by SQL LIKE pattern(s), e.g. --fname '%%track1%%'")
    p.add_argument("--pair", nargs="+", metavar="DIR", default=[],
                   help="Compare only these dirs, in this order (first = baseline); a single dir is the baseline for all others")
    p.add_argument("--cdf", action="store_true", help="ONLY the summary table and the CDF plot")
    p.add_argument("--nopairwise", action="store_true")
    p.add_argument("--top", type=int, default=15, help="Rows in speedup/slowdown tables")
    p.add_argument("--timeout", type=int, default=TIMEOUT)
    args = p.parse_args()
    TIMEOUT = args.timeout

    os.makedirs(TMP_DIR, exist_ok=True)
    fl = ""
    if args.fname:
        fl = " AND (" + " OR ".join("fname LIKE '" + x.replace("'", "''") + "'" for x in args.fname) + ")"

    con = sqlite3.connect(DB)
    dirs = get_dirs(con, args.dirs)
    if not dirs:
        print("No matching dirs in the DB")
        return

    print_header("solved counts and PAR2")
    print_summary(con, dirs, fl)
    if not args.cdf:
        print_failures(con, dirs, fl)
    if not args.cdf and not args.nopairwise:
        # d1 is the baseline, d2 the new one: gained/lost is from d2's point of view
        pairs = list(itertools.combinations(sorted(dirs, reverse=True), 2))
        if args.pair:
            sel = [d for x in args.pair for d in dirs if (d + "/").startswith(x)]
            pairs = [(sel[0], d) for d in dirs if d != sel[0]] if len(sel) == 1 else list(itertools.combinations(sel, 2))
        print_header("pairwise comparisons")
        for d1, d2 in pairs:
            compare_pair(con, d1, d2, fl, args.top)
            scatter_plot(con, d1, d2, fl)
    cdf_plot(con, dirs, fl)


if __name__ == "__main__":
    main()
