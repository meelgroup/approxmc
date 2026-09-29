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


def print_table(headers, rows, wrap=None):
    """wrap = {column index: width} hard-wraps those cells onto extra lines"""
    wrap = wrap or {}
    rows = [["" if v is None else str(v) for v in r] for r in rows]
    lines = []
    for r in rows:
        cells = [[c[j:j + wrap[i]] for j in range(0, len(c), wrap[i])] or [""] if i in wrap else [c]
                 for i, c in enumerate(r)]
        lines += [[c[k] if k < len(c) else "" for c in cells] for k in range(max(len(c) for c in cells))]
    widths = [max(len(h), max((len(r[i]) for r in lines), default=0)) for i, h in enumerate(headers)]
    sep = "+-" + "-+-".join("-" * w for w in widths) + "-+"
    fmt = "| " + " | ".join(f"{{:<{w}}}" for w in widths) + " |"
    print(sep)
    print(fmt.format(*headers))
    print(sep)
    for r in lines:
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


# validated colorblind-safe categorical order
COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#4a3aa7"]
GREY = "#b9bec8"


def gp_rgb(hexcol):
    return int(hexcol[1:], 16)


def gp_output(f, pdf_size, png_size, pdf_file, png_file):
    # yields once per terminal so the same plot commands are written for both
    for term, out in [(f"pdfcairo size {pdf_size}", pdf_file), (f"pngcairo size {png_size}", png_file)]:
        f.write(f'set terminal {term}\nset output "{out}"\n')
        yield


def load_steps(con, d, fl):
    steps = {}
    for r in con.execute(f"SELECT fname, round, hashes, dur, finished, above FROM steps"
                         f" WHERE dirname=?{fl} ORDER BY fname, idx", (d,)):
        steps.setdefault(r[0], []).append(r[1:])
    return steps


def round_answers(st):
    """Hash count each finished round settled on: the lowest count with fewer than threshold+1 solutions."""
    ans = {}
    for i, (rnd, h, _, fin, above) in enumerate(st):
        if fin and (i + 1 == len(st) or st[i + 1][0] != rnd):
            ans[rnd] = h if above == 0 else h + 1
    return ans


def stuck_info(st):
    rnd, cur = st[-1][0], st[-1][1]
    rs = [x for x in st if x[0] == rnd]
    above = [x[1] for x in rs if x[4] == 1]
    below = [x[1] for x in rs if x[4] == 0]
    lo = max(above) + 1 if above else None
    hi = min(below) if below else None
    phase = "recount" if rnd > 0 else ("gallop" if hi is None else "binary search")
    return rnd, cur, lo, hi, phase


def time_breakdown(con, d, fl, steps):
    buckets = ["Arjun preprocessing", "Arjun minimization", "Round 0", "Round 1", "Round 2+", "Other"]
    res = {"all": dict.fromkeys(buckets, 0.0), "solved": dict.fromkeys(buckets, 0.0)}
    for fname, solved, tot, arj, bw in con.execute(
            f"SELECT fname, solved, timeout_t, arjun_time, backward_t FROM data WHERE dirname=?{fl}", (d,)):
        b = dict.fromkeys(buckets, 0.0)
        if arj is None:
            b["Other"] = tot
        else:
            b["Arjun preprocessing"] = arj - (bw or 0)
            b["Arjun minimization"] = bw or 0
            used = arj
            for rnd, _, dur, _, _ in steps.get(fname, []):
                b[f"Round {rnd}" if rnd < 2 else "Round 2+"] += dur or 0
                used += dur or 0
            b["Other"] = max(0.0, tot - used)
        for k in ("all", "solved") if solved else ("all",):
            for x in buckets:
                res[k][x] += b[x]
    return buckets, res


def pie_plot(d, buckets, res):
    for k, title in (("all", "All instances, timeouts at full runtime"), ("solved", "Solved instances only")):
        gp_file, pdf_file, png_file = (f"{TMP_DIR}/pie-{k}-{short(d)}.{e}" for e in ("gnuplot", "pdf", "png"))
        total = sum(res[k].values())
        with open(gp_file, "w") as f:
            for _ in gp_output(f, "18cm,10cm", "700,390", pdf_file, png_file):
                f.write(f'set title "{title}: {total / 3600:.1f} CPU hours"\nset size ratio -1\nunset border\n'
                        'unset tics\nset xrange [-1.05:3.1]\nset yrange [-1.05:1.05]\n'
                        "set style fill solid 1 border lc rgb 'white'\n"
                        "set key at 1.25,0.8 left top Left reverse samplen 1 spacing 1.4\n")
                parts, a = [], 90.0
                for i, b in enumerate(buckets):
                    ang = 360 * res[k][b] / total if total else 0
                    parts.append(f"'+' u (0):(0):(1):({a - ang}):({a}) every ::0::0 w circles lc rgb '{COLORS[i]}'"
                                 f" title '{b} ({100 * res[k][b] / total if total else 0:.1f}%)'")
                    a -= ang
                f.write("plot " + ", \\\n  ".join(parts) + "\n\n")
        print(f"\n{BLUE}Pie: CPU time per phase, {title.lower()}{RESET}")
        run_gnuplot(gp_file, pdf_file, png_file)


def timeline_plot(d, rows):
    safe = short(d)
    dat, gp_file, pdf_file, png_file = (f"{TMP_DIR}/timeline-{safe}.{e}" for e in ("dat", "gnuplot", "pdf", "png"))
    ytics = []
    with open(dat, "w") as f:
        for y, r in enumerate(rows):
            ytics.append(f'"{gp_str(r["fname"].replace(".cnf", "")[:30])}" {y}')
            f.write(f"0\t{r['arjT']}\t{y}\t{gp_rgb(GREY)}\n")
            t = r["arjT"]
            for i, (_, _, dur, _, _) in enumerate(r["st"]):
                last = i + 1 == len(r["st"])
                f.write(f"{t}\t{t + (dur or 0)}\t{y}\t{gp_rgb(COLORS[1] if last else COLORS[0])}\n")
                t += dur or 0
    h = 120 + 16 * len(rows)
    with open(gp_file, "w") as f:
        for _ in gp_output(f, f"25cm,{h / 38:.1f}cm", f"1000,{h}", pdf_file, png_file):
            f.write(f'set title "Timeouts: time spent per bounded\\\\_sol\\\\_count call"\n'
                    f'set xrange [0:{TIMEOUT}]\nset yrange [-0.7:{len(rows) - 0.3}] reverse\n'
                    f'set xlabel "CPU time (s)"\nset ytics ({", ".join(ytics)}) font ",8"\nset grid x\n'
                    'set style fill solid 1 border lc rgb "white"\nset key outside top center horizontal\n')
            f.write(f'plot "{dat}" u (($1+$2)/2):3:1:2:($3-0.4):($3+0.4):4 w boxxyerror lc rgb variable notitle, \\\n'
                    f'  NaN w boxes lc rgb "{GREY}" title "Arjun", \\\n'
                    f'  NaN w boxes lc rgb "{COLORS[0]}" title "finished call", \\\n'
                    f'  NaN w boxes lc rgb "{COLORS[1]}" title "call running at timeout"\n\n')
    print(f"\n{BLUE}Timeline of timeouts: each bar is one timed-out instance from start to the timeout, split into"
          f" Arjun (grey), each bounded_sol_count call that finished (blue) and the call still running when"
          f" the timeout hit (orange), sorted by how much of the run that last call took.{RESET}")
    run_gnuplot(gp_file, pdf_file, png_file)


def dist_plot(d, rows):
    safe = short(d)
    dat, gp_file, pdf_file, png_file = (f"{TMP_DIR}/dist-{safe}.{e}" for e in ("dat", "gnuplot", "pdf", "png"))
    with open(dat, "w") as f:
        for label, n, tot, *_ in rows:
            col = COLORS[0] if label.startswith(("-", "<=")) else COLORS[1]
            f.write(f'"{label}"\t{tot / 3600}\t{gp_rgb(col)}\n')
    with open(gp_file, "w") as f:
        for _ in gp_output(f, "20cm,12cm", "800,480", pdf_file, png_file):
            f.write('set title "Solved instances: call time by distance to the round\'s final hash count"\n'
                    'set xlabel "hash count of call - final hash count"\nset ylabel "CPU hours"\n'
                    'set style fill solid 1 border lc rgb "white"\nset boxwidth 0.8\nset grid y\n'
                    'set yrange [0:*]\nset key top left\n')
            f.write(f'plot "{dat}" u 0:2:3:xtic(1) w boxes lc rgb variable notitle, \\\n'
                    f'  NaN w boxes lc rgb "{COLORS[0]}" title ">= threshold+1 solutions", \\\n'
                    f'  NaN w boxes lc rgb "{COLORS[1]}" title "< threshold+1 (cell exhausted)"\n\n')
    print(f"\n{BLUE}Call time vs. distance to the final hash count{RESET}")
    run_gnuplot(gp_file, pdf_file, png_file)


def hashes_plot(d, solved_pts, unsolved):
    safe = short(d)
    dat, dat2, gp_file, pdf_file, png_file = (f"{TMP_DIR}/hashes-{safe}.{e}" for e in ("dat", "uns.dat", "gnuplot", "pdf", "png"))
    with open(dat, "w") as f:
        for sz, a in solved_pts:
            f.write(f"{sz}\t{a}\n")
    with open(dat2, "w") as f:
        for r in unsolved:
            lo, hi = r["lo"] or r["cur"], r["hi"] or r["cur"]
            f.write(f"{r['sz']}\t{lo}\t{hi}\n")
    with open(gp_file, "w") as f:
        for _ in gp_output(f, "18cm,14cm", "700,560", pdf_file, png_file):
            f.write('set title "Hash count needed (round 0) vs. sampling set size"\nset logscale xy\nset grid\n'
                    'set xlabel "sampling set size after Arjun"\nset ylabel "hash count"\nset key top left\n')
            f.write(f'plot "{dat}" u 1:2 w points pt 7 ps 0.4 lc rgb "{COLORS[0]}" title "solved", \\\n'
                    f'  "{dat2}" u 1:2:(0):($3-$2) w vectors nohead lw 3 lc rgb "{COLORS[1]}" title "unsolved: from lo to hi (or stuck-at)", \\\n'
                    f'  x w lines dt 2 lc rgb "gray50" title "hashes = sampling set size"\n\n')
    print(f"\n{BLUE}Hash count needed vs. sampling set size{RESET}")
    run_gnuplot(gp_file, pdf_file, png_file)


def calls_plot(d, hist):
    safe = short(d)
    dat, gp_file, pdf_file, png_file = (f"{TMP_DIR}/calls-{safe}.{e}" for e in ("dat", "gnuplot", "pdf", "png"))
    with open(dat, "w") as f:
        for k in range(1, max(hist) + 1):
            f.write(f"{k}\t{hist.get(k, 0)}\n")
    with open(gp_file, "w") as f:
        for _ in gp_output(f, "18cm,11cm", "700,430", pdf_file, png_file):
            f.write('set title "bounded\\\\_sol\\\\_count calls in round 0 (solved)"\nset xlabel "calls"\n'
                    'set ylabel "instances"\nset style fill solid 1 border lc rgb "white"\nset boxwidth 0.85\n'
                    'set grid y\nset yrange [0:*]\nset xrange [0.3:*]\nset xtics 1\nunset key\n')
            f.write(f'plot "{dat}" u 1:2 w boxes lc rgb "{COLORS[0]}"\n\n')
    print(f"\n{BLUE}Calls per round 0{RESET}")
    run_gnuplot(gp_file, pdf_file, png_file)


def boundary_table(con, d, fl):
    if not con.execute("SELECT name FROM sqlite_master WHERE name='rounds'").fetchone():
        return
    fl_r = fl.replace("fname", "r.fname")
    groups = {"round 0": "r.round = 0", "rounds 1+": "r.round > 0", "all": "1"}
    rows = []
    for label, cond in groups.items():
        n, bs, bu, os_, ou = con.execute(
            f"SELECT COUNT(DISTINCT r.fname), SUM(r.bnd_sat_t), SUM(r.bnd_unsat_t), SUM(r.oth_sat_t), SUM(r.oth_unsat_t)"
            f" FROM rounds r JOIN data USING (dirname, fname) WHERE r.dirname=? AND data.solved=1 AND {cond}{fl_r}",
            (d,)).fetchone()
        if not n:
            continue
        tot = bs + bu + os_ + ou
        rows.append((label, n, ft(bs / 3600, 2), ft(bu / 3600, 2), ft(os_ / 3600, 2), ft(ou / 3600, 2),
                     ft(100 * (bs + bu) / tot if tot else 0), ft(100 * (bu + ou) / tot if tot else 0)))
    print(f"\n{BLUE}Solved instances: solve() time in each round, boundary (the two hash counts the round ended on)"
          f" vs. the rest, split into finding solutions and the final 'no more solutions' call.{RESET}")
    if not rows:
        print("  No per-round times in this dir: needs an approxmc that prints 'round N solve T:' at --verb 1")
        return
    print_table(["rounds", "inst", "bnd sols CPU h", "bnd UNSAT CPU h", "other sols CPU h", "other UNSAT CPU h",
                 "bnd %", "UNSAT %"], rows)


def analyze(con, d, fl, graphs):
    if not con.execute("SELECT name FROM sqlite_master WHERE name='steps'").fetchone():
        print(f"{RED}No 'steps' table, rerun get_data_approxmc.py --reparse{RESET}")
        return
    print_header(f"analysis of {short(d)}")
    steps = load_steps(con, d, fl)
    data = {r[0]: r for r in con.execute(
        f"SELECT fname, solved, signal, mem_out, timeout_t, arjun_time, appmc_time, sampl_sz, orig_sampl_sz, mem_MB"
        f" FROM data WHERE dirname=?{fl}", (d,))}
    others = {}
    for fn, dn, sv in con.execute(f"SELECT fname, dirname, solved FROM data WHERE dirname!=?{fl}", (d,)):
        others.setdefault(fn, []).append((dn, sv))

    # unsolved
    unsolved = []
    for fn, r in data.items():
        if r[1]:
            continue
        u = {"fname": fn, "sz": r[7], "orig": r[8], "arjT": r[5], "t": r[4], "mem": r[9], "memout": r[3], "st": steps.get(fn)}
        if u["st"] and r[5] is not None:
            u["rnd"], u["cur"], u["lo"], u["hi"], u["phase"] = stuck_info(u["st"])
            u["stuck"] = u["st"][-1][2] or 0
            u["share"] = u["stuck"] / max(r[4] - r[5], 1e-9)
        unsolved.append(u)
    timeouts = [u for u in unsolved if "phase" in u]
    other_runs = sorted({dn for v in others.values() for dn, _ in v})

    def others_col(fn):
        return " ".join("Y" if sv else "-" for _, sv in sorted(others.get(fn, [])))

    print(f"\n{BLUE}Unsolved instances. bracket = [lo, hi]: lo-1 hashes still gave >= threshold+1 solutions,"
          f" hi gave fewer (in the round being worked on). share = stuck call / ApproxMC time.{RESET}")
    if other_runs:
        print(f"{BLUE}solved-elsewhere column order: {' '.join(short(x) for x in other_runs)}{RESET}")
    rows = []
    for u in sorted(unsolved, key=lambda u: (-u.get("share", 2), u["fname"])):
        if "phase" in u:
            rows.append((u["fname"], u["sz"], u["orig"], ft(u["arjT"]), len(u["st"]), u["rnd"], u["cur"],
                         f"[{'?' if u['lo'] is None else u['lo']}, {'?' if u['hi'] is None else u['hi']}]",
                         ft(u["stuck"], 0), f"{100 * u['share']:.0f}%", u["phase"], others_col(u["fname"])))
        else:
            why = "memout" if u["memout"] else "no appmc output"
            rows.append((u["fname"], u["sz"], u["orig"], ft(u["arjT"]), "", "", "", "", "", "",
                         f"{why} ({ft(u['t'], 0)}s, {ft((u['mem'] or 0) / 1000)}GB)", others_col(u["fname"])))
    print_table(["fname", "sampl", "orig", "arjT", "calls", "round", "stuck at", "bracket", "stuck T", "share",
                 "phase", "solved elsewhere"], rows, wrap={0: 25})

    if timeouts:
        shares = sorted(u["share"] for u in timeouts)
        print(f"\n{BLUE}Timeouts by round reached and search phase  [{len(timeouts)}]{RESET}")
        cnt = {}
        for u in timeouts:
            cnt[(u["rnd"], u["phase"])] = cnt.get((u["rnd"], u["phase"]), 0) + 1
        print_table(["round", "phase", "count"], [(k[0], k[1], v) for k, v in sorted(cnt.items())])
        print(f"  stuck call share of ApproxMC time: median {100 * median(shares):.0f}%,"
              f" >90% on {sum(s > .9 for s in shares)} of {len(shares)}")

    # time breakdown
    buckets, res = time_breakdown(con, d, fl, steps)
    ta, ts = sum(res["all"].values()), sum(res["solved"].values())
    print(f"\n{BLUE}CPU time per phase. Inter-round simplification is counted in the round it follows.{RESET}")
    print_table(["phase", "all inst: CPU h", "all inst: %", "solved inst: CPU h", "solved inst: %"],
                [(b, ft(res["all"][b] / 3600, 2), ft(100 * res["all"][b] / ta), ft(res["solved"][b] / 3600, 2),
                  ft(100 * res["solved"][b] / ts)) for b in buckets] +
                [("TOTAL", ft(ta / 3600, 2), "100.0", ft(ts / 3600, 2), "100.0")])
    print("  'all inst' sums over every instance with timeouts counted at their full runtime, while"
          " 'solved inst' sums only over instances that got a count, i.e. where a successful run spends its time.")
    boundary_table(con, d, fl)

    # distance to the answer, on solved instances with non-trivial ApproxMC time
    order = ["<=-5", "-4", "-3", "-2", "-1", "0", "+1", "+2", "+3", "+4", ">=+5"]
    durs = {k: [] for k in order}
    hard = 0
    solved_pts, calls_hist, calls_later, top_share = [], {}, {}, []
    for fn, r in data.items():
        st = steps.get(fn)
        if not r[1] or not st:
            continue
        ans = round_answers(st)
        if 0 in ans and ans[0] > 0 and r[7]:
            solved_pts.append((r[7], ans[0]))
        for rnd in ans:
            n = sum(1 for x in st if x[0] == rnd)
            h = calls_hist if rnd == 0 else calls_later
            h[n] = h.get(n, 0) + 1
        if (r[6] or 0) < 10:
            continue
        hard += 1
        top_share.append(max(x[2] or 0 for x in st) / r[6])
        last_rnd = max(ans)
        for i, (rnd, h, dur, _, _) in enumerate(st):
            ends_round = i + 1 == len(st) or st[i + 1][0] != rnd
            if ends_round and rnd != last_rnd:
                continue  # includes the inter-round simplify
            dd = h - ans[rnd]
            k = "<=-5" if dd <= -5 else ">=+5" if dd >= 5 else f"{dd:+d}".replace("+0", "0")
            durs[k].append(dur or 0)
    tot = sum(sum(v) for v in durs.values())
    drows = []
    for k in order:
        v = sorted(durs[k])
        drows.append((k, len(v), sum(v), 100 * sum(v) / tot if tot else 0,
                      median(v), v[int(0.9 * (len(v) - 1))] if v else None, v[-1] if v else None))
    print(f"\n{BLUE}Solved with ApproxMC T >= 10s [{hard}]: call time by (hash count - round's final hash count).{RESET}")
    print_table(["distance", "calls", "total T", "% of T", "median T", "p90 T", "max T"],
                [(k, n, ft(t, 0), ft(p), ft(m, 2), ft(q, 2), ft(x, 1)) for k, n, t, p, m, q, x in drows])
    print("  Every bounded_sol_count call is binned by how far its hash count was from the count its round settled on:"
          " negative means the cell still had >= threshold+1 solutions, 0 or more means it had fewer and so had to be"
          " enumerated to the end, which needs a final UNSAT call.")
    print("  The time columns show where ApproxMC's time goes, and it concentrates at -1, 0 and +1,"
          " where the cell holds close to threshold solutions.")
    if top_share:
        print(f"  most expensive single call: median {100 * median(top_share):.0f}% of ApproxMC time,"
              f" >50% on {sum(s > .5 for s in top_share)} of {len(top_share)}")

    print(f"\n{BLUE}bounded_sol_count calls per round (solved){RESET}")
    ks = sorted(set(calls_hist) | set(calls_later))
    print_table(["calls"] + [str(k) for k in ks], [["round 0"] + [calls_hist.get(k, 0) for k in ks],
                                                   ["rounds 1+"] + [calls_later.get(k, 0) for k in ks]])
    fin = sorted(a for _, a in solved_pts)
    if fin:
        print(f"  round-0 final hash count over {len(fin)} solved: median {median(fin)},"
              f" p90 {fin[int(0.9 * (len(fin) - 1))]}, max {fin[-1]}")

    if not graphs:
        return
    pie_plot(d, buckets, res)
    if timeouts:
        timeline_plot(d, sorted(timeouts, key=lambda u: -u["share"]))
    dist_plot(d, drows)
    hashes_plot(d, solved_pts, timeouts)
    if calls_hist:
        calls_plot(d, calls_hist)


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
    p.add_argument("--analyze", metavar="PREFIX",
                   help="ONLY the summary table and a per-call analysis of this dir: unsolved, time per phase, hashes")
    p.add_argument("--nographs", action="store_true", help="With --analyze: tables only")
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
    if args.analyze:
        sel = [d for d in dirs if (d + "/").startswith(args.analyze) or short(d) == args.analyze]
        if len(sel) != 1:
            print(f"{RED}--analyze must match exactly one dir, matched: {sel}{RESET}")
            return
        analyze(con, sel[0], fl, not args.nographs)
        return
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
