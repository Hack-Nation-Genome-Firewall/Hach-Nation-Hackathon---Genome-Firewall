"""
Parallel AMRFinderPlus over a cohort, streaming FASTAs (download -> run -> delete).
Each worker keeps only the small TSV. Resumes: skips genomes whose TSV exists.

  python run_cohort_parallel.py --ids cohort_target_ids.json --db <db> \
         --amrfinder <path> --outdir cohort_tsv --workers 6 --amr-threads 2
"""
import argparse, json, os, subprocess, time, urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed

def download_fasta(gid, outdir, timeout=180):
    url = (f"https://www.bv-brc.org/api/genome_sequence/?eq(genome_id,{gid})"
           f"&sort(+sequence_id)&limit(10000)&http_accept=application/dna+fasta")
    req = urllib.request.Request(url, headers={"Accept": "application/dna+fasta"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        data = r.read()
    p = os.path.join(outdir, f"{gid}.fasta")
    open(p, "wb").write(data)
    return p

def process(gid, amrfinder, db, outdir, fasta_dir, threads, organism="Klebsiella_pneumoniae"):
    tsv = os.path.join(outdir, f"{gid}.tsv")
    if os.path.exists(tsv) and os.path.getsize(tsv) > 0:
        return gid, "skip"
    fa = None
    try:
        fa = download_fasta(gid, fasta_dir)
        subprocess.run([amrfinder, "-n", fa, "--organism", organism, "--plus",
                        "--database", db, "-o", tsv, "--threads", str(threads)],
                       check=True, capture_output=True, text=True)
        return gid, "ok"
    except Exception as e:
        return gid, f"ERR {type(e).__name__} {str(e)[:80]}"
    finally:
        if fa and os.path.exists(fa):
            os.remove(fa)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ids", required=True)
    ap.add_argument("--db", required=True)
    ap.add_argument("--amrfinder", required=True)
    ap.add_argument("--outdir", default="cohort_tsv")
    ap.add_argument("--fasta-dir", default="fasta_tmp")
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--amr-threads", type=int, default=2)
    args = ap.parse_args()
    if args.ids.endswith(".csv"):
        import csv
        with open(args.ids) as fh:
            ids = [r["genome_id"] for r in csv.DictReader(fh)]
    else:
        ids = json.load(open(args.ids))
    os.makedirs(args.outdir, exist_ok=True)
    os.makedirs(args.fasta_dir, exist_ok=True)
    todo = [g for g in ids if not (os.path.exists(os.path.join(args.outdir, f"{g}.tsv"))
                                   and os.path.getsize(os.path.join(args.outdir, f"{g}.tsv")) > 0)]
    print(f"{len(ids)} total, {len(ids)-len(todo)} already done, {len(todo)} to run", flush=True)
    t0 = time.time(); n_ok = n_err = 0
    with ThreadPoolExecutor(max_workers=args.workers) as ex:
        futs = {ex.submit(process, g, args.amrfinder, args.db, args.outdir,
                          args.fasta_dir, args.amr_threads): g for g in todo}
        for i, fut in enumerate(as_completed(futs), 1):
            gid, status = fut.result()
            if status == "ok": n_ok += 1
            elif status.startswith("ERR"): n_err += 1; print(f"  {gid}: {status}", flush=True)
            if i % 25 == 0:
                el = time.time() - t0
                print(f"  {i}/{len(todo)}  ok={n_ok} err={n_err}  {el:.0f}s  "
                      f"({el/i:.1f}s/genome, eta {el/i*(len(todo)-i)/60:.0f}m)", flush=True)
    print(f"DONE ok={n_ok} err={n_err} in {time.time()-t0:.0f}s", flush=True)

if __name__ == "__main__":
    main()
