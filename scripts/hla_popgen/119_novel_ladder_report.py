"""119 build the novel-allele ladder report (HTML, self-contained) from 118's ladder_data JSON.

  python3 scripts/hla_popgen/119_novel_ladder_report.py --data ~/s10/ladder_data.private.json --out ~/s10/S08_NOVEL_LADDER_EXACT.html     # on the VM; exact counts, never leaves the Workbench
  python3 scripts/hla_popgen/119_novel_ladder_report.py --data ladder_data.pub.json --out reports/sprint_reports/S08_REPORT.html            # cloud-safe (counts 1-19 are classes)
  python3 scripts/hla_popgen/119_novel_ladder_report.py --out reports/sprint_reports/S08_REPORT.html                                          # no data: static sections + a file loader
Refuses to embed a data file whose meta.masked is false unless --out ends in .private.html or --allow-exact is given, so an exact build cannot be dropped into reports/ by accident.
"""
import argparse
import datetime
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))


def build(data, release="3.66.0"):
    tpl = open(os.path.join(HERE, "119_report_template.html"), encoding="utf-8").read()
    rel = (data or {}).get("meta", {}).get("ipd_release", release)
    return (tpl.replace("__BUILT__", datetime.date.today().isoformat()).replace("__REL__", rel)
               .replace("__DATA_JSON__", json.dumps(data, separators=(",", ":")) if data else "null"))


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--data")
    p.add_argument("--out", required=True)
    p.add_argument("--allow-exact", action="store_true")
    a = p.parse_args(argv)
    data = json.load(open(os.path.expanduser(a.data))) if a.data else None
    if data and not data["meta"].get("masked", True) and not (a.allow_exact or a.out.endswith(".private.html")):
        sys.exit("refusing: exact-count data; use --out *.private.html (stays on the VM) or --allow-exact")
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    open(a.out, "w", encoding="utf-8").write(build(data))
    print("wrote", a.out)


if __name__ == "__main__":
    main()
