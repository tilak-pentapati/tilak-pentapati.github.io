"""Build assets/data/projects.json for the home-page map and skills network.

Single source of truth: each project's index.qmd front matter (title, date,
categories) plus its folder's _metadata.yml (summary, study_area), and
projects/external.yml for projects hosted elsewhere. Adding a project only
needs that metadata; re-run this script and commit the JSON:

    /Users/jay/anaconda3/envs/envs563/bin/python scripts/build_projects_json.py

Needs only PyYAML and no raw data, so it is safe to run anywhere.
"""

import json
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "assets" / "data" / "projects.json"


def front_matter(qmd):
    text = qmd.read_text(encoding="utf-8")
    if not text.startswith("---"):
        return {}
    return yaml.safe_load(text.split("---", 2)[1]) or {}


def read_yaml(path):
    if not path.exists():
        return {}
    return yaml.safe_load(path.read_text(encoding="utf-8")) or {}


def site_projects():
    for qmd in sorted(ROOT.glob("projects/*/*/index.qmd")):
        fm = front_matter(qmd)
        meta = read_yaml(qmd.parent / "_metadata.yml")
        summary = meta.get("summary", {})
        yield {
            "id": qmd.parent.name,
            "title": fm["title"],
            "short_title": meta.get("short_title", fm["title"]),
            "href": qmd.relative_to(ROOT).with_suffix(".html").as_posix(),
            "external": False,
            "date": str(fm.get("date", "")),
            "blurb": summary.get("problem") or fm.get("description", ""),
            "categories": fm.get("categories", []),
            "stack": summary.get("stack", []),
            "study_area": meta.get("study_area"),
        }


def external_projects():
    for item in read_yaml(ROOT / "projects" / "external.yml") or []:
        yield {
            "id": item["path"].rstrip("/").rsplit("/", 1)[-1],
            "title": item["title"],
            "short_title": item.get("short_title", item["title"]),
            "href": item["path"],
            "external": True,
            "date": str(item.get("date", "")),
            "blurb": item.get("description", ""),
            "categories": item.get("categories", []),
            "stack": item.get("stack", []),
            "study_area": item.get("study_area"),
        }


def main():
    projects = list(site_projects()) + list(external_projects())
    for p in projects:
        if not p["study_area"]:
            raise SystemExit(f"{p['id']}: missing study_area in its metadata")
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps({"projects": projects}, indent=1, ensure_ascii=False) + "\n",
                   encoding="utf-8")
    print(f"wrote {OUT.relative_to(ROOT)} ({len(projects)} projects)")


if __name__ == "__main__":
    main()
