#!/usr/bin/env python3
"""Builds the HTML manual from the Markdown sources listed in docs/manual/manifest.json.

    python3 tools/build_manual.py            writes docs/manual/html/
    python3 tools/build_manual.py --check    exits with 1 if docs/manual/html/ is not up to date (used by the CI and the tests),
                                             with 2 if the installed `markdown` is not the version that generated it

The HTML is committed, so building ligint needs neither Python nor this script: it is embedded in the GUI (Help menu) through
docs/manual/html/manual.qrc and installed as plain files. Only the people who edit the manual need `pip install -r tools/requirements-manual.txt`.

Pages are made of Markdown files, some of them shared with the rest of the documentation (docs/*.md), and of sections taken
from other files with a line of its own such as   <!-- include: README.md "What it recognizes (default thresholds)" -->   so that numbers
that appear in more than one place have a single source. The markup is kept simple (tables, lists, code blocks, images) so that
it displays well both in the Qt text browser of the GUI and in a web browser.
"""
import argparse, html, json, os, re, shutil, struct, sys, tempfile

try:
    import markdown
except ImportError:
    sys.exit("This script needs the Python module 'markdown' (pip install -r tools/requirements-manual.txt)")

MARKDOWN_VERSION = "3.10.2"      # the version that the committed HTML was generated with (tools/requirements-manual.txt)
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MANUAL = os.path.join(ROOT, "docs", "manual")
OUT = os.path.join(MANUAL, "html")

CSS = """
body { margin: 12px 18px; }
h1 { font-size: 26px; margin-top: 6px; }
h2 { font-size: 20px; margin-top: 22px; border-bottom: 1px solid #adb5bd; }
h3 { font-size: 16px; margin-top: 16px; }
h4 { font-size: 14px; }
pre { background-color: #f1f3f5; color: #212529; padding: 6px; }
code { font-family: monospace; }
th { background-color: #dee2e6; color: #212529; }
blockquote { margin-left: 18px; }
p.nav { background-color: #e9ecef; color: #212529; padding: 4px; }
p.foot { font-size: 11px; color: #868e96; margin-top: 24px; }
a { color: #1c7ed6; }
"""

def read_version():
    with open(os.path.join(ROOT, "include", "ligint", "ligint.hpp"), encoding="utf-8") as f:
        m = re.search(r'#define LIGINT_VERSION "([^"]+)"', f.read())
    return m.group(1) if m else "?"

def read(path):
    with open(os.path.join(ROOT, path), encoding="utf-8") as f:
        return f.read()

def section(text, heading):
    """The part of `text` from the heading `heading` up to the next heading of the same or a higher level."""
    lines = text.split("\n")
    start = level = None
    in_code = False
    for i, l in enumerate(lines):
        if l.startswith("```"):
            in_code = not in_code
        m = re.match(r"^(#+) (.*)$", l) if not in_code else None
        if start is None:
            if m and m.group(2).strip() == heading:
                start, level = i, len(m.group(1))
        elif m and len(m.group(1)) <= level:
            return "\n".join(lines[start:i]).rstrip() + "\n"
    if start is None:
        raise KeyError(heading)
    return "\n".join(lines[start:]).rstrip() + "\n"

def expand_includes(text, errors, where):
    """Replaces the include lines with the sections they name. The relative links and images of the included text are rewritten so that
    they are correct from the including page (they were written relative to the file the section comes from)."""
    def relocate(chunk, origin):
        def fix(m):
            target = m.group(2)
            if re.match(r"^(https?:|mailto:|#)", target):
                return m.group(0)
            path, sep, frag = target.partition("#")
            moved = os.path.relpath(os.path.normpath(os.path.join(os.path.dirname(origin), path)), os.path.dirname(where)).replace(os.sep, "/")
            return "%s(%s%s%s)" % (m.group(1), moved, sep, frag)
        return re.sub(r"(\]|\!\[[^\]]*\])\(([^)\s]+)\)", fix, chunk)
    def sub(m):
        try:
            return relocate(section(read(m.group(1)), m.group(2)), m.group(1))
        except (KeyError, OSError):
            errors.append("%s: include not found: %s \"%s\"" % (where, m.group(1), m.group(2)))
            return ""
    return re.sub(r'^<!-- include: (\S+) "(.+?)" -->[ \t]*$', sub, text, flags=re.M)         # only a line of its own

def exists_exact(path):
    """os.path.exists with the case of every component checked, so that Windows and macOS give the same answer as Linux."""
    full = os.path.join(ROOT, path)
    if not os.path.exists(full):
        return False
    rel = os.path.relpath(full, ROOT)
    parent = ROOT
    for part in rel.split(os.sep):
        if part not in (os.curdir, os.pardir) and part not in os.listdir(parent):
            return False
        parent = os.path.join(parent, part)
    return True

def png_size(path):
    with open(path, "rb") as f:
        head = f.read(24)
    if head[:8] == b"\x89PNG\r\n\x1a\n":
        return struct.unpack(">II", head[16:24])
    return None

class Builder:
    def __init__(self, manifest):
        self.m = manifest
        self.pages = manifest["pages"]
        self.sources = {os.path.normpath(p["src"]): p for p in self.pages}
        self.ids = {p["id"] for p in self.pages}
        self.errors = []
        self.images = {}          # output name -> source path
        self.version = read_version()

    # ---- links and images ------------------------------------------------------------------------------------------
    def rewrite_href(self, href, src):
        if re.match(r"^(https?:|mailto:|#)", href):
            return href
        path, _, frag = href.partition("#")
        target = os.path.normpath(os.path.join(os.path.dirname(src), path))
        base = os.path.basename(path)
        if os.path.dirname(path) == "" and base.endswith((".md", ".html")) and base.rsplit(".", 1)[0] in self.ids and not exists_exact(target):
            return base.rsplit(".", 1)[0] + ".html" + ("#" + frag if frag else "")     # a link to a page of the manual, by its id
        if target in self.sources:
            return self.sources[target]["id"] + ".html" + ("#" + frag if frag else "")
        if os.path.exists(os.path.join(ROOT, target)):
            kind = "tree" if os.path.isdir(os.path.join(ROOT, target)) else "blob"
            return "%s/%s/main/%s" % (self.m["repo_url"], kind, target.replace(os.sep, "/"))
        self.errors.append("%s: link to a missing file: %s" % (src, href))
        return href

    def rewrite_img(self, tag, src):
        m = re.search(r'src="([^"]+)"', tag)
        if not m or re.match(r"^https?:", m.group(1)):
            return tag
        target = os.path.normpath(os.path.join(os.path.dirname(src), m.group(1)))
        if target.endswith(".svg") and os.path.exists(os.path.join(ROOT, target[:-4] + ".png")):
            target = target[:-4] + ".png"            # the text browser of the GUI does not always have an SVG decoder
        if not os.path.exists(os.path.join(ROOT, target)):
            target = os.path.join("docs", "img", os.path.basename(m.group(1)))      # the pages of docs/manual find the images in docs/img
        if not os.path.exists(os.path.join(ROOT, target)):
            self.errors.append("%s: image not found: %s" % (src, m.group(1)))
            return tag
        name = "img/" + os.path.basename(target)
        if self.images.get(name, target) != target:
            self.errors.append("two different images are called %s" % name)
        self.images[name] = target
        tag = tag.replace(m.group(0), 'src="%s"' % name)
        size = png_size(os.path.join(ROOT, target))
        if size and "width=" not in tag:
            tag = tag.replace("<img ", '<img width="%d" ' % min(size[0], 760), 1)
        return tag

    # ---- one page --------------------------------------------------------------------------------------------------
    def render(self, index):
        page = self.pages[index]
        src = page["src"]
        text = expand_includes(read(src), self.errors, src)
        md = markdown.Markdown(extensions=["tables", "fenced_code", "toc", "sane_lists"], extension_configs={"toc": {"toc_depth": "2-3"}})
        body = md.convert(text)
        body = re.sub(r'href="([^"]*)"', lambda m: 'href="%s"' % self.rewrite_href(html.unescape(m.group(1)), src), body)
        body = re.sub(r"<img [^>]*>", lambda m: self.rewrite_img(m.group(0), src), body)
        body = body.replace("<table>", '<table border="1" cellspacing="0" cellpadding="5" width="100%">')
        # anchors that the Qt text browser also understands
        body = re.sub(r'<(h[1-4]) id="([^"]+)">', r'<\1 id="\2"><a name="\2"></a>', body)
        heads = []
        def walk(tokens):
            for t in tokens:
                heads.append({"level": t["level"], "text": t["name"], "id": t["id"]})
                walk(t["children"])
        walk(md.toc_tokens)
        title = page["title"]
        prev_p = self.pages[index - 1] if index > 0 else None
        next_p = self.pages[index + 1] if index + 1 < len(self.pages) else None
        nav = ['<a href="index.html">Home</a>', '<a href="index.html#contents">Contents</a>']
        if prev_p: nav.append('&laquo; <a href="%s.html">%s</a>' % (prev_p["id"], html.escape(prev_p["title"])))
        if next_p: nav.append('<a href="%s.html">%s</a> &raquo;' % (next_p["id"], html.escape(next_p["title"])))
        nav_html = '<p class="nav">' + " &middot; ".join(nav) + "</p>"
        doc = ('<!DOCTYPE html>\n<html lang="en">\n<head>\n<meta charset="utf-8">\n<title>%s - ligint manual</title>\n<style>%s</style>\n</head>\n<body>\n%s\n%s\n%s\n'
               '<p class="foot">ligint %s manual</p>\n</body>\n</html>\n') % (html.escape(title), CSS, nav_html, body, nav_html, self.version)
        return doc, {"id": page["id"], "title": title, "file": page["id"] + ".html", "headings": [h for h in heads if h["level"] >= 2]}

    def build(self, out):
        os.makedirs(os.path.join(out, "img"), exist_ok=True)
        toc = []
        for i, page in enumerate(self.pages):
            doc, entry = self.render(i)
            toc.append(entry)
            with open(os.path.join(out, page["id"] + ".html"), "w", encoding="utf-8", newline="\n") as f:
                f.write(doc)
        for name, target in sorted(self.images.items()):
            shutil.copyfile(os.path.join(ROOT, target), os.path.join(out, name))
        with open(os.path.join(out, "toc.json"), "w", encoding="utf-8", newline="\n") as f:
            json.dump({"title": self.m["title"], "version": self.version, "pages": toc}, f, indent=1)
            f.write("\n")
        files = sorted(["toc.json"] + [p["id"] + ".html" for p in self.pages] + list(self.images))
        with open(os.path.join(out, "manual.qrc"), "w", encoding="utf-8", newline="\n") as f:
            f.write('<RCC>\n    <qresource prefix="/manual">\n')
            for name in files:
                f.write('        <file alias="%s">%s</file>\n' % (name, name))
            f.write("    </qresource>\n</RCC>\n")

def tree(path):
    out = {}
    for d, _, files in os.walk(path):
        for n in files:
            full = os.path.join(d, n)
            with open(full, "rb") as f:
                out[os.path.relpath(full, path).replace(os.sep, "/")] = f.read()
    return out

def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--check", action="store_true", help="do not write: fail if docs/manual/html is out of date")
    args = ap.parse_args()
    with open(os.path.join(MANUAL, "manifest.json"), encoding="utf-8") as f:
        manifest = json.load(f)
    b = Builder(manifest)
    tmp = tempfile.mkdtemp(prefix="ligint_manual_")
    try:
        b.build(tmp)
        if b.errors:
            print("\n".join(sorted(set(b.errors))), file=sys.stderr)
            return 1
        if args.check:
            if markdown.__version__ != MARKDOWN_VERSION:
                print("cannot check: the HTML was generated with markdown %s and this is %s (pip install -r tools/requirements-manual.txt)" %
                      (MARKDOWN_VERSION, markdown.__version__), file=sys.stderr)
                return 2
            new, old = tree(tmp), tree(OUT) if os.path.isdir(OUT) else {}
            stale = sorted(n for n in set(new) | set(old) if new.get(n) != old.get(n))
            if stale:
                print("docs/manual/html is out of date (run python3 tools/build_manual.py): " + ", ".join(stale[:8]), file=sys.stderr)
                return 1
            print("manual up to date: %d files" % len(new))
            return 0
        if os.path.isdir(OUT):
            shutil.rmtree(OUT)
        shutil.copytree(tmp, OUT)
        print("manual written: %d pages, %d images -> %s" % (len(manifest["pages"]), len(b.images), os.path.relpath(OUT, ROOT)))
        return 0
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

if __name__ == "__main__":
    sys.exit(main())
