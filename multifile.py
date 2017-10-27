import re

from cStringIO import StringIO


_FASTA_FIRST_LINE_RE = re.compile(r'>gi\|\d+\|ref\|')
def read(fobj):
    s = StringIO(fobj.read())
    first = s.readline()
    if _FASTA_FIRST_LINE_RE.match(first):
        so = StringIO()
        for line in s:
            so.write(line.strip())
        so.reset()
        return so
    s.reset()
    return s
