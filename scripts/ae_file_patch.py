"""Changes one field of an all-electron data file (ae_file_m.f90), in place,
for ae_file_error_test.sh: FIELD version reverses the bytes of the version
(another byte order); zz, fa (the last occupancy), rr (the first mesh
point after the first) and c_light are multiplied by FACTOR.

usage: ae_file_patch.py FILE FIELD [FACTOR]
"""

import struct
import sys

path, field = sys.argv[1:3]
with open(path, "rb") as f:
    buf = bytearray(f.read())
zz = 8 + 4 + 4 + 2
nc, nv = struct.unpack_from("=ii", buf, zz + 8)
nt = nc + nv
fa = zz + 8 + 12 + 8 * nt
rr = fa + 8 * nt + 4
(mmax,) = struct.unpack_from("=i", buf, rr - 4)
at = {"zz": zz, "fa": fa + 8 * (nt - 1), "rr": rr + 8, "c_light": rr + 8 * mmax}
if field == "version":
    buf[8:12] = buf[8:12][::-1]
else:
    (x,) = struct.unpack_from("=d", buf, at[field])
    struct.pack_into("=d", buf, at[field], x * float(sys.argv[3]))
with open(path, "wb") as f:
    f.write(buf)
