"""
Draw N-glycan cartoons in the Symbol Nomenclature for Glycans (SNFG) with matplotlib.

Topologies follow the common glycan names used in the paper (G0F, A2G2S1, A3G3S2, Man5, ...).
Residues whose attachment site was not determined are drawn outside a bracket, and targets named
only by composition (for example N4H6F1G1) are drawn as compositions with generic symbols.
Linkages are not shown. The reducing end is on the right.
"""

import re

from matplotlib.patches import Circle, Polygon, Rectangle

# SNFG colors (Varki et al., Glycobiology 2015; Neelamegham et al., Glycobiology 2019)
BLUE, GREEN, YELLOW, RED = "#0072BC", "#00A651", "#FFD400", "#ED1C24"
PURPLE, LIGHTBLUE, WHITE = "#A54399", "#8FCCE9", "#FFFFFF"
STYLE = {"GlcNAc": ("square", BLUE), "Man": ("circle", GREEN), "Gal": ("circle", YELLOW),
         "Fuc": ("triangle", RED), "NeuAc": ("diamond", PURPLE), "NeuGc": ("diamond", LIGHTBLUE),
         "HexNAc": ("square", WHITE), "Hex": ("circle", WHITE), "dHex": ("triangle", WHITE)}
R = 0.29  # half-size of a symbol, in residue spacings
LW = 0.5


def symbol(ax, kind, x, y, r=R):
    shape, color = STYLE[kind]
    kw = dict(facecolor=color, edgecolor="black", lw=LW, zorder=3)
    if shape == "square":
        p = Rectangle((x - r, y - r), 2 * r, 2 * r, **kw)
    elif shape == "circle":
        p = Circle((x, y), r * 1.08, **kw)
    elif shape == "triangle":
        p = Polygon([(x - r * 1.1, y - r), (x + r * 1.1, y - r), (x, y + r)], closed=True, **kw)
    else:
        p = Polygon([(x - r * 1.2, y), (x, y + r * 1.2), (x + r * 1.2, y), (x, y - r * 1.2)], closed=True, **kw)
    ax.add_patch(p)


def bond(ax, p, q):
    ax.plot([p[0], q[0]], [p[1], q[1]], color="black", lw=LW, zorder=1, solid_capstyle="butt")


def parse(name, composition=None):
    """Return a drawing spec for a target name; fall back to the composition."""
    base = re.sub(r"_(a|b|early|late)$", "", name)
    if base == "Man5":
        return {"type": "man5"}
    m = re.fullmatch(r"G([012])(F?)", base)
    if m:
        return {"type": "complex", "ant": 2, "gal": int(m.group(1)), "sia": 0, "fuc": bool(m.group(2))}
    m = re.fullmatch(r"(F?)A([23])G([0-3])(F?)(?:S([0-4]))?", base)
    if m:
        return {"type": "complex", "ant": int(m.group(2)), "gal": int(m.group(3)), "sia": int(m.group(5) or 0),
                "fuc": bool(m.group(1) or m.group(4))}
    if composition is None:
        raise ValueError(f"no structure or composition for {name}")
    return {"type": "composition", "counts": [(k, int(v)) for k, v in re.findall(r"([NHFSG])(\d+)", composition)]}


def _core(ax, fuc):
    """Chitobiose core and the branching mannose; returns the position of the branching mannose."""
    pts = [(0, 0), (-1, 0), (-2, 0)]
    bond(ax, pts[0], pts[1])
    bond(ax, pts[1], pts[2])
    symbol(ax, "GlcNAc", *pts[0])
    symbol(ax, "GlcNAc", *pts[1])
    symbol(ax, "Man", *pts[2])
    if fuc:
        bond(ax, (0, 0), (0, -0.95))
        symbol(ax, "Fuc", 0, -0.95)
    return pts[2]


def _bracket(ax, x, y0, y1, items):
    """Draw a bracket at x spanning y0..y1 and the chains in items to its left."""
    pad = 0.45
    ax.plot([x + 0.18, x, x, x + 0.18], [y1 + pad, y1 + pad, y0 - pad, y0 - pad], color="black", lw=LW, zorder=1)
    n = len(items)
    ys = [((n - 1) / 2 - i) * 0.85 + (y0 + y1) / 2 for i in range(n)]
    left = x
    for chain, y in zip(items, ys):
        xs = [x - 0.6 - k for k in range(len(chain))]
        for k in range(len(chain) - 1):
            bond(ax, (xs[k], y), (xs[k + 1], y))
        for kind, xx in zip(chain, xs):
            symbol(ax, kind, xx, y)
        left = min(left, xs[-1])
    return left


def _complex(ax, spec):
    man = _core(ax, spec["fuc"])
    if spec["ant"] == 2:
        arms = [((-3, 0.8), [0.8]), ((-3, -0.8), [-0.8])]
    else:
        arms = [((-3, 0.9), [0.9]), ((-3, -0.75), [-0.3, -1.2])]
    ant_y = []
    for pos, ys in arms:
        bond(ax, man, pos)
        symbol(ax, "Man", *pos)
        for y in ys:
            bond(ax, pos, (-4, y))
            ant_y.append(y)
    n, gal, sia = spec["ant"], spec["gal"], spec["sia"]
    full_gal = gal == n
    full_sia = full_gal and sia >= n
    left = -4
    for y in ant_y:
        chain = ["GlcNAc"] + (["Gal"] if full_gal else []) + (["NeuAc"] if full_sia else [])
        for k in range(1, len(chain)):
            bond(ax, (-4 - k + 1, y), (-4 - k, y))
        for k, kind in enumerate(chain):
            symbol(ax, kind, -4 - k, y)
        left = min(left, -4 - len(chain) + 1)
    items = []
    if not full_gal:
        items = [["Gal", "NeuAc"]] * min(sia, gal) + [["Gal"]] * (gal - min(sia, gal))
    elif not full_sia:
        items = [["NeuAc"]] * sia
    else:
        items = [["NeuAc"]] * (sia - n)
    if items:
        left = _bracket(ax, left - 0.6, min(ant_y), max(ant_y), items)
    return left


def _man5(ax):
    man = _core(ax, False)
    up, down = (-3, 0.8), (-3, -0.8)
    for p in (up, down):
        bond(ax, man, p)
        symbol(ax, "Man", *p)
    for p in ((-4, 1.25), (-4, 0.35)):
        bond(ax, up, p)
        symbol(ax, "Man", *p)
    return -4


def _composition(ax, counts):
    kinds = {"N": "HexNAc", "H": "Hex", "F": "dHex", "S": "NeuAc", "G": "NeuGc"}
    x = 0
    for k, v in reversed(counts):
        ax.text(x + 0.05, -0.32, str(v), fontsize=5.2, ha="left", va="center", zorder=4)
        symbol(ax, kinds[k], x - 0.35, 0)
        x -= 1.35
    return x + 0.6


def draw(ax, name, composition=None):
    """Draw the cartoon for one target on ax (data units: one residue spacing). Returns the left edge."""
    spec = parse(name, composition)
    if spec["type"] == "man5":
        return _man5(ax)
    if spec["type"] == "complex":
        return _complex(ax, spec)
    return _composition(ax, spec["counts"])


def save_png(path, name, composition=None, unit_in=0.11, ylim=(-1.75, 1.6), dpi=600):
    """Write a single cartoon as a tight PNG (used for the Supporting Information tables)."""
    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax = fig.add_axes([0, 0, 1, 1])
    left = draw(ax, name, composition)
    ax.set_xlim(left - 0.5, 0.75)
    ax.set_ylim(*ylim)
    ax.set_aspect("equal")
    ax.axis("off")
    fig.set_size_inches((0.75 - left + 0.5) * unit_in, (ylim[1] - ylim[0]) * unit_in)
    fig.savefig(path, dpi=dpi, transparent=False, facecolor="white")
    plt.close(fig)
