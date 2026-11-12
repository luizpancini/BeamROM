import argparse

import numpy as np
from pyNastran.bdf.bdf import BDF
from pyNastran.bdf.case_control_deck import CaseControlDeck

from common import case_dir, write_bdf_if_changed, write_meta


def _node_id(i, j, k, nL, nw):
    return 1 + i + j * (nL + 1) + k * (nL + 1) * (nw + 1)


def _build_solid_beam(L, w, t, nL, nw, nt, E, nu, rho):
    """Mesh, material and property shared by all analyses."""
    model = BDF()

    # Material and property
    model.add_mat1(1, E=E, G=None, nu=nu, rho=rho)
    model.add_psolid(1, 1)

    xs = np.linspace(0.0, L, nL + 1)
    ys = np.linspace(0.0, w, nw + 1)
    zs = np.linspace(0.0, t, nt + 1)

    for k, z in enumerate(zs):
        for j, y in enumerate(ys):
            for i, x in enumerate(xs):
                model.add_grid(_node_id(i, j, k, nL, nw), [x, y, z])

    # CHEXA8 elements (bottom face counter-clockwise, then top face)
    def nid(i, j, k):
        return _node_id(i, j, k, nL, nw)

    eid = 1
    for k in range(nt):
        for j in range(nw):
            for i in range(nL):
                nodes = [
                    nid(i,     j,     k),     nid(i + 1, j,     k),
                    nid(i + 1, j + 1, k),     nid(i,     j + 1, k),
                    nid(i,     j,     k + 1), nid(i + 1, j,     k + 1),
                    nid(i + 1, j + 1, k + 1), nid(i,     j + 1, k + 1),
                ]
                model.add_chexa(eid, 1, nodes)
                eid += 1

    return model


def solid_beam_static(L, w, t, nL, nw, nt, E, nu, rho, filename="solid_beam_static.bdf"):
    """Solid beam with six static load cases (unit forces and moments at the tip centre).
    L: length (x), w: width (y), t: thickness (z), in meters.
    nL, nw, nt: number of elements in x, y, z."""

    model = _build_solid_beam(L, w, t, nL, nw, nt, E, nu, rho)

    # Load point: extra node at the centre of the tip section (x = L), tied to the
    # tip face with an RBE3 (interpolation only, adds no stiffness)
    n_el = nL * nw * nt
    load_node = _node_id(nL, nw, nt, nL, nw) + 1
    model.add_grid(load_node, [L, w / 2, t / 2])
    tip_face = [_node_id(nL, j, k, nL, nw) for k in range(nt + 1) for j in range(nw + 1)]
    model.add_rbe3(n_el + 1, load_node, "123456", weights=[1.0], comps=["123"], Gijs=[tip_face])

    # Root clamp: node at the centre of the root section (x = 0), rigidly tied to the root face (RBE2) and fixed in all DOF, so its SPC force is the total root reaction (F and M)
    root_node = load_node + 1
    model.add_grid(root_node, [0., w / 2, t / 2])
    root_face = [_node_id(0, j, k, nL, nw) for k in range(nt + 1) for j in range(nw + 1)]
    model.add_rbe2(n_el + 2, root_node, "123", root_face)
    model.add_spc1(1, "123456", [root_node])

    _add_static_loads(model, root_node, load_node)

    write_bdf_if_changed(model, filename)

    return model, root_node, load_node


def _add_static_loads(model, root_node, load_node):
    """Six static subcases (SOL 101), each with one unit force or moment (N and N*m) at the load node, with the root node fixed (SPC 1)."""
    loads = [
        ("FORCE", [1., 0., 0.]), ("FORCE", [0., 1., 0.]), ("FORCE", [0., 0., 1.]),
        ("MOMENT", [1., 0., 0.]), ("MOMENT", [0., 1., 0.]), ("MOMENT", [0., 0., 1.]),
    ]
    case_lines = [
        "ECHO = NONE",
        f"SET 2 = {root_node}",
        "SPCFORCES = 2",  # root reaction (F and M) at the root centre node
    ]
    for sid, (kind, direction) in enumerate(loads, start=1):
        if kind == "FORCE":
            model.add_force(sid, load_node, 1.0, direction)
        else:
            model.add_moment(sid, load_node, 1.0, direction)
        case_lines += [
            f"SUBCASE {sid}",
            f"  LOAD = {sid}",
            "  SPC = 1",
            "  DISPLACEMENT(PLOT) = ALL",
        ]

    # Static analysis
    model.sol = 101
    model.case_control_deck = CaseControlDeck(case_lines)
    model.add_param("POST", [-1])  # generates .op2


def shell_beam_static(L, w, t, nL, nw, E, nu, rho, filename="shell_beam_static.bdf"):
    """Shell beam (CQUAD4 mid-surface in the x-y plane, thickness t along z) with the same six static load cases as the solid beam. 
    L: length (x), w: width (y), t: thickness (z), in meters.
    nL, nw: number of elements in x, y."""

    # Create a new BDF model for the shell beam
    model = BDF()

    # Define material with given properties (E, nu, rho)
    model.add_mat1(1, E=E, G=None, nu=nu, rho=rho)

    # Membrane, bending and transverse shear (shear correction factor 5/6)
    model.add_pshell(1, mid1=1, t=t, mid2=1, mid3=1, tst=5 / 6)

    # Function to compute node ID based on i, j indices in the mesh
    def nid(i, j):
        return 1 + i + j * (nL + 1)

    # Add grid points for the shell mesh
    for j, y in enumerate(np.linspace(0.0, w, nw + 1)):
        for i, x in enumerate(np.linspace(0.0, L, nL + 1)):
            model.add_grid(nid(i, j), [x, y, 0.0])

    # Add shell elements (CQUAD4) for the mesh
    eid = 1
    for j in range(nw):
        for i in range(nL):
            model.add_cquad4(eid, 1, [nid(i, j), nid(i + 1, j), nid(i + 1, j + 1), nid(i, j + 1)])
            eid += 1
    n_el = eid - 1

    # Load point at the centre of the tip edge, tied to it with an RBE3. Moments about x and z are carried by the translations (force couples) and moments about y by the edge rotations about y (the nodes are collinear in the plane, so translations cannot carry them). The drilling rotation (about z) must not be used: it has almost no stiffness and would absorb Mz
    load_node = nid(nL, nw) + 1
    model.add_grid(load_node, [L, w / 2, 0.0])
    tip_edge = [nid(nL, j) for j in range(nw + 1)]
    model.add_rbe3(n_el + 1, load_node, "123456", weights=[1.0, 1.0], comps=["123", "5"], Gijs=[tip_edge, tip_edge])

    # Root clamp: node at the centre of the root edge, rigidly tied to it and fixed
    root_node = load_node + 1
    model.add_grid(root_node, [0., w / 2, 0.0])
    root_edge = [nid(0, j) for j in range(nw + 1)]
    model.add_rbe2(n_el + 2, root_node, "123456", root_edge)
    model.add_spc1(1, "123456", [root_node])

    _add_static_loads(model, root_node, load_node)

    write_bdf_if_changed(model, filename)

    return model, root_node, load_node


def line_beam_static(L, w, t, nL, E, nu, rho, filename="line_beam_static.bdf"):
    """Beam elements (CBAR) along x with a w x t rectangular section (w along y, t along z) and the same six static load cases as the solid beam. nL: number of elements."""

    # Create a new BDF model for the line beam
    model = BDF()

    # Define material with given properties (E, nu, rho)
    model.add_mat1(1, E=E, G=None, nu=nu, rho=rho)

    # Section properties. Plane 1 is x-y (bending about z), plane 2 is x-z (bending about y). Torsion constant of a rectangle (w >= t), shear correction factor 5/6 in both planes
    A = w * t
    I1 = t * w**3 / 12
    I2 = w * t**3 / 12
    J = w * t**3 * (1 / 3 - 0.21 * (t / w) * (1 - t**4 / (12 * w**4)))
    model.add_pbar(1, 1, A=A, i1=I1, i2=I2, j=J, k1=5 / 6, k2=5 / 6)

    # Add grid points along the beam axis (x-direction)
    for i, x in enumerate(np.linspace(0.0, L, nL + 1)):
        model.add_grid(i + 1, [x, 0.0, 0.0])

    # Orientation vector along y: element axes coincide with the global axes
    for i in range(nL):
        model.add_cbar(i + 1, 1, [i + 1, i + 2], x=[0., 1., 0.], g0=None)

    root_node, load_node = 1, nL + 1
    model.add_spc1(1, "123456", [root_node])

    _add_static_loads(model, root_node, load_node)

    write_bdf_if_changed(model, filename)

    return model, root_node, load_node


def solid_beam_modal(L, w, t, nL, nw, nt, E, nu, rho, nd=10, filename="solid_beam_modal.bdf"):
    """Same solid beam as solid_beam_static, with a modal analysis (SOL 103).
    nd: number of modes requested."""

    model = _build_solid_beam(L, w, t, nL, nw, nt, E, nu, rho)

    # Clamped at the root (x = 0): solids only have 3 translational DOF
    root = [_node_id(0, j, k, nL, nw) for k in range(nt + 1) for j in range(nw + 1)]
    model.add_spc1(1, "123", root)

    # Modal analysis
    model.add_eigrl(10, nd=nd, norm="MASS")
    model.sol = 103
    model.case_control_deck = CaseControlDeck([
        "ECHO = NONE",
        "SUBCASE 1",
        "  SPC = 1",
        "  METHOD = 10",
        "  DISPLACEMENT(PLOT) = ALL",
    ])
    model.add_param("POST", [-1])  # generates .op2

    write_bdf_if_changed(model, filename)

    return model


# Element types: case directory name (nastran/simple_beams/<name>/mesh<N>), static model builder and the number of elements along the length (nL), width (nw) and thickness (nt) of each mesh
ELEMENTS = {
    "solid": dict(name="solid_beam", static=solid_beam_static, meshes={
        1: dict(nL=100, nw=20, nt=5),
        2: dict(nL=100, nw=25, nt=10),
        3: dict(nL=100, nw=20, nt=10),
        4: dict(nL=200, nw=20, nt=5),
    }),
    "shell": dict(name="shell_beam", static=shell_beam_static, meshes={
        1: dict(nL=50, nw=10),
        2: dict(nL=100, nw=20),
        3: dict(nL=200, nw=40),
    }),
    "line": dict(name="line_beam", static=line_beam_static, meshes={
        1: dict(nL=10),
        2: dict(nL=20),
        3: dict(nL=50),
        4: dict(nL=100),
    }),
}


def main():
    parser = argparse.ArgumentParser(description="Generate the beam BDF files into nastran/simple_beams/<element>_beam/mesh<N>/")
    parser.add_argument("--element", choices=sorted(ELEMENTS), default="solid", help="solid: CHEXA, shell: CQUAD4, line: CBAR")
    parser.add_argument("--mesh", type=int, required=True)
    parser.add_argument("--analyses", nargs="+", choices=["static", "modal"], help="default: static and modal for solid, static for shell and line")
    parser.add_argument("--L", type=float, default=0.55, help="length [m]")
    parser.add_argument("--w", type=float, default=60e-3, help="width [m]")
    parser.add_argument("--t", type=float, default=2.25e-3, help="thickness [m]")
    parser.add_argument("--E", type=float, default=71e9, help="Young's modulus [Pa]")
    parser.add_argument("--nu", type=float, default=0.33, help="Poisson's ratio")
    parser.add_argument("--rho", type=float, default=2.795e3, help="density [kg/m^3]")
    args = parser.parse_args()

    element = ELEMENTS[args.element]
    if args.mesh not in element["meshes"]:
        parser.error(f"--mesh must be one of {sorted(element['meshes'])} for {args.element} elements")
    analyses = args.analyses or (["static", "modal"] if args.element == "solid" else ["static"])
    if "modal" in analyses and args.element != "solid":
        parser.error("modal analysis is only available for solid elements")

    props = dict(L=args.L, w=args.w, t=args.t, **element["meshes"][args.mesh], E=args.E, nu=args.nu, rho=args.rho)
    name = element["name"]
    directory = case_dir(name, args.mesh)
    directory.mkdir(parents=True, exist_ok=True)

    if "modal" in analyses:
        solid_beam_modal(**props, filename=directory / f"{name}_modal.bdf")
    if "static" in analyses:
        _, root_node, load_node = element["static"](**props, filename=directory / f"{name}_static.bdf")
        write_meta(directory, mesh=args.mesh, root_node=root_node, load_node=load_node, **props)
        print("Root node:", root_node)
        print("Load node:", load_node)
    print("BDF files in", directory)


if __name__ == "__main__":
    main()