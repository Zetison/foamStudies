import click
from collections import OrderedDict
from math import ceil, log, sqrt, pi
import numpy as np
from lxml import etree as xml
import sys

from splipy import curve_factory as cf, surface_factory as sf
from splipy.io import G2
from splipy.utils.refinement import geometric_refine
from splipy.volume_factory import extrude


def graded_space(start, step, factor, N):
    def gen_graded_space(start, step, factor, N):
        for _ in range(N):
            yield start
            start += step
            step *= factor
    return list(gen_graded_space(float(start), float(step), float(factor), N))


def find_factor(initial, total, N, tol=1e-7):
    # Solve initial * (1 - a^N) / (1 - a) = total
    overshoot = lambda a: initial * (1 - (1 + a)**N) / (1 - (1 + a)) - total
    l, u = 0.5, 0.5
    while overshoot(l) > 0:
        l /= 2
    while overshoot(u) < 0:
        u *= 2
    while True:
        m = (u + l) / 2
        s = overshoot(m)
        if abs(s) < tol:
            break
        if s > 0:
            u = m
        else:
            l = m
    return 1 + m


class PatchDict(OrderedDict):

    def __init__(self, dim, *args, **kwargs):
        super(PatchDict, self).__init__(*args, **kwargs)
        self.dim = dim
        self.masters = {}
        self.boundaries = {}
        self.periodics = []

    def add(self, *args):
        names, patches = args[:-1], args[-1]
        for name, patch in zip(names, patches):
            self[name] = patch

    def connect(self, *args):
        for master, medge, slave, sedge, *rest in args:
            rev = 'rev' in rest
            per = 'per' in rest

            if master == slave and per:
                self.periodics.append((master, max(medge, sedge) // 2))
            elif master != slave:
                self.masters[(master, medge)] = (slave, sedge, rev, per)
            else:
                raise Exception('What?')

    def boundary(self, name, patch, number, dim=-1, add=0):
        kind = {
            2: ['vertex', 'edge'],
            3: ['vertex', 'edge', 'face'],
        }[self.dim][dim]
        number += add
        self.boundaries.setdefault(name, {}).setdefault(kind, []).append((patch, number))

    def write(self, fn, order=4):
        pids, patches = {}, []
        for i, (name, patch) in enumerate(self.items()):
            pids[name] = i
            diff = [o - order for o in patch.order()]
            patches.append(patch.lower_order(*diff))

        with G2(fn + '.g2') as f:
            f.write(patches)

        dom = xml.Element('geometry')
        dom.attrib['dim'] = str(patches[0].pardim)
        xml.SubElement(dom, 'patchfile').text = fn + '.g2'

        topology = xml.SubElement(dom, 'topology')
        for (master, medge), (slave, sedge, rev, periodic) in self.masters.items():
            mid = pids[master] + 1
            sid = pids[slave] + 1
            if mid > sid:
                mid, sid = sid, mid
                medge, sedge = sedge, medge
            xml.SubElement(topology, 'connection').attrib.update({
                'master': str(mid),
                'midx': str(medge),
                'slave': str(sid),
                'sidx': str(sedge),
                'reverse': 'true' if rev else 'false',
                'periodic': 'true' if periodic else 'false',
            })
        for patch, direction in self.periodics:
            pid = pids[patch] + 1
            xml.SubElement(topology, 'periodic').attrib.update({
                'patch': str(pid),
                'dir': str(direction),
            })

        topsets = xml.SubElement(dom, 'topologysets')
        for name, kinds in self.boundaries.items():
            for kind, items in kinds.items():
                topset = xml.SubElement(topsets, 'set')
                topset.attrib.update({'name': name, 'type': kind})
                for patch, number in items:
                    item = xml.SubElement(topset, 'item')
                    item.attrib['patch'] = str(pids[patch] + 1)
                    item.text = str(number)

        xml.ElementTree(dom).write(fn + '.xinp', encoding='utf-8', xml_declaration=True, pretty_print=True, standalone=False)


@click.command()
@click.option('--diam', default=1.0)
@click.option('--width', default=20.0)
@click.option('--front', default=20.0)
@click.option('--back', default=40.0)
@click.option('--side', default=20.0)
@click.option('--height', default=0.0)
@click.option('--Re', default=100.0)
@click.option('--grad', type=float, required=False)
@click.option('--nel-bndl', default=10)
@click.option('--inner-elsize', type=float, required=False)
@click.option('--nel-side', type=int, required=False)
@click.option('--nel-circ', default=40)
@click.option('--nel-height', default=10)
@click.option('--order', default=4)
@click.option('--outer-graded/--no-outer-graded', default=True)
@click.option('--thickness', default=4.0)
@click.option('--out', default='out')
def cylinder(diam, width, front, back, side, height, re, grad, inner_elsize,
             nel_side, nel_bndl, nel_circ, nel_height, order, out, outer_graded, thickness):
    assert all(f >= width for f in [front, back, side])

    rad_cyl = diam / 2
    width *= rad_cyl
    back = back * rad_cyl - width
    front = front * rad_cyl - width
    side = side * rad_cyl - width
    elsize = width * 2 / nel_circ

    dim = 2 if height == 0.0 else 3
    vx_add = 8 if dim == 3 else 0
    patches = PatchDict(dim)

    if inner_elsize and nel_side:
        # Calculate grading factor based on first element size,
        # total length and number of elements
        dr = rad_cyl * inner_elsize
        grad = find_factor(dr, width - rad_cyl, nel_side)

    elif nel_bndl and grad:
        # Calculate first element size based on total length
        # and number of elements

        # We want nel_bndl elements inside the boundary layer
        # Calculate how small the inner element must be
        size_bndl = 1 / sqrt(re) * diam
        dr = (1 - grad) / (1 - grad ** nel_bndl) * size_bndl

        # Potentially reduce element size so we get a whole number of elements
        # on either side of the cylinder
        #nel_side = int(ceil(log(1 - 1/dr * (1 - grad) * (width - rad_cyl)) / log(grad)))
        nel_side = int(round(log(1 - 1/dr * (1 - grad) * (width - rad_cyl)) / log(grad)))
        dr = (1 - grad) / (1 - grad ** nel_side) * (width - rad_cyl)
    else:
        print('Specify (inner-elsize and nel-side) or (nel-bndl and grad)', file=sys.stderr)
        sys.exit(1)

    # Graded radial space from cylinder to edge of domain
    radial_kts = graded_space(rad_cyl, dr, grad, nel_side) + [width]

    # Create a radial and divide it
    radial = cf.cubic_curve(np.matrix(radial_kts).T, boundary=cf.Boundary.NATURAL)
    radial.set_dimension(3)
    radial_kts = radial.knots('u')
    # middle = radial_kts[len(radial_kts) // 2]
    middle = next(k for k in radial_kts if k>=thickness*diam)
    radial_inner, radial_outer = radial.split(middle, 'u')
    radial_inner.rotate(pi/4)
    dl = np.linalg.norm(radial_outer(radial_outer.knots('u')[-2]) - radial_outer.section(u=-1)) * grad

    # Revolve the inner radial and divide it
    radials = [radial_inner.clone().rotate(v) for v in np.linspace(0, 2*pi, 4*nel_circ + 1)]
    inner = sf.loft(radials)
    ikts = inner.knots('v')
    inner.insert_knot((ikts[0] + ikts[1]) / 2, 'v')
    inner.insert_knot((ikts[-1] + ikts[-2]) / 2, 'v')

    ikts = inner.knots('v')
    patches.add('iu', 'il', 'id', 'ir', inner.split([ikts[k*nel_circ] for k in range(5)][1:-1], 'v'))
    patches.connect(
        ('ir', 4, 'iu', 3),
        ('iu', 4, 'il', 3),
        ('il', 4, 'id', 3),
        ('id', 4, 'ir', 3),
    )
    patches.boundary('cylinder', 'ir', 1)
    patches.boundary('cylinder', 'iu', 1)
    patches.boundary('cylinder', 'il', 1)
    patches.boundary('cylinder', 'id', 1)

    # Create an outer section
    rc = radial_outer.section(u=0)
    alpha = (sqrt(2) * width - rc[0]) / (width - rc[0])
    right = ((radial_outer - rc) * alpha + rc).rotate(pi/4)
    left = right.clone().rotate(pi/2).reverse()
    outer = cf.line((-width, width, 0), (width, width, 0))
    outer.set_order(4).refine(nel_circ - 1)
    inner = patches['iu'].section(u=-1)
    outer = sf.edge_curves(right, outer, left, inner).reverse('v')

    patches.add('ou', 'ol', 'od', 'or', [outer.clone().rotate(v) for v in [0, pi/2, pi, 3*pi/2]])
    patches.connect(
        ('or', 4, 'ou', 3),
        ('ou', 4, 'ol', 3),
        ('ol', 4, 'od', 3),
        ('od', 4, 'or', 3),
        ('or', 1, 'ir', 2),
        ('ou', 1, 'iu', 2),
        ('ol', 1, 'il', 2),
        ('od', 1, 'id', 2),
    )

    if front > 0:
        la = patches['ol'].section(u=-1)
        lb = la.clone() - (front, 0, 0)
        front_srf = sf.edge_curves(lb, la).set_order(4,4).swap().reverse('v')
        nel = int(ceil(log(1 - 1/dl * (1 - grad) * front) / log(grad)))
        geometric_refine(front_srf, grad, nel - 1, reverse=True)
        patches['fr'] = front_srf
        patches.connect(('fr', 2, 'ol', 2, 'rev'))
        patches.boundary('inflow', 'fr', 1)
    else:
        patches.boundary('inflow', 'ol', 2)
        patches.boundary('inflow', 'ou', 4, dim=-2, add=vx_add)
        patches.boundary('inflow', 'od', 2, dim=-2, add=vx_add)

    if back > 0:
        la = patches['or'].section(u=-1)
        lb = la.clone() + (back, 0, 0)
        back_srf = sf.edge_curves(la, lb).set_order(4,4).swap()
        if outer_graded:
            nel = int(ceil(log(1 - 1/dl * (1 - grad) * back) / log(grad)))
            geometric_refine(back_srf, grad, nel - 1)
        else:
            #nel = int(ceil(back / dl))
            nel = int(round(back / (2*width) * nel_circ))
            back_srf.refine(nel - 1, direction='u')
        patches['ba'] = back_srf
        patches.connect(('ba', 1, 'or', 2))
        patches.boundary('outflow', 'ba', 2)
    else:
        patches.boundary('outflow', 'or', 2)

    if side > 0:
        la = patches['ou'].section(u=-1).reverse()
        lb = la + (0, side, 0)
        patches['up'] = sf.edge_curves(la, lb).set_order(4,4)
        patches.connect(('up', 3, 'ou', 2, 'rev'), ('dn', 4, 'od', 2))
        patches.boundary('top', 'up', 4)
        patches.boundary('bottom', 'dn', 3)

        if 'fr' in patches:
            btm = front_srf.section(v=-1)
            right = patches['up'].section(u=0)
            top = (btm + (0, side, 0)).reverse()
            left = (right - (front, 0, 0)).reverse()
            patches['upfr'] = sf.edge_curves(btm, right, top, left)
            patches.connect(
                ('upfr', 3, 'fr', 4), ('upfr', 2, 'up', 1),
                ('dnfr', 4, 'fr', 3), ('dnfr', 2, 'dn', 1),
            )
            patches.boundary('wall', 'upfr', 4)
            patches.boundary('inflow', 'upfr', 1)
            patches.boundary('wall', 'dnfr', 3)
            patches.boundary('inflow', 'dnfr', 1)
        else:
            patches.boundary('inflow', 'up', 1)
            patches.boundary('inflow', 'dn', 1)

        if 'ba' in patches:
            btm = back_srf.section(v=-1)
            left = patches['up'].section(u=-1).reverse()
            top = (btm + (0, side, 0)).reverse()
            right = (left + (back, 0, 0)).reverse()
            patches['upba'] = sf.edge_curves(btm, right, top, left)
            patches.connect(
                ('upba', 3, 'ba', 4), ('upba', 1, 'up', 2),
                ('dnba', 4, 'ba', 3), ('dnba', 1, 'dn', 2),
            )
            patches.boundary('wall', 'upba', 4)
            patches.boundary('outflow', 'upfr', 2)
            patches.boundary('wall', 'dnba', 3)
            patches.boundary('outflow', 'dnba', 2)
        else:
            patches.boundary('outflow', 'up', 2)
            patches.boundary('outflow', 'dn', 2)

        nel = int(ceil(log(1 - 1/dl * (1 - grad) * side) / log(grad)))
        for uk in {'up', 'upfr', 'upba'} & patches.keys():
            dk = 'dn' + uk[2:]
            patches[dk] = patches[uk] - (0, side + 2 * width, 0)
            geometric_refine(patches[uk], grad, nel - 1, direction='v')
            geometric_refine(patches[dk], grad, nel - 1, direction='v', reverse=True)

    else:
        patches.boundary('wall', 'ou', 2)
        patches.boundary('wall', 'od', 2)
        patches.boundary('wall', 'or', 4, dim=-2, add=vx_add)
        patches.boundary('wall', 'or', 2, dim=-2, add=vx_add)

        if 'fr' in patches:
            patches.boundary('wall', 'fr', 4)
            patches.boundary('wall', 'fr', 3)
            patches.boundary('wall', 'ol', 2, dim=-2, add=vx_add)
            patches.boundary('wall', 'ol', 4, dim=-2, add=vx_add)

        if 'ba' in patches:
            patches.boundary('wall', 'ba', 4)
            patches.boundary('wall', 'ba', 3)


    if height > 0.0:
        names = ['iu', 'il', 'id', 'ir', 'ou', 'ol', 'od', 'or',
                 'fr', 'ba', 'up', 'upba', 'upfr', 'dn', 'dnba', 'dnfr']
        for pname in names:
            if pname not in patches:
                continue
            patch = patches[pname]
            patch = extrude(patch, (0, 0, height))
            patch.raise_order(0, 0, 2)
            patch.refine(nel_height-1, direction='w')
            patches[pname] = patch
            patches.boundary('zup', pname, 6)
            patches.boundary('zdown', pname, 5)
            patches.connect((pname, 5, pname, 6, 'per'))

    patches.write(out, order=order)


if __name__ == '__main__':
    cylinder()
