"""Defect-one block tilings for complete degree lists and ordered prefixes."""


def defect_tiles(counts, p, *, last=None, slots=0):
    """Can the closed runs be tiled and their outstanding reservations fit?

    counts maps denominators to multiplicities. With suffix slots available,
    the run at last is still open. Only denominators prime to p have defect
    one. At 3 and 5, choose the least number of mixed blocks at each closed
    run. A surplus of p mixed blocks can be replaced by equal blocks at
    both denominators. At 5 the eight-label block is two four-label blocks
    as a multiset of degrees.
    """
    if p not in (2, 3, 5):
        raise ValueError("missing defect-one tiling rule")
    owed = {}
    for z, copies in sorted(counts.items()):
        if slots and z == last:
            break
        if z % p == 0:
            continue
        remaining = copies - owed.pop(z, 0)
        if remaining < 0:
            return False
        if p == 2:
            if remaining % 2:
                return False
        else:
            mixed = remaining % 3 if p == 3 else (3*remaining) % 5
            if p == 5 and 2*mixed > remaining:
                return False
            if mixed:
                owed[2*z] = 2*mixed
    if not slots:
        return not owed
    if any(z < last for z in owed):
        return False
    # Outstanding copies at the open run may be paid before using the suffix.
    prepaid = min(owed.get(last, 0), counts.get(last, 0))
    return sum(owed.values()) - prepaid <= slots
