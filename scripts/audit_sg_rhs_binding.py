"""Check complete element-edge support and its signed continuity RHS ledger.

The caller supplies equally normalized terms from one immutable RHS snapshot.
No missing support is inferred from a node total; zero-weight edges count.
"""
from fractions import Fraction


def finite(value):
    if isinstance(value, bool):
        raise ValueError("Numeric operands required")
    try:
        return Fraction(value)
    except (ValueError, TypeError, OverflowError, ZeroDivisionError) as exc:
        raise ValueError("Finite operands required") from exc


def bind_row(node, cells, edges, source_additions, solver_rhs, snapshot):
    """Return exact sums after checking full support and snapshot identity.

    cells maps cell IDs to three distinct vertex IDs. Edges contain cell,
    nodes (directed pair), product, factor (+/-1), and snapshot. Source
    additions contain value and snapshot. Solver RHS is the negative of the
    assembled continuity residual, as explicitly chosen by this interface.
    """
    if not isinstance(snapshot,str) or not snapshot:
        raise ValueError("Explicit RHS snapshot identity required")
    expected=set()
    for cell, nodes in cells.items():
        if len(nodes)!=3 or len(set(nodes))!=3:
            raise ValueError("Distinct Tri3 vertices required")
        if node in nodes:
            expected.update((str(cell),min(node,n),max(node,n))for n in nodes if n!=node)
    if not expected:
        raise ValueError("Target has no mesh support")
    observed=set();terms=[]
    for edge in edges:
        if edge['snapshot']!=snapshot:
            raise ValueError("Mixed RHS snapshots")
        a,b=edge['nodes']
        if a==b:raise ValueError("Degenerate edge")
        factor=finite(edge['factor'])
        if factor not in [-1,1]:raise ValueError("Signed assembly factor required")
        product=finite(edge['product'])
        if node not in (a,b):continue
        key=(str(edge['cell']),min(a,b),max(a,b))
        if key in observed or key not in expected:
            raise ValueError("Duplicate or foreign support")
        observed.add(key);terms.append(factor*product*(1 if a==node else -1))
    if observed!=expected:
        raise ValueError("Incomplete support, including zero terms")
    sources=[]
    for source in source_additions:
        if source['snapshot']!=snapshot:raise ValueError("Mixed source/RHS snapshot")
        sources.append(finite(source['value']))
    transport=sum(terms,Fraction());source=sum(sources,Fraction())
    absolute=sum(map(abs,terms),Fraction())+sum(map(abs,sources),Fraction())
    difference=abs(transport+source+finite(solver_rhs))
    return {'transport':transport,'source':source,'absolute':absolute,
            'difference':difference,'zero_budget':absolute==0,
            'error':difference/absolute if absolute else (Fraction()if not difference else None),
            'support_terms':len(terms)}
