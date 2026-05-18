# -*- coding: utf-8 -*-
"""
Fragment — a reactive sub-region of a Molspace with an initiator atom, a set
of edge atoms, the NTAs of each edge's severed (external) neighbors, and an
optional set of atom ids excluded from charge remapping.

For LAMMPS ``fix bond/react`` template matching, edge atoms must each have
exactly one bond into the fragment interior; the constructor warns if any edge
violates this. The fragment ``mol`` carries the original atom ids so the
caller's ``breaks`` / ``makes`` / ``changed_typ`` lists keep referring to the
ids they already know.

Construction modes (each accepts a :class:`~mooonpy.Molspace` *or* another
:class:`Fragment` as the source — when a Fragment is passed, its stored
``edge_external_ntas`` are inherited for any edge that survives into the new
fragment):

* :meth:`Fragment.from_initiator_and_edges` — caller supplies an initiator and
  a list of edge ids; BFS walks outward from the initiator and stops at the
  edges.
* :meth:`Fragment.from_interior` — caller supplies the interior atom set;
  edges are derived as the boundary atoms.
* :meth:`Fragment.from_initiator_steps` — caller supplies an initiator and an
  integer ``n_steps``; BFS walks exactly ``n_steps`` outward.
* :meth:`Fragment.assembled` — wrap an already-built Molspace (e.g. a combined
  reaction template) without slicing; edges and external NTAs are supplied
  explicitly. ``initiator`` is ``None`` for this whole-template form.
"""
import warnings

from mooonpy.molspace.graph_theory.interface import (
    find_reachable_neighbors, steps_from_initiator,
)


def _resolve_input(src):
    """
    Accept a Molspace or a Fragment as a construction source.

    :return: ``(mol, inherited_edge_ntas)`` where ``mol`` is the Molspace to
             walk and ``inherited_edge_ntas`` is the source Fragment's stored
             ``edge_external_ntas`` (empty dict when ``src`` is a Molspace).
    """
    if isinstance(src, Fragment):
        return src.mol, dict(src.edge_external_ntas)
    return src, {}


class Fragment:
    """
    Reactive sub-region of a Molspace.

    :param mol: Molspace containing the fragment. Sliced (on a copy) down to
                ``interior_ids`` without renumbering, so original atom ids
                survive.
    :type mol: Molspace
    :param initiator: Atom id of the initiator inside the fragment, or
                      ``None`` for an assembled whole-template Fragment.
    :type initiator: int or None
    :param interior_ids: Atom ids that belong to the fragment.
    :type interior_ids: iterable[int]
    :param edges: Atom ids that bound the fragment.
    :type edges: iterable[int]
    :param charge_edge_depth: If ``>= 0``, atoms within this many bonds of any
                              edge are added to :attr:`no_remap_ids`. ``-1``
                              disables.
    :type charge_edge_depth: int
    :param inherited_edge_ntas: When constructing from another Fragment, the
        parent's ``edge_external_ntas``. Used to keep the severed-neighbor NTA
        list for any edge whose external neighbors were already pruned from
        the immediate ``mol``.
    :type inherited_edge_ntas: dict[int, list[str]] or None
    :param edge_external_ntas: Explicit ``{edge_id: [nta, ...]}`` override
        (used by :meth:`assembled`). When given, the parent-graph computation
        is skipped entirely.
    :type edge_external_ntas: dict[int, list[str]] or None
    :param no_remap_ids: Explicit set of charge-preserved ids (used by
        :meth:`assembled`). When given, ``charge_edge_depth`` is ignored.
    :type no_remap_ids: set[int] or None

    Attributes
    ----------
    mol : Molspace
        Sliced Molspace with original atom ids.
    initiator : int or None
    edges : set[int]
    edge_external_ntas : dict[int, list[str]]
        ``{edge_id: [nta_of_severed_neighbor, ...]}``.
    no_remap_ids : set[int]
        Atom ids whose charges should *not* be remapped during the reaction.
    charge_edge_depth : int
    """

    def __init__(self, mol, initiator, interior_ids, edges,
                 charge_edge_depth=-1, inherited_edge_ntas=None,
                 edge_external_ntas=None, no_remap_ids=None):
        interior = set(interior_ids)
        edge_set = set(edges)

        if initiator is not None and initiator not in interior:
            raise ValueError(
                f"initiator={initiator} is not in the interior set"
            )
        if not edge_set.issubset(interior):
            stray = edge_set - interior
            raise ValueError(
                f"edges {sorted(stray)} are not in the interior set"
            )

        if edge_external_ntas is not None:
            # Assembled form: caller already knows the severed-neighbor NTAs.
            self.edge_external_ntas = dict(edge_external_ntas)
        else:
            # Compute each edge's external-neighbor NTAs from the parent's
            # bond graph *before* slicing. When an edge's externals were
            # already pruned (the parent was itself a sliced Fragment.mol),
            # fall back to the inherited list so the info isn't lost.
            inherited = inherited_edge_ntas or {}
            parent_graph = mol.generate_graph()
            self.edge_external_ntas = {}
            for edge_id in edge_set:
                externals = [nbr for nbr in parent_graph[edge_id]
                             if nbr not in interior]
                computed = [mol.atoms[nbr].comment for nbr in externals]
                if not computed and edge_id in inherited:
                    self.edge_external_ntas[edge_id] = list(inherited[edge_id])
                else:
                    self.edge_external_ntas[edge_id] = computed

        self.mol = mol.copy()
        # Slice the working copy without renumbering so original ids persist.
        self.mol.slice(interior, renumber=False)

        self.initiator = initiator
        self.edges = edge_set
        self.charge_edge_depth = int(charge_edge_depth)

        # The single-interior-bond rule is about a *sliced* fragment's
        # boundary anchors. An assembled whole-template (initiator is None)
        # keeps every atom, so its edges legitimately have many in-template
        # bonds — skip the check there.
        if initiator is not None:
            self._validate_edge_topology()
        if no_remap_ids is not None:
            self.no_remap_ids = set(no_remap_ids)
        else:
            self.no_remap_ids = self._compute_no_remap_ids()

    # ------------------------------------------------------------------
    # Constructors
    # ------------------------------------------------------------------

    @classmethod
    def from_initiator_and_edges(cls, src, initiator, edges,
                                 max_steps=10, charge_edge_depth=-1):
        """
        BFS outward from ``initiator``; stop at any atom in ``edges``.

        :param src: Molspace or Fragment to walk.
        :param initiator: Atom id to start from.
        :param edges: Atom ids that bound the fragment.
        :param max_steps: Safety cutoff on BFS depth. A warning is emitted if
                          BFS hits this limit before fully enclosing the
                          fragment (caller probably forgot an edge id).
        :param charge_edge_depth: See :class:`Fragment`.
        """
        mol, inherited = _resolve_input(src)
        edge_set = set(edges)
        if initiator in edge_set:
            raise ValueError("initiator cannot also be in edges")
        if initiator not in mol.atoms:
            raise KeyError(f"initiator={initiator} not in mol.atoms")

        graph = mol.generate_graph()
        interior = {initiator}
        frontier = {initiator}
        steps = 0
        while frontier and steps < max_steps:
            next_frontier = set()
            for id_ in frontier:
                if id_ in edge_set:
                    continue  # don't walk past an edge
                for nbr in graph[id_]:
                    if nbr in interior:
                        continue
                    interior.add(nbr)
                    next_frontier.add(nbr)
            frontier = next_frontier
            steps += 1

        if frontier:
            warnings.warn(
                f"Fragment.from_initiator_and_edges: BFS hit max_steps={max_steps} "
                f"without enclosing the fragment from initiator {initiator}; "
                f"some edge ids may be missing.",
                stacklevel=2,
            )

        reached_edges = edge_set & interior
        return cls(mol, initiator, interior, reached_edges,
                   charge_edge_depth=charge_edge_depth,
                   inherited_edge_ntas=inherited)

    @classmethod
    def from_interior(cls, src, initiator, interior_ids, charge_edge_depth=-1):
        """
        Use ``interior_ids`` as-is; derive edges as the boundary atoms
        (interior atoms with at least one neighbor outside).
        """
        mol, inherited = _resolve_input(src)
        interior = set(interior_ids)
        if initiator not in interior:
            raise ValueError(
                f"initiator={initiator} is not in interior_ids"
            )
        graph = mol.generate_graph()
        edges = {id_ for id_ in interior
                 if any(nbr not in interior for nbr in graph[id_])}
        return cls(mol, initiator, interior, edges,
                   charge_edge_depth=charge_edge_depth,
                   inherited_edge_ntas=inherited)

    @classmethod
    def from_initiator_steps(cls, src, initiator, n_steps, charge_edge_depth=-1):
        """
        Walk exactly ``n_steps`` bonds outward from ``initiator``. The
        outermost shell becomes the edge set.
        """
        mol, inherited = _resolve_input(src)
        if initiator not in mol.atoms:
            raise KeyError(f"initiator={initiator} not in mol.atoms")
        graph = mol.generate_graph()
        interior, boundary = steps_from_initiator(graph, initiator, n_steps)
        return cls(mol, initiator, interior, boundary,
                   charge_edge_depth=charge_edge_depth,
                   inherited_edge_ntas=inherited)

    @classmethod
    def assembled(cls, mol, edges, edge_external_ntas, no_remap_ids=None):
        """
        Wrap an already-built Molspace (e.g. a combined reaction template)
        as a Fragment without slicing. ``initiator`` is ``None``.

        :param mol: The whole-template Molspace.
        :param edges: Edge atom ids (in ``mol``'s numbering).
        :param edge_external_ntas: ``{edge_id: [nta, ...]}`` severed-neighbor
                                   NTAs, already in ``mol``'s numbering.
        :param no_remap_ids: Charge-preserved ids, or ``None``.
        """
        return cls(mol, None, set(mol.atoms), edges,
                   edge_external_ntas=edge_external_ntas,
                   no_remap_ids=no_remap_ids)

    # ------------------------------------------------------------------
    # Internals
    # ------------------------------------------------------------------

    def _validate_edge_topology(self):
        """
        Each edge atom must have exactly one bond into the fragment interior.
        ``fix bond/react`` matches templates by graph isomorphism with the
        edges as anchors; multi-bond edges break that match.
        """
        graph = self.mol.generate_graph()
        interior = set(self.mol.atoms)
        for edge_id in self.edges:
            n_interior = sum(1 for nbr in graph[edge_id] if nbr in interior)
            if n_interior != 1:
                warnings.warn(
                    f"Fragment edge atom {edge_id} has {n_interior} bonds "
                    f"into the fragment interior; fix bond/react requires "
                    f"exactly 1 for template matching.",
                    stacklevel=3,
                )

    def _compute_no_remap_ids(self):
        if self.charge_edge_depth < 0:
            return set()
        graph = self.mol.generate_graph()
        keep = set()
        for edge_id in self.edges:
            keep.update(find_reachable_neighbors(graph, edge_id, self.charge_edge_depth))
        return keep
