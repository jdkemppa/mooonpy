# -*- coding: utf-8 -*-
from itertools import product
from collections import defaultdict
from typing import Dict, List, Tuple, Union, Any, Optional, Iterator, Iterable
from dataclasses import dataclass

from mooonpy.tools.file_utils import Path


@dataclass
class DimensionEntry:
    """Represents a single entry in a dimension with its parsing information."""
    cart_key: str
    condition: Optional[str]
    name: str
    iterable: List[Any]


@dataclass
class DimensionInfo:
    """Information about a parsed dimension including base and conditional entries."""
    base: DimensionEntry
    conditionals: List[DimensionEntry]
    parallels: List[DimensionEntry]
    size: int


@dataclass
class IndexInfo:
    """Index information for a generated permutation."""
    global_idx: int
    dimensions: Dict[str, int]
    combined: Optional[Dict[Tuple[str, ...], int]] = None


class CartesianResult(dict):
    """Result dictionary with index information as attributes."""

    def __init__(self, values: Dict[str, Any],
                 global_idx: Optional[int] = None,
                 dimension_indices: Optional[Dict[str, int]] = None,
                 combined_indices: Optional[Dict[Tuple[str, ...], int]] = None):
        """
        Initialize CartesianResult.

        Args:
            values: The main result dictionary {name: value}
            global_idx: Global permutation index
            dimension_indices: Per-dimension indices
            combined_indices: Combined dimension indices
        """
        super().__init__(values)
        self.global_idx = global_idx
        self.dimension_indices = dimension_indices or {}
        self.combined_indices = combined_indices or {}

    def string_format(self, strings: Union[str | Path, List[str | Path], Dict[Any, str | Path]]) -> Union[
        str | Path, List[str | Path], Dict[Any, str | Path]]:
        """
        Format strings with {var} using the CartesianResult dict.
        Strings may be a single string, list of strings, or dict of strings.
        Path objects are preserved as paths.
        """
        if isinstance(strings, str):
            if isinstance(strings, Path):
                return Path(strings.format(**self))
            else:
                return strings.format(**self)

        elif isinstance(strings, list):
            out = []
            for string in strings:
                formatted = string.format(**self)
                if isinstance(string, Path):
                    out.append(Path(formatted))
                else:
                    out.append(formatted)
            return out

        elif isinstance(strings, dict):
            out = {}
            for key, string in strings.items():
                formatted = string.format(**self)
                if isinstance(string, Path):
                    out[key] = Path(formatted)
                else:
                    out[key] = formatted
            return out
        else:
            raise TypeError('strings must be str or Path, and datastructure must be str, list or dict')


class Cartesian:
    """Enhanced cartesian product generator with index tracking and conditional dimensions."""

    def __init__(self, cartesian_sets: Dict[Union[str, int], Iterable[Any]], delim: str = '.'):
        """
        Initialize Cartesian generator.

        Args:
            cartesian_sets: Dictionary of cart_key: iterable pairs
            delim: Delimiter for parsing dimension.condition.name syntax
        """
        self.cartesian_sets = cartesian_sets
        self.delim = delim
        self.dimensions: Dict[str, DimensionInfo] = {}
        self.dimension_sizes: Dict[str, int] = {}
        self._parse_cartesian_sets()

    def _parse_cartesian_sets(self) -> None:
        """Parse cartesian_sets into organized dimension structure."""
        # Group by dimensions and parse conditions
        dim_groups: Dict[str, List[DimensionEntry]] = defaultdict(list)

        for cart_key, iterable in self.cartesian_sets.items():
            # Ensure iterable
            if isinstance(iterable, str):
                pass  # keep as a str
            elif not hasattr(iterable, '__iter__'):
                if not isinstance(iterable, (list, tuple)):
                    iterable = [iterable]
            else:
                iterable = list(iterable)  # Convert to list for consistent indexing

            # Parse cart_key
            parts = str(cart_key).split(self.delim)
            if len(parts) == 1:
                # Simple case: dimension and name are the same
                dimension = parts[0]
                condition = None
                name = parts[0]
            else:
                dimension = parts[0]
                condition = parts[1] if len(parts) > 2 else None
                name = parts[-1]

            entry = DimensionEntry(
                cart_key=str(cart_key),
                condition=condition,
                name=name,
                iterable=iterable
            )
            dim_groups[dimension].append(entry)

        # Sort dimensions for consistent ordering
        sorted_dims = sorted(dim_groups.keys(), key=str)

        # Build dimension structure
        for dimension in sorted_dims:
            entries = dim_groups[dimension]

            # Find base entry (no condition) and conditional entries
            base_entry: Optional[DimensionEntry] = None
            conditional_entries: List[DimensionEntry] = []
            parallel_entries: List[DimensionEntry] = []
            for entry in entries:
                if entry.condition is not None:
                    conditional_entries.append(entry)
                elif entry.name != entry.cart_key:
                    parallel_entries.append(entry)
                else:
                    base_entry = entry

            if base_entry is None:
                raise ValueError(f"Dimension '{dimension}' must have at least one unconditional entry")

            dim_info = DimensionInfo(
                base=base_entry,
                conditionals=conditional_entries,
                parallels=parallel_entries,
                size=len(base_entry.iterable)
            )

            self.dimensions[dimension] = dim_info
            self.dimension_sizes[dimension] = len(base_entry.iterable)

    def __iter__(self, return_combined_indices: Optional[List[Tuple[str, ...]]] = None) -> Iterator[CartesianResult]:
        """Call generator function. """
        return self.generate(return_combined_indices)

    def generate(self,
                 return_combined_indices: Optional[List[Tuple[str, ...]]] = None) -> Iterator[CartesianResult]:
        """
        Generate cartesian product permutations.

        Args:
            return_combined_indices: List of dimension name tuples for combined indexing
                                   e.g., [('R', 'dir')] to get combined R+dir index

        Yields:
            CartesianResult containing {name: value} pairs with index attributes populated
        """
        # Get dimension names in sorted order
        dim_names = sorted(self.dimensions.keys(), key=str)

        # Validate combined index requests
        if return_combined_indices:
            for dim_combo in return_combined_indices:
                if not all(dim in self.dimensions for dim in dim_combo):
                    raise ValueError(f"Unknown dimensions in combination: {dim_combo}")

        # Generate all combinations of dimension indices
        dim_ranges = [range(self.dimensions[dim].size) for dim in dim_names]

        global_idx = 0
        for dim_indices in product(*dim_ranges):
            # Build dimension index mapping
            dim_index_map = dict(zip(dim_names, dim_indices))

            # Build result dictionary
            result: Dict[str, Any] = {}

            # Track if any conditionals were triggered
            conditional_results: List[Dict[str, Any]] = []

            for dim_name, dim_idx in zip(dim_names, dim_indices):
                dim_info = self.dimensions[dim_name]

                # Add base values
                base_entry = dim_info.base
                result[base_entry.name] = base_entry.iterable[dim_idx]

                # Add parallel values
                for para in self.dimensions[dim_name].parallels:
                    result[para.name] = para.iterable[dim_idx]

                # Check for conditional values
                triggered_conditionals: List[DimensionEntry] = []
                for cond_entry in dim_info.conditionals:
                    condition = cond_entry.condition
                    if condition is not None and str(dim_idx) == str(condition):
                        triggered_conditionals.append(cond_entry)

                # If conditionals are triggered, generate combinations
                if triggered_conditionals:
                    if not conditional_results:
                        conditional_results = [{}]

                    new_conditional_results: List[Dict[str, Any]] = []
                    for cond_entry in triggered_conditionals:
                        for existing_result in conditional_results:
                            for cond_value in cond_entry.iterable:
                                new_result = existing_result.copy()
                                new_result[cond_entry.name] = cond_value
                                new_conditional_results.append(new_result)
                    conditional_results = new_conditional_results

            # Yield results
            if conditional_results:
                # Yield each conditional combination
                for cond_result in conditional_results:
                    final_result = {**result, **cond_result}
                    yield self._create_result(
                        final_result,
                        global_idx,
                        dim_index_map,
                        return_combined_indices
                    )
                    global_idx += 1
            else:
                # Yield main result (no conditionals triggered)
                yield self._create_result(
                    result,
                    global_idx,
                    dim_index_map,
                    return_combined_indices
                )
                global_idx += 1

    def _create_result(self,
                       result: Dict[str, Any],
                       global_idx: int,
                       dim_index_map: Dict[str, int],
                       return_combined_indices: Optional[List[Tuple[str, ...]]]) -> CartesianResult:
        """Create CartesianResult with index information."""
        # Build combined indices if requested
        combined_indices: Optional[Dict[Tuple[str, ...], int]] = None
        if return_combined_indices:
            combined_indices = {}
            for dim_combo in return_combined_indices:
                # Calculate combined index
                combined_idx = 0
                multiplier = 1

                # Process dimensions in reverse order for proper indexing
                for dim in reversed(dim_combo):
                    combined_idx += dim_index_map[dim] * multiplier
                    multiplier *= self.dimension_sizes[dim]

                combined_indices[dim_combo] = combined_idx

        # Always create result with index information
        return CartesianResult(
            values=result,
            global_idx=global_idx,
            dimension_indices=dim_index_map.copy(),
            combined_indices=combined_indices
        )

    def get_combined_index(self, dimension_combo: Tuple[str, ...], value_dict: Dict[str, Any]) -> int:
        """
        Calculate combined index for given dimensions from a value dictionary.

        Args:
            dimension_combo: Tuple of dimension names
            value_dict: Dictionary of {name: value} pairs

        Returns:
            Combined index for the specified dimensions
        """
        if not all(dim in self.dimensions for dim in dimension_combo):
            raise ValueError(f"Unknown dimensions: {dimension_combo}")

        combined_idx = 0
        multiplier = 1

        for dim in reversed(dimension_combo):
            dim_info = self.dimensions[dim]
            base_entry = dim_info.base

            # Find index of value in dimension
            try:
                dim_idx = list(base_entry.iterable).index(value_dict[base_entry.name])
            except (ValueError, KeyError):
                raise ValueError(f"Value not found for dimension {dim}")

            combined_idx += dim_idx * multiplier
            multiplier *= self.dimension_sizes[dim]

        return combined_idx

    def generate_list(self, **kwargs) -> List[CartesianResult]:
        """Generate full list of permutations (convenience method)."""
        return list(self.generate(**kwargs))

    def __len__(self) -> int:
        """Compute Length By generating all results"""
        return len(self.generate_list())

    def __repr__(self) -> str:
        """Get information about parsed dimensions."""
        info = f"Cartesian Product with {self.dimension_sizes} dimension sizes"
        return info
