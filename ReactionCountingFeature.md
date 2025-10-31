# Reaction Counting Per Cell Feature

## Overview

This feature adds comprehensive tracking of reactions that occur in each cell during binary collisions in WarpX. The system tracks:
- **Reaction counts**: Total number of reactions per cell
- **Reaction rates**: Reactions per unit volume (normalized by cell volume)
- **Reaction subtypes**: Support for distinguishing between different reaction types (e.g., multiple excitation/ionization levels in DSMC)

## Implementation Details

### Data Structure

- **MultiFab Arrays**: Each collision that produces particles maintains two arrays of MultiFab objects, one per refinement level:
  - `m_reaction_counts_per_cell`: Raw reaction counts
  - `m_reaction_rates_per_cell`: Reaction rates (counts per unit volume)
- **Multi-component Support**: Each MultiFab can have multiple components to track different reaction subtypes
- **Volume Normalization**: Rates are automatically normalized by cell volume, accounting for geometric factors in RZ, RCYLINDER, and RSPHERE geometries
- **Reset**: Counters are reset to zero at the start of each collision call

### Key Components

1. **Member Variables** (`BinaryCollision.H`):
   ```cpp
   CollisionType m_collision_type;
   amrex::Vector<std::string> m_reaction_names;
   std::map<int, int> m_mask_to_component;
   amrex::Vector<std::unique_ptr<amrex::MultiFab>> m_reaction_counts_per_cell;
   amrex::Vector<std::unique_ptr<amrex::MultiFab>> m_reaction_rates_per_cell;
   ```

2. **Initialization** (in `doCollisions`):
   - MultiFabs are created with the same BoxArray and DistributionMapping as the particle containers
   - N components per MultiFab where N = number of reaction subtypes
   - Zero ghost cells (cell-centered data)
   - Reaction names are extracted from collision configuration

3. **Counting** (in `doCollisionsWithinTile`):
   - After particle creation, the code loops through each cell
   - Reads `p_mask` values to determine which reactions occurred
   - For DSMC: mask values indicate specific scattering process types
   - For other collisions: mask is boolean (0 or non-zero)
   - Counts reactions by subtype using local accumulators
   - Computes rates by dividing counts by cell volume (including geometric factors)
   - Atomically updates MultiFabs using `Gpu::Atomic::Add`

### Accessor Methods

```cpp
// Get the collision type
CollisionType getCollisionType() const;

// Get the reaction counts MultiFab for a specific level
amrex::MultiFab* getReactionCountsAtLevel(int lev) const;

// Get the reaction rates MultiFab for a specific level
amrex::MultiFab* getReactionRatesAtLevel(int lev) const;

// Get the total number of reactions across all cells at a level
// comp = -1: sum all reaction types, comp >= 0: specific reaction type
long getTotalReactionCountAtLevel(int lev, int comp = -1) const;

// Get reaction names/labels
const amrex::Vector<std::string>& getReactionNames() const;

// Get number of reaction subtypes being tracked
int getNumReactionTypes() const;
```

## Usage Example

```cpp
// In your collision handler or diagnostic code:
BinaryCollision* collision = /* get your collision object */;

// Get the collision type
CollisionType type = collision->getCollisionType();

// Get reaction names
const auto& reaction_names = collision->getReactionNames();
int n_reaction_types = collision->getNumReactionTypes();

std::cout << "Tracking " << n_reaction_types << " reaction types:" << std::endl;
for (int i = 0; i < n_reaction_types; ++i) {
    std::cout << "  [" << i << "] " << reaction_names[i] << std::endl;
}

// Get reaction counts for level 0
amrex::MultiFab* reaction_counts = collision->getReactionCountsAtLevel(0);
amrex::MultiFab* reaction_rates = collision->getReactionRatesAtLevel(0);

if (reaction_counts && reaction_rates) {
    // Access per-cell reaction data
    for (amrex::MFIter mfi(*reaction_counts); mfi.isValid(); ++mfi) {
        const amrex::FArrayBox& counts_fab = (*reaction_counts)[mfi];
        const amrex::FArrayBox& rates_fab = (*reaction_rates)[mfi];
        const amrex::Box& box = mfi.validbox();

        // Process data for each reaction type
        for (int comp = 0; comp < n_reaction_types; ++comp) {
            const amrex::Real* counts = counts_fab.dataPtr(comp);
            const amrex::Real* rates = rates_fab.dataPtr(comp);

            // Access data for each cell in the box
            // counts[i] = number of reactions of type 'comp' in cell i
            // rates[i] = reactions per unit volume in cell i
        }
    }

    // Get total reactions across all cells
    long total_reactions_all_types = collision->getTotalReactionCountAtLevel(0);

    // Get total for a specific reaction type (e.g., component 0)
    long total_reactions_type0 = collision->getTotalReactionCountAtLevel(0, 0);

    std::cout << "Total reactions (all types): " << total_reactions_all_types << std::endl;
    std::cout << "Total " << reaction_names[0] << ": " << total_reactions_type0 << std::endl;
}
```

## Supported Collision Types

This feature works with any binary collision that produces particles, including:
- Nuclear fusion reactions (D-T, D-D, D-He3, p-B11)
- DSMC collisions with ionization
- Compton scattering
- Breit-Wheeler pair production
- Any collision with `product_species` defined

## Performance Considerations

- **GPU Compatible**: All counting operations use GPU-compatible atomic operations
- **Thread Safe**: Atomic adds ensure correct counting in parallel execution
- **Minimal Overhead**: Counting is done once per cell after collision processing using local accumulators
- **Memory**: Adds N×2 `Real` values per cell (N components × 2 MultiFabs: counts + rates)
- **Scalability**: Supports up to 10 reaction subtypes per collision (compile-time limit for stack allocation)
- **Volume Normalization**: Geometric factors (RZ, RCYLINDER, RSPHERE) computed once per cell

## Integration with Diagnostics

The reaction count MultiFabs can be integrated with WarpX's diagnostic system by:
1. Adding a diagnostic output that writes the MultiFab to file
2. Computing derived quantities (reaction rate per volume, fusion power, etc.)
3. Tracking time evolution of reactions in specific regions

## Future Enhancements

Possible extensions to this feature:

### Already Implemented ✓
- ✓ Track reactions by reaction subtype (multi-component MultiFabs)
- ✓ Reaction rate diagnostics (reactions per unit volume)
- ✓ Volume normalization with geometric factors

### Potential Future Work
- **Enhanced DSMC Support**: Dynamically extract scattering process names from DSMC configuration instead of generic labels
- **Process-Specific Mapping**: Runtime detection of which DSMC scattering processes are actually used, with proper mask-to-component mapping
- **Separate Product Species Tracking**: Individual counters for each product species (e.g., track electrons vs ions separately in ionization)
- **Temporal Integration**: Time-integrated reaction counts over multiple timesteps
- **Spatial Reduction Operations**: Built-in sum over regions, max/min finding, centroid calculations
- **Energy Spectrum**: Track energy distribution of reactions
- **Expand Component Limit**: Increase the 10-component limit for collisions with many scattering processes

### DSMC Integration Note

The current implementation provides basic DSMC support by counting all reactions together. For full per-process tracking in DSMC:
1. Add a public accessor method to `DSMCFunc` to expose `m_scattering_processes`
2. Extract process names and types in `initializeReactionInfo()`
3. Create a proper mapping from `ScatteringProcessType` enum values to component indices
4. Update the counting loop to use this mapping

This would allow separate tracking of elastic, ionization, excitation, and other scattering processes.
