# Reaction Counting Per Cell Feature

## Overview

This feature adds the capability to track the number of reactions that occur in each cell during binary collisions in WarpX. The reaction counts are stored per refinement level and can be accessed for diagnostic purposes or further analysis.

## Implementation Details

### Data Structure

- **MultiFab**: Each collision that produces particles (i.e., has `product_species` defined) now maintains a `Vector` of `MultiFab` objects, one per refinement level
- **Storage**: One `Real` value per cell representing the total number of reactions in that cell during the current timestep
- **Reset**: Counters are reset to zero at the start of each collision call

### Key Components

1. **Member Variables** (`BinaryCollision.H`):
   ```cpp
   CollisionType m_collision_type;
   amrex::Vector<std::unique_ptr<amrex::MultiFab>> m_reaction_counts_per_cell;
   ```

2. **Initialization** (in `doCollisions`):
   - MultiFabs are created with the same BoxArray and DistributionMapping as the particle containers
   - One component per MultiFab (reaction count)
   - Zero ghost cells (cell-centered data)

3. **Counting** (in `doCollisionsWithinTile`):
   - After particle creation, the code loops through each cell
   - Counts the number of `true` values in the `p_mask` array for that cell's pairs
   - Atomically adds the count to the MultiFab using `Gpu::Atomic::Add`

### Accessor Methods

```cpp
// Get the collision type
CollisionType getCollisionType() const;

// Get the reaction counts MultiFab for a specific level
amrex::MultiFab* getReactionCountsAtLevel(int lev) const;

// Get the total number of reactions across all cells at a level
long getTotalReactionCountAtLevel(int lev) const;
```

## Usage Example

```cpp
// In your collision handler or diagnostic code:
BinaryCollision* collision = /* get your collision object */;

// Get the collision type
CollisionType type = collision->getCollisionType();

// Get reaction counts for level 0
amrex::MultiFab* reaction_counts = collision->getReactionCountsAtLevel(0);

if (reaction_counts) {
    // Access per-cell reaction counts
    for (amrex::MFIter mfi(*reaction_counts); mfi.isValid(); ++mfi) {
        const amrex::FArrayBox& fab = (*reaction_counts)[mfi];
        const amrex::Box& box = mfi.validbox();

        // Process reaction counts for this box
        // fab contains the reaction count for each cell
    }

    // Or get total reactions across all cells
    long total_reactions = collision->getTotalReactionCountAtLevel(0);
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
- **Minimal Overhead**: Counting is done once per cell after collision processing
- **Memory**: Adds one `Real` value per cell per collision type at each level

## Integration with Diagnostics

The reaction count MultiFabs can be integrated with WarpX's diagnostic system by:
1. Adding a diagnostic output that writes the MultiFab to file
2. Computing derived quantities (reaction rate per volume, fusion power, etc.)
3. Tracking time evolution of reactions in specific regions

## Future Enhancements

Possible extensions to this feature:
- Track reactions by reaction subtype (for DSMC with multiple scattering processes)
- Separate counters for different product species
- Reaction rate diagnostics (reactions per unit time per unit volume)
- Spatial reduction operations (sum over regions, max/min finding)
