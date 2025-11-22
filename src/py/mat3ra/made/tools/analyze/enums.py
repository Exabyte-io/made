POSSIBLE_TRANSFORMATION_MATRICES = [
    [[1, 0, 0], [0, 1, 0], [0, 0, 1]],  # Identity - no transformation
    # Direct swaps: new[i] = old[j]
    [[0, 1, 0], [1, 0, 0], [0, 0, 1]],
    [[0, 0, 1], [0, 1, 0], [1, 0, 0]],
    [[1, 0, 0], [0, 0, 1], [0, 1, 0]],
    # Swaps with sign flips
    [[0, 0, 1], [1, 0, 0], [0, -1, 0]],
    [[0, 0, -1], [1, 0, 0], [0, 1, 0]],
    [[0, 1, 0], [0, 0, 1], [-1, 0, 0]],
    [[0, -1, 0], [0, 0, 1], [1, 0, 0]],
    # Inverted swaps
    [[0, -1, 0], [-1, 0, 0], [0, 0, 1]],
    [[0, 0, -1], [0, 1, 0], [-1, 0, 0]],
    [[1, 0, 0], [0, 0, -1], [0, -1, 0]],
    # Rotations around x-axis: -90 and +90 degrees
    [[1, 0, 0], [0, 0, 1], [0, -1, 0]],  # -90° around x
    [[1, 0, 0], [0, 0, -1], [0, 1, 0]],  # +90° around x
    # Rotations around y-axis: -90 and +90 degrees
    [[0, 0, -1], [0, 1, 0], [1, 0, 0]],  # -90° around y
    [[0, 0, 1], [0, 1, 0], [-1, 0, 0]],  # +90° around y
    # Rotations around z-axis: -90 and +90 degrees
    [[0, 1, 0], [-1, 0, 0], [0, 0, 1]],  # -90° around z
    [[0, -1, 0], [1, 0, 0], [0, 0, 1]],  # +90° around z
    # 180° rotations around axes
    [[1, 0, 0], [0, -1, 0], [0, 0, -1]],  # 180° around x
    [[-1, 0, 0], [0, 1, 0], [0, 0, -1]],  # 180° around y
    [[-1, 0, 0], [0, -1, 0], [0, 0, 1]],  # 180° around z
    # Mirrors (reflections)
    [[1, 0, 0], [0, 1, 0], [0, 0, -1]],  # Mirror along z
    [[1, 0, 0], [0, -1, 0], [0, 0, 1]],  # Mirror along y
    [[-1, 0, 0], [0, 1, 0], [0, 0, 1]],  # Mirror along x
]
