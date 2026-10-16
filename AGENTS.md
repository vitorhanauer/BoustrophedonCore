# Overview

Implement the Boustrophedon algorithm in C to traverse a polygon and return an array of points. BoustrophedonCore provides the core logic for wrappers in other programming languages and must also be callable directly from C.

The algorithm has four stages:

1. Divide the input area into smaller connected polygons.
2. Treat those polygons as graph nodes and find a short route between them.
3. Generate a coverage path inside each node.
4. Join the node paths, connecting the end of one to the beginning of the next.

# Entry point

Provide a single entry point for callers. Keep implementation responsibilities in separate functions. Conceptual signature:

```text
generateWaypoints(Start Point, End Point, Area Vertices, Obstacle Vertices, Coverage Width, Overlap);
```

- **Start point:** two coordinates, such as `[10, 10]` or `{x: 10, y: 10}`, where the autonomous agent starts the route.
- **End point:** two coordinates where the agent finishes the route.
- **Area vertices:** a sequence of coordinate pairs describing the polygon to cover, such as `[[10, 10], [10, 9], ...]`.
- **Obstacle vertices:** a sequence of polygons, each represented by its vertices, such as `[[[10, 10], [10, 9], ...], ...]`.
- **Coverage width:** the full sensing or working width, including both the left and right sides of the agent.
- **Overlap:** the percentage by which adjacent coverage swaths overlap.

Lane spacing is `width * (1 - overlapPercent / 100)`. The API accepts overlap in `[0, 100)`; exactly 100% would produce zero spacing.

# Coding practices

- Use Conventional Commits and commit each completed functional unit of work.
- Use readable names for variables, functions and other identifiers.
- Use TDD: write a failing functional test before implementing behavior.

# Resource constraints

The target test configuration is 512 MB of RAM and a 1 GHz CPU. Use the following square sizes and coverage widths:

| Square size | Coverage widths |
|---|---|
| 10 × 10 | 1, 5 |
| 100 × 100 | 1, 5, 10, 50 |
| 1,000 × 1,000 | 1, 5, 10, 50 |
| 10,000 × 10,000 | 1, 5, 10, 50 |
