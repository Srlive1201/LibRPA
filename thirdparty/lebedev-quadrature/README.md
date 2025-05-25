# Source

Here we describe the implementation details of the algorithm. 
A set of Lebedev quadrature points of any order is composed of sets of points which are invariant under the action of the [Octahedral group](https://www.wikiwand.com/en/Octahedral_symmetry).
Because of this invariance, all the points in a particular set share the same weight. 

One way to generate a given set is to start with a point and a weight (tabulated in the original papers) and then operate on it with each of the members of the Octahedral group.
Since there are 48 members of the Octahedral group, a set generated in this way can contain up to 48 different points.
However, certain points are unmoved by the action of certain elements of the Octahedral group.
For example, the point $(1, 0, 0)$ is unmoved by any rotations about the $x$-axis. 
Indeed, applying all the members of the Octahedral group to this point only gives 6 distinct points.

We wrap up this logic in the `GeneratorPoint` class, whose members are coordinates of a generating point $(a, b, c)$, the corresponding weight, and a rule for generating.  
In the case of $(1, 0, 0)$, ,the generating rule is `OctahedralPointGeneration::point_6` so named because it produces 6 points in all. 
All of the generating rules actually correspond to applying each of the members of the Octahedral group, but we streamline the process for points which are invariant under certain transformations.

To actually obtain a quadrature rule, one must use a `QuadratureOrder` enum which defines the order of quadrature to use.
This, along with several helper functions to tell when rules are available, and what precision they have (that is, what degree polynomial they can exactly integrate) are in the `quadrature_order.hpp` and `quadrature_order.inl` files.

Finally, to create all of the generator points which are used to create the quadrature points, there is a templated helper function `make_generator_points` which use the tabulated coordinates and weights to create the generator points.
Each of these generator points are then used to create sets of quadrature points according to their rules.
This is all done under the hood in the constructor of the `QuadraturePoints` object.
