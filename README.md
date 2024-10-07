# micos - Minimal Constraint Solver

micos (Minimalistic Constraint Solver) is a minimallistic accurate constraint-based physics solver designed for 2D simulations. 

![micos Demo](images/examples_micos.gif)

[Test Systems](https://jaimelopezgarcia.github.io/micos.js/dist/index.html#)


### Example usage

```html
<script src="https://cdn.statically.io/gh/jaimelopezgarcia/micos.js@master/dist/micos.bundle.js"></script>
```

```js
let PLAYER = new micos.Player(SVG);

let stateDoublePendulum = {
  "xs": [
    [0, 1],
    [0, 0]
  ],
  "vs": [
    [0, 0],
    [1, 0]
  ],
  "masses": [1, 1],
  "constraints_distance": [
    [0, 1, 1]
  ],
  "constraints_pin": [
    [0, [0, 2], 1]
  ],
  "gravity": [
    1,
    [0, -1]
  ],
  "time": 0,
  "damping": [
    [0, 0.01],
    [1, 0.01]
  ]
};

let options = { dt: 0.005 };

PLAYER.play(stateDoublePendulum, options);
```
### Mathematical Problem

The solver aims to minimize the difference between the unconstrained and constrained accelerations, following the **principle of least constraint**. Mathematically, this is expressed as the optimization problem:

    minimize || a_const - a_unconst ||^2_(M^-1)

subject to the following constraints:

- **Equality Constraints**: C(a_const) = 0
- **Inequality Constraints**: C_contact(a_const) <= 0

Here, `a_unconst` represents the unconstrained accelerations, `M` is the mass matrix, and `a_const` represents the constrained accelerations. The solver ensures that all constraints are respected while minimizing the deviation from the unconstrained solution.


### Solver Approach

- **Active Set Method**: Inequality constraints are handled using an active set method, which iteratively determines the active (binding) and inactive (non-binding) constraints.
- **Direct Solver for Constraints**: A direct LU method to resolve constraints equation.

### Numerical Integration

The solver uses a **semi-implicit Euler** method for time integration.


And relies on "pseudo-barrier" **Baumgarte stabilization**, to absorb inellastic collision energy and correct drift.
