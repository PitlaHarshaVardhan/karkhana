Code Structure:

The Python script defines a MobiusStrip class to model a Möbius strip and compute its geometric properties. The class is initialized with parameters R (radius), w (width), and n (resolution). 
It includes:
__init__: Sets up parameters and creates a meshgrid for parametric coordinates u and v.
parametric_equations: Generates a 3D mesh of (x, y, z) points using the 

parametric equations:
x(u,v) = (R + v·cos(u/2))·cos(u)
y(u,v) = (R + v·cos(u/2))·sin(u)
z(u,v) = v·sin(u/2)
compute_surface_area: Calculates the surface area numerically.
compute_edge_length: Computes the total length of the strip’s boundaries (v = ±w/2).
plot: Visualizes the strip in 3D using matplotlib and saves it as mobius_strip.png.

Surface Area Calculation:
The surface area is computed using the formula ∫∫ ||∂r/∂u × ∂r/∂v|| du dv, where ∂r/∂u and ∂r/∂v are partial derivatives of the position vector r(u,v) = (x(u,v), y(u,v), z(u,v)). I derived these derivatives analytically
∂r/∂u involves terms like -(R + v·cos(u/2))·sin(u) - (v/2)·sin(u/2)·cos(u) for x.
∂r/∂v includes terms like cos(u/2)·cos(u) for x. The cross product’s magnitude is calculated at each mesh point, and the double integral is approximated using NumPy’s trapz (trapezoidal rule) over u ∈ [0, 2π] and v ∈ [-w/2, w/2]. This numerical approach ensures accuracy for the non-orientable surface.


Challenges:
Analytical Derivatives: Deriving partial derivatives was tricky due to nested trigonometric functions, requiring careful chain rule application to avoid errors.
Numerical Stability: High resolution (n=100) was needed for accurate integration, but increased computation time. I balanced this by testing smaller n values.
3D Plotting: Ensuring the 3D plot was clear required adjusting transparency (alpha=0.8) and resolution. In VS Code, interactive display needed a compatible backend (e.g., TkAgg), which was addressed by modifying the plot method to include plt.show().
