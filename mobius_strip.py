import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

class MobiusStrip:
    def __init__(self, R=1.0, w=0.5, n=100):
        """
        Initialize Möbius strip with given parameters.
        
        Args:
            R (float): Radius from center to strip midline
            w (float): Width of the strip
            n (int): Resolution (number of points in each dimension)
        """
        self.R = R
        self.w = w
        self.n = n
        self.u = np.linspace(0, 2 * np.pi, n)
        self.v = np.linspace(-w/2, w/2, n)
        self.U, self.V = np.meshgrid(self.u, self.v)
        
    def parametric_equations(self):
        """
        Compute 3D coordinates (x, y, z) using parametric equations.
        
        Returns:
            tuple: (x, y, z) arrays of shape (n, n)
        """
        x = (self.R + self.V * np.cos(self.U / 2)) * np.cos(self.U)
        y = (self.R + self.V * np.cos(self.U / 2)) * np.sin(self.U)
        z = self.V * np.sin(self.U / 2)
        return x, y, z
    
    def compute_surface_area(self):
        """
        Compute surface area using numerical integration.
        Surface area = ∫∫ ||∂r/∂u × ∂r/∂v|| du dv
        """
        # Partial derivatives
        du = self.u[1] - self.u[0]
        dv = self.v[1] - self.v[0]
        
        # Parametric equations
        x = lambda u, v: (self.R + v * np.cos(u/2)) * np.cos(u)
        y = lambda u, v: (self.R + v * np.cos(u/2)) * np.sin(u)
        z = lambda u, v: v * np.sin(u/2)
        
        # Partial derivatives w.r.t u and v
        dx_du = lambda u, v: (-(self.R + v * np.cos(u/2)) * np.sin(u) - 
                            (v/2) * np.sin(u/2) * np.cos(u))
        dy_du = lambda u, v: ((self.R + v * np.cos(u/2)) * np.cos(u) - 
                            (v/2) * np.sin(u/2) * np.sin(u))
        dz_du = lambda u, v: (v/2) * np.cos(u/2)
        
        dx_dv = lambda u, v: np.cos(u/2) * np.cos(u)
        dy_dv = lambda u, v: np.cos(u/2) * np.sin(u)
        dz_dv = lambda u, v: np.sin(u/2)
        
        # Cross product magnitude (norm of ∂r/∂u × ∂r/∂v)
        integrand = np.zeros_like(self.U)
        for i in range(self.n):
            for j in range(self.n):
                u, v = self.U[i, j], self.V[i, j]
                # Cross product components
                cp_x = dy_du(u, v) * dz_dv(u, v) - dz_du(u, v) * dy_dv(u, v)
                cp_y = dz_du(u, v) * dx_dv(u, v) - dx_du(u, v) * dz_dv(u, v)
                cp_z = dx_du(u, v) * dy_dv(u, v) - dy_du(u, v) * dx_dv(u, v)
                integrand[i, j] = np.sqrt(cp_x**2 + cp_y**2 + cp_z**2)
        
        # Numerical integration using trapezoidal rule
        surface_area = np.trapz(np.trapz(integrand, dx=dv, axis=0), dx=du, axis=0)
        return surface_area
    
    def compute_edge_length(self):
        """
        Compute the length of the Möbius strip's boundary (v = ±w/2).
        
        Returns:
            float: Total edge length
        """
        # Parametric equations at boundaries (v = ±w/2)
        v_bound = [self.w/2, -self.w/2]
        edge_length = 0.0
        
        for v in v_bound:
            # Parametric equations for the edge
            x = (self.R + v * np.cos(self.u/2)) * np.cos(self.u)
            y = (self.R + v * np.cos(self.u/2)) * np.sin(self.u)
            z = v * np.sin(self.u/2)
            
            # Derivatives w.r.t u
            dx_du = np.gradient(x, self.u)
            dy_du = np.gradient(y, self.u)
            dz_du = np.gradient(z, self.u)
            
            # Integrand for arc length
            integrand = np.sqrt(dx_du**2 + dy_du**2 + dz_du**2)
            
            # Numerical integration
            edge_length += np.trapz(integrand, dx=self.u[1] - self.u[0])
        
        return edge_length
    
    def plot(self):
        """
        Visualize the Möbius strip in 3D and save the plot.
        """
        x, y, z = self.parametric_equations()

        import matplotlib
        matplotlib.use('TkAgg')
        
        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111, projection='3d')
        ax.plot_surface(x, y, z, cmap='viridis', alpha=0.8)
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.set_title('Möbius Strip')
        plt.savefig('mobius_strip.png')
        plt.show()

if __name__ == "__main__":
    # Create Möbius strip with default parameters
    mobius = MobiusStrip(R=1.0, w=0.5, n=100)
    
    # Compute properties
    surface_area = mobius.compute_surface_area()
    edge_length = mobius.compute_edge_length()
    
    # Print results
    print(f"Surface Area: {surface_area:.4f}")
    print(f"Edge Length: {edge_length:.4f}")
    
    # Generate plot
    mobius.plot()