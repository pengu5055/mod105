"""
Contains the base class HydrogenFusion, which is responsible for 
running the simulation of the second task.
"""
import numpy as np
from scipy.integrate import solve_ivp
from typing import Iterable

class HydrogenFusion:
    def __init__(self,
                 init_conc: Iterable[float],
                 t_eval: Iterable[float],
                 rates: Iterable[float],
                 ) -> None:
        """
        Initializes the HydrogenFusion class.
        """
        self.init_conc = init_conc
        self.t_eval = t_eval
        self.rates = rates

    def rate_eq(self, t_eval, state):
        """
        Wrapper that combines the rate equations into one function.
        """
        return np.array([self._du(t_eval, state),
                         self._dv(t_eval, state),
                         self._dx(t_eval, state),
                         self._dy(t_eval, state),
                         self._dz(t_eval, state)])
    
    def solve(self):
        """
        Solve the rate equations using scipy's solve_ivp function.
        """
        sol = solve_ivp(self.rate_eq,
                        t_span=(self.t_eval[0], self.t_eval[-1]),
                        y0=self.init_conc,
                        t_eval=self.t_eval,)
        return sol
        

    def _du(self, t_eval, state):
        # Unpack the state vector
        u, v, x, y, z = state

        # Get rate coefficients
        p, q, r, s, t = self.rates

        return s*x*y - r*u*z
    
    def _dv(self, t_eval, state):
        # Unpack the state vector
        u, v, x, y, z = state

        # Get rate coefficients
        p, q, r, s, t = self.rates

        return q*z**2 - p*v - t*v*y
    
    def _dx(self, t_eval, state):
        # Unpack the state vector
        u, v, x, y, z = state

        # Get rate coefficients
        p, q, r, s, t = self.rates

        return t*v*y - s*x*y + r*u*z
    
    def _dy(self, t_eval, state):
        # Unpack the state vector
        u, v, x, y, z = state

        # Get rate coefficients
        p, q, r, s, t = self.rates

        return r*u*z - t*v*y - s*x*y
    
    def _dz(self, t_eval, state):
        # Unpack the state vector
        u, v, x, y, z = state

        # Get rate coefficients
        p, q, r, s, t = self.rates

        return t*v*y + p*v - q*z**2
    
class HydrogenFusionStationaryH(HydrogenFusion):
    def __init__(self, init_conc: Iterable[float], t_eval: Iterable[float], rates: Iterable[float]) -> None:
        """
        A stationary solution for the rate equations where the concentration of
        hydrogen is constant. Thus the rate equations for hydrogen are set to zero.
        That way we can write the term:

        s*x*y = r*u*z - t*v*y
        """
        super().__init__(init_conc, t_eval, rates)

    def _du(self, t_eval, state):
        # Unpack the state vector
        u, v, x, y, z = state

        # Get rate coefficients
        p, q, r, s, t = self.rates

        return -t*v*y 
    
    def _dv(self, t_eval, state):
        # Unpack the state vector
        u, v, x, y, z = state

        # Get rate coefficients
        p, q, r, s, t = self.rates

        return q*z**2 - p*v - t*v*y
    
    def _dx(self, t_eval, state):
        # Unpack the state vector
        u, v, x, y, z = state

        # Get rate coefficients
        p, q, r, s, t = self.rates

        return 2*t*v*y
    
    def _dy(self, t_eval, state):
        # Unpack the state vector
        u, v, x, y, z = state

        # Get rate coefficients
        p, q, r, s, t = self.rates

        return 0
    
    def _dz(self, t_eval, state):
        # Unpack the state vector
        u, v, x, y, z = state

        # Get rate coefficients
        p, q, r, s, t = self.rates

        return t*v*y + p*v - q*z**2
    
class HydrogenFusionStationaryBr(HydrogenFusion):
    def __init__(self, init_conc: Iterable[float], t_eval: Iterable[float], rates: Iterable[float]) -> None:
        """
        A stationary solution for the rate equations where the concentration of
        bromine is constant. Thus the rate equations for bromine are set to zero.
        That way we can write the term:

        t*v*y = q*z**2 - p*v
        """
        super().__init__(init_conc, t_eval, rates)

    def _du(self, t_eval, state):
        # Unpack the state vector
        u, v, x, y, z = state

        # Get rate coefficients
        p, q, r, s, t = self.rates

        return s*x*y - r*u*z
    
    def _dv(self, t_eval, state):
        # Unpack the state vector
        u, v, x, y, z = state

        # Get rate coefficients
        p, q, r, s, t = self.rates

        return 0
    
    def _dx(self, t_eval, state):
        # Unpack the state vector
        u, v, x, y, z = state

        # Get rate coefficients
        p, q, r, s, t = self.rates
        
        return q*z**2 - p*v - s*x*y + r*u*z
    
    def _dy(self, t_eval, state):
        # Unpack the state vector
        u, v, x, y, z = state

        # Get rate coefficients
        p, q, r, s, t = self.rates

        return r*u*z - s*x*y - q*z**2 + p*v
    
    def _dz(self, t_eval, state):
        # Unpack the state vector
        u, v, x, y, z = state
        
        # Get rate coefficients
        p, q, r, s, t = self.rates

        return 0

