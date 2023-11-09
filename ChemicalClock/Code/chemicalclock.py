"""
Implementation of the Chemical Clock model via
cubic autocatalysis.
"""
import numpy as np
from scipy.integrate import solve_ivp

class ChemicalClock:
    def __init__(self,
                 init_conc: np.ndarray,
                 t: np.ndarray,
                 rates: np.ndarray,
                 p_fast_factor: float = 100,
                 q_fast_factor: float = 100,
                 ) -> None:
        """
        A simple attempt at modeling the chemical clock
        based on iodine and iodide concentrations.
        """
        self.init_conc = init_conc
        self.t = t
        self.rates = rates
        self.p_slow = rates[0]
        self.q_slow = rates[1]
        self.p_fast = rates[0] * p_fast_factor
        self.q_fast = rates[1] * q_fast_factor
        pass

    def model_eq(self, t: float, state: np.ndarray) -> np.ndarray:
        """
        The model equations for the chemical clock.
        """
        # Unpack the state vector
        m, n, u, v, w, x, y, z = state

        # Define the rate constants
        p_slow = self.p_slow
        p_fast = self.p_fast
        q_slow = self.q_slow
        q_fast = self.q_fast

        # Define the model equations
        dmdt = -p_slow*m*u - p_fast*m*v + q_slow*x*n + q_fast*x*y
        dndt = p_fast*m*v - q_slow*x*n
        dudt = -p_slow*m*u
        dvdt = -p_fast*m*v + p_slow*m*u
        dwdt = p_fast*m*v
        dxdt = -q_slow*x*n - q_fast*x*y
        dydt = q_slow*x*n - q_fast*x*y
        dzdt = q_fast*x*y

        # Return the model equations
        return np.array([dmdt, dndt, dudt, dvdt, dwdt, dxdt, dydt, dzdt])

    def solve(self):
        """
        Solve the model equations.
        """
        sol = solve_ivp(self.model_eq, (self.t[0], self.t[-1]), 
                        self.init_conc, method='LSODA', t_eval=self.t)
        return sol
    

class SimpleClock(ChemicalClock):
    def __init__(self, init_conc: np.ndarray,
                  t: np.ndarray, rates: np.ndarray, 
                  p_fast_factor: float = 100, q_fast_factor: float = 100) -> None:
        # super().__init__(init_conc, t, rates, p_fast_factor, q_fast_factor)
        self.init_conc = init_conc
        self.t = t
        self.rate = rates[0]


    def model_eq(self, t: float, state: np.ndarray) -> np.ndarray:
        """
        The model equations for the chemical clock.
        """
        # Unpack the state vector
        gamma, beta = state

        # Define the rate constants
        rate_factor = self.rate

        # Define the model equations
        dgamma = - rate_factor * gamma * beta**2
        dbeta = rate_factor * gamma * beta**2

        # Return the model equations
        return np.array([dgamma, dbeta])
    


