package main

func Pressure(rho, e, u, gamma float64) float64 {
	P := rho * (gamma - 1.0) * (e + .5*u*u)
	return P
}
