package main

import (
	"encoding/csv"
	"fmt"
	"math"
	"os"
)

// #### Problem constants ####//
// Environmental parameters
const g float64 = 9.81    // [m / s^2]] gravitational acceleration
const rho float64 = 1.091 // [kg / m^3] air density

// Rocket parameters
const r float64 = 0.5             // [m] rocket radius
const A float64 = math.Pi * r * r // [m^2]] rocket cross-section
const ve float64 = 325            // [m/s] propellant exit velocity
const Cd float64 = 0.15           // Drag coefficient
const ms float64 = 50             // [kg] Weight of rocket shell
const mp0 float64 = 100           // [kg] Initial weight of fuel

type state struct {
	h map[float64]float64
	v map[float64]float64
}

// Propellant time parameters
func mp_dot(t float64) float64 { // [kg/s] propellant burn rate as a function of time
	if (t >= 0.0) && (t <= 5.0) {
		return 20.0
	} else {
		return 0
	}
}

func mp(t float64) float64 { // This can be found by integrating mp_dot, and doing so is trivial
	if t < 0 {
		return mp0
	} else if t <= 5 {
		return mp0 - (20.0 * t)
	} else {
		return 0
	}
}

// Simulation parameters
const dt float64 = 0.1 // [s]

// drag calculation
func drag(v float64) float64 {
	D := 0.5 * rho * v * math.Abs(v) * Cd * A
	return D
}

func FindForces(t, v float64) float64 {
	D := drag(v)          //drag
	T := mp_dot(t) * ve   //thrust
	G := g * (ms + mp(t)) //gravity

	return (T - D - G)
}

func main() {
	var t float64 = 0
	Rocket := new(state)

	// Declare Simulation Outputs
	// Set initial altitude and velocity to 0

	Rocket.h = make(map[float64]float64)
	Rocket.v = make(map[float64]float64)
	Rocket.h[0] = 0
	Rocket.v[0] = 0

	// h := make(map[float64]float64)
	// v := make(map[float64]float64)

	var a float64
	i := 0
	// Euler step until ground impact
	for Rocket.h[t] >= 0 {
		a = FindForces(t, Rocket.v[t]) / (ms + mp(t))
		Rocket.v[t+dt] = Rocket.v[t] + a*dt
		Rocket.h[t+dt] = Rocket.h[t] + Rocket.v[t]

		i += 1
		t += dt
	}

	f, _ := os.Create("Trajectory.csv")
	defer f.Close()
	w := csv.NewWriter(f)
	w.Write([]string{"time", "altitude", "velocity"})
	for time := range Rocket.h {
		r := make([]string, 0, 3)
		r = append(r, fmt.Sprintf("%f", time))
		r = append(r, fmt.Sprintf("%f", Rocket.h[time]))
		r = append(r, fmt.Sprintf("%f", Rocket.v[time]))
		w.Write(r)
	}
	w.Flush()
}
