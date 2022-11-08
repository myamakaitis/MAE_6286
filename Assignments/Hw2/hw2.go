package main

import (
	"encoding/csv"
	"fmt"
	"os"
	"time"
)

const VmaxA float64 = 90  // [km / hr]
const VmaxB float64 = 130 // [km / hr]
const L float64 = 25      // [km]
const Rmax float64 = 100  // [cars / km]

const dt float64 = 1e-3 // [hr]

const Nx int32 = 101 // Discrete Points

func LinSpace(start, stop float64, n int32) []float64 {
	space := make([]float64, n)
	increment := (stop - start) / (float64(n) - 1)

	for ii := 0; int32(ii) < n; ii++ {
		space[ii] = start + float64(ii)*increment
	}
	space[n-1] = stop

	return space
}

func timeTrack(start time.Time, name string) {
	elapsed := time.Since(start).Nanoseconds() / 1e6
	fmt.Printf("%s took %d ms", name, elapsed)

}

type Highway struct {
	L    float64 // Length
	N    int32   // Number of points
	Vmax float64 // max velocity
	Rmax float64 // max density

	dL float64 // Point Spacing

	X  []float64 // coordinates
	R  []float64 // density
	Rc []float64 // copy of density
	V  []float64 // Velocity
	F  []float64 // Flux

	dR []float64 // dr wrt t
	dF []float64 // dF wrt x

}

func (way *Highway) init(rho0 func(float64) float64) { // road structure to hold all the properties of a road
	way.X = LinSpace(0, way.L, way.N)
	way.dL = way.X[1] - way.X[0]

	way.R = make([]float64, way.N)
	way.Rc = make([]float64, way.N) // copy of R for runge-kutta 2nd order

	way.V = make([]float64, way.N)
	way.F = make([]float64, way.N)
	way.dF = make([]float64, way.N)
	way.dR = make([]float64, way.N)

	for i, x := range way.X {
		way.R[i] = rho0(x)
	}

	way.VelFlux()

}

func (way *Highway) VelFlux() { //calculate velocity and flux
	for i, r := range way.R {
		way.V[i] = way.Vmax * (1 - r/way.Rmax)
		way.F[i] = r * way.V[i]
	}
}

func (way *Highway) D() { // calculate dF/dx and dR/dt
	for i := int32(1); i < (way.N - 1); i++ {
		way.dF[i] = (way.F[i] - way.F[i-1]) / way.dL
		way.dR[i] = -way.dF[i]
	}
}

func (way *Highway) dRdt_rk2(dt float64) {
	copy(way.Rc, way.R)
	way.D()                 // calculate dF/dx, dR/dt, using current R
	way.updateVFR(0.5 * dt) // update R using dR/dt and half time step, recalculates V and F
	way.D()                 // re-calculate dF/dx, dR/dt at the half time-step
	copy(way.R, way.Rc)     // re-set the density to be at the beginning of the step
	way.updateVFR(dt)       // update the density using the derivatives at the half step, recalculate V and F for next caclulation
}

func (way *Highway) updateVFR(dt float64) {
	for i := int32(1); i < (way.N - 1); i++ {
		way.R[i] += way.dR[i] * dt
	}
	way.VelFlux()
}

func (way *Highway) write_start(t float64, Vfile *csv.Writer, Rfile *csv.Writer) {
	str_start := make([]string, 0, (1 + way.N))
	str_start = append(str_start, "t")
	for i := int32(0); i < (way.N); i++ {
		str_start = append(str_start, fmt.Sprintf("%f", way.X[i]))
	}
	Vfile.Write(str_start)
	Rfile.Write(str_start)

	way.writeVR(t, Vfile, Rfile)
}

func (way *Highway) writeVR(t float64, Vfile *csv.Writer, Rfile *csv.Writer) {
	strV := make([]string, 0, (1 + way.N))
	strR := make([]string, 0, (1 + way.N))

	strV = append(strV, fmt.Sprintf("%f", t))
	strR = append(strR, fmt.Sprintf("%f", t))

	for i := int32(0); i < (way.N); i++ {
		strV = append(strV, fmt.Sprintf("%f", way.V[i]))
		strR = append(strR, fmt.Sprintf("%f", way.R[i]))
	}

	Vfile.Write(strV)
	Rfile.Write(strR)
}

func R0A(x float64) (r float64) {
	if x >= 2 && x <= 4.2 {
		r = 50
	} else {
		r = 10
	}
	return
}

func R0B(x float64) (r float64) {
	if x >= 2 && x <= 4.2 {
		r = 50
	} else {
		r = 20
	}
	return
}

func main() {
	road := new(Highway)
	road.L = L
	road.N = Nx
	road.Vmax = VmaxA
	road.Rmax = Rmax

	road.init(R0A)

	VhistA, _ := os.Create("VelocityA.csv")
	RhistA, _ := os.Create("DensityA.csv")
	defer VhistA.Close()
	defer RhistA.Close()
	Vfile := csv.NewWriter(VhistA)
	Rfile := csv.NewWriter(RhistA)

	var t float64 = 0
	road.write_start(t, Vfile, Rfile)

	for t = 0; t <= .55; t += dt {
		road.dRdt_rk2(dt)
		road.writeVR(t, Vfile, Rfile)

	}

	road.Vmax = VmaxB
	road.init(R0B)

	VhistB, _ := os.Create("VelocityB.csv")
	RhistB, _ := os.Create("DensityB.csv")
	defer VhistB.Close()
	defer RhistB.Close()
	Vfile = csv.NewWriter(VhistB)
	Rfile = csv.NewWriter(RhistB)

	t = 0
	road.write_start(t, Vfile, Rfile)

	defer timeTrack(time.Now(), "Run B time")

	for t = 0; t <= .55; t += dt {
		road.dRdt_rk2(dt)
		road.writeVR(t, Vfile, Rfile)

	}

}
