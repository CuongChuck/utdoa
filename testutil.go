package utdoa

import (
	"fmt"
	"math"
	"slices"
	"testing"

	"gonum.org/v1/gonum/mat"
)

var (
	_ MeasurementBlock = (*AzimuthMeasurement)(nil)
	_ MeasurementBlock = (*ZenithMeasurement)(nil)
	_ MeasurementBlock = (*TdoaMeasurement)(nil)
	_ UncertaintyShape = (*UncertaintyEllipse2D)(nil)
	_ UncertaintyShape = (*UncertaintyEllipsoid3D)(nil)
)

const (
	chi2P90Dof2 = 4.605170185988092
	chi2P90Dof3 = 6.251388631170325
)

const commonClockOffsetUnits = 12345.0

func synthMeas(ue []float64, trps [][]float64, rsrp []uint8, tdoaSigma float64) []TrpMeas {
	ueVec := newVec(slices.Clone(ue)...)
	meas := make([]TrpMeas, len(trps))
	for i, p := range trps {
		pos := newVec(slices.Clone(p)...)
		dist := math.Sqrt(squaredDist(ueVec, pos))
		units := dist/LIGHTSPEED/BASE_TIME_UNIT + commonClockOffsetUnits

		dx, dy := ue[0]-p[0], ue[1]-p[1]
		azimuthDeg := math.Atan2(dy, dx) * DEGREE_PER_PI / math.Pi
		zenithDeg := 0.0
		if len(ue) == 3 {
			zenithDeg = math.Atan2(math.Hypot(dx, dy), ue[2]-p[2]) * DEGREE_PER_PI / math.Pi
		}

		meas[i] = TrpMeas{
			id:             i + 1,
			rtoa3gpp:       units - LOWEST_RTOA,
			factor:         0,
			rsrp:           rsrp[i],
			isHighAcc:      false,
			pos3gpp:        pos,
			azimuth:        (azimuthDeg + DEGREE_PER_PI) / ANGLE_RESOLUTION,
			zenith:         zenithDeg / ANGLE_RESOLUTION,
			tdoaQuality:    MeasQuality{quality: tdoaSigma, res: 1},
			azimuthQuality: MeasQuality{quality: 1, res: 1},
			zenithQuality:  MeasQuality{quality: 1, res: 1},
		}
	}
	return meas
}

var (
	square2D = [][]float64{{0, 0}, {100, 0}, {0, 100}, {100, 100}}
	rsrp2D   = []uint8{50, 90, 70, 60}

	trps3D = [][]float64{
		{0, 0, 10},
		{100, 0, 30},
		{0, 100, 20},
		{100, 100, 40},
		{60, 20, 5},
	}
	rsrp3D = []uint8{50, 90, 70, 60, 80}
)

func distSlice(pos1, pos2 []float64) float64 {
	dim := len(pos1)
	dist := 0.0
	for idx := range dim {
		distDim := pos1[idx] - pos2[idx]
		dist += distDim * distDim
	}
	return math.Sqrt(dist)
}

func gcsTo3gpp(pos []float64, isHighFactor bool) *mat.VecDense {
	dim := len(pos)
	result := make([]float64, dim)
	if isHighFactor {
		result[1] = pos[1] / LONGITUDE_CONST * math.Pow(2.0, HIGH_FACTOR)
		result[0] = pos[0] / LATITUDE_CONST * math.Pow(2.0, HIGH_FACTOR)
		if dim == 3 {
			result[2] = pos[2] / math.Pow(2.0, -7.0)
		}
	} else {
		result[1] = pos[1] / LONGITUDE_CONST * math.Pow(2.0, DIRECT_FACTOR)
		result[0] = pos[0] / LATITUDE_CONST * math.Pow(2.0, DIRECT_FACTOR)
		if dim == 3 {
			result[2] = pos[2]
		}
	}
	return mat.NewVecDense(dim, result)
}

func pcsToGcsSlice(pos []float64, origin []float64) []float64 {
	dim := len(pos)
	result := make([]float64, dim)
	lat0Rad := origin[0] * DEG_TO_RAD
	result[0] = origin[0] + (pos[0]/EARTH_RADIUS)/DEG_TO_RAD
	result[1] = origin[1] + (pos[1]/(EARTH_RADIUS*math.Cos(lat0Rad)))/DEG_TO_RAD
	if dim == 3 {
		result[2] = pos[2]
	}
	return result
}

func gcsToRtoa(trpPosPcs, uePosPcs []float64, factor float64) float64 {
	dist := distSlice(trpPosPcs, uePosPcs)
	rtoa := (dist/LIGHTSPEED/BASE_TIME_UNIT - LOWEST_RTOA) / math.Pow(BASE_TIME, factor)
	return rtoa
}

func gcsToPcsSlice(pos []float64, origin []float64) []float64 {
	dim := len(pos)
	result := make([]float64, dim)
	lat0Rad := origin[0] * DEG_TO_RAD
	dLong := (pos[1] - origin[1]) * DEG_TO_RAD
	dLat := (pos[0] - origin[0]) * DEG_TO_RAD
	result[1] = dLong * EARTH_RADIUS * math.Cos(lat0Rad)
	result[0] = dLat * EARTH_RADIUS
	if dim == 3 {
		result[2] = pos[2]
	}
	return result
}

func pcsToAzimuth3gpp(uePos, trpPos []float64) float64 {
	azimuth := math.Atan2(uePos[1]-trpPos[1], uePos[0]-trpPos[0])
	return (azimuth/DEG_TO_RAD + DEGREE_PER_PI) / ANGLE_RESOLUTION
}

func pcsToZenith3gpp(uePos, trpPos []float64) float64 {
	dX := uePos[0] - trpPos[0]
	dY := uePos[1] - trpPos[1]
	dZ := uePos[2] - trpPos[2]
	zenith := math.Atan2(math.Hypot(dX, dY), dZ)
	return zenith / DEG_TO_RAD / ANGLE_RESOLUTION
}

func assertNear(t *testing.T, name string, got, want, tol float64) {
	t.Helper()
	if math.IsNaN(got) || math.Abs(got-want) > tol {
		t.Errorf("%s got %.12g, want %.12g (±%g)", name, got, want, tol)
	}
}

func assertSlice(t *testing.T, name string, got, want []float64, tol float64) {
	t.Helper()
	if len(got) != len(want) {
		t.Fatalf("%s has length %d, want %d", name, len(got), len(want))
	}
	for i := range want {
		assertNear(t, fmt.Sprintf("%s[%d]", name, i), got[i], want[i], tol)
	}
}

func assertSliceSlice(t *testing.T, name string, got, want [][]float64, tol float64) {
	t.Helper()
	if len(got) != len(want) {
		t.Fatalf("%s has length %d, want %d", name, len(got), len(want))
	}
	for i := range want {
		for j := range want[i] {
			assertNear(t, fmt.Sprintf("%s[%d]", name, i), got[i][j], want[i][j], tol)
		}
	}
}

func assertVec(t *testing.T, name string, got mat.Vector, want []float64, tol float64) {
	t.Helper()
	if got.Len() != len(want) {
		t.Fatalf("%s has length %d, want %d", name, got.Len(), len(want))
	}
	for i, w := range want {
		assertNear(t, fmt.Sprintf("%s[%d]", name, i), got.AtVec(i), w, tol)
	}
}

func assertMatrix(t *testing.T, name string, got mat.Matrix, want [][]float64, tol float64) {
	t.Helper()
	r, c := got.Dims()
	wantC := 0
	if len(want) > 0 {
		wantC = len(want[0])
	}
	if r != len(want) || c != wantC {
		t.Fatalf("%s is %dx%d, want %dx%d", name, r, c, len(want), wantC)
	}
	for i := range want {
		for j := range want[i] {
			assertNear(t, fmt.Sprintf("%s[%d][%d]", name, i, j), got.At(i, j), want[i][j], tol)
		}
	}
}

func assertMatricesNear(t *testing.T, name string, got, want mat.Matrix, tol float64) {
	t.Helper()
	r, c := got.Dims()
	wr, wc := want.Dims()
	if r != wr || c != wc {
		t.Fatalf("%s is %dx%d, want %dx%d", name, r, c, wr, wc)
	}
	for i := range r {
		for j := range c {
			assertNear(t, fmt.Sprintf("%s[%d][%d]", name, i, j), got.At(i, j), want.At(i, j), tol)
		}
	}
}

func assertSliceNearMatrix(t *testing.T, name string, got [][]float64, want mat.Matrix, tol float64) {
	t.Helper()
	r, c := len(got), len(got[0])
	wr, wc := want.Dims()
	if r != wr || c != wc {
		t.Fatalf("%s is %dx%d, want %dx%d", name, r, c, wr, wc)
	}
	for i := range r {
		for j := range c {
			assertNear(t, fmt.Sprintf("%s[%d][%d]", name, i, j), got[i][j], want.At(i, j), tol)
		}
	}
}

func newVec(v ...float64) *mat.VecDense {
	return mat.NewVecDense(len(v), v)
}

func identitySlice(n int) [][]float64 {
	m := make([][]float64, n)
	for i := range n {
		m[i] = make([]float64, n)
		m[i][i] = 1
	}
	return m
}

func newTrps(positions ...[]float64) []TrpMeas {
	out := make([]TrpMeas, len(positions))
	for i, p := range positions {
		out[i] = TrpMeas{
			id:  i,
			pos: mat.NewVecDense(len(p), slices.Clone(p)),
		}
	}
	return out
}

func distancesTo(ue *mat.VecDense, measSet []TrpMeas) []float64 {
	d := make([]float64, len(measSet))
	for i := range measSet {
		d[i] = math.Sqrt(squaredDist(ue, measSet[i].pos))
	}
	return d
}

func allTrp(n int) *PositioningContext {
	kept := make([]int, n)
	for i := range kept {
		kept[i] = i
	}
	return &PositioningContext{KeptTrp: kept}
}

func numericJacobian(
	t *testing.T,
	block MeasurementBlock,
	ue *mat.VecDense,
	measSet []TrpMeas,
) *mat.Dense {
	t.Helper()
	const h = 1e-6

	residualAt := func(x *mat.VecDense) []float64 {
		_, r := block.JacobianResidual(
			x, measSet, distancesTo(x, measSet), allTrp(len(measSet)),
		)
		return r
	}

	n := ue.Len()
	var out *mat.Dense
	for k := range n {
		plus := mat.VecDenseCopyOf(ue)
		minus := mat.VecDenseCopyOf(ue)
		plus.SetVec(k, ue.AtVec(k)+h)
		minus.SetVec(k, ue.AtVec(k)-h)
		rp, rm := residualAt(plus), residualAt(minus)
		if out == nil {
			out = mat.NewDense(len(rp), n, nil)
		}
		for i := range len(rp) {
			out.Set(i, k, (rp[i]-rm[i])/(2*h))
		}
	}
	return out
}

type mockBlock struct {
	dim      int
	calls    int
	seenKept [][]int
	jr       func(ue *mat.VecDense, ctx *PositioningContext) ([][]float64, []float64)
}

func (m *mockBlock) Calculate([]TrpMeas)                      {}
func (m *mockBlock) CovarianceMatrix([]TrpMeas) *mat.SymDense { return nil }
func (m *mockBlock) SetDim(dim int)                           { m.dim = dim }

func (m *mockBlock) JacobianResidual(
	ue *mat.VecDense,
	_ []TrpMeas,
	_ []float64,
	ctx *PositioningContext,
) ([][]float64, []float64) {
	m.calls++
	m.seenKept = append(m.seenKept, slices.Clone(ctx.KeptTrp))
	return m.jr(ue, ctx)
}

func (m *mockBlock) FilterJacobianResidual(
	jacobian [][]float64,
	residual []float64,
	ctx *PositioningContext,
) ([][]float64, []float64) {
	return jacobian, residual
}

func fixedBlock(j [][]float64, r []float64) *mockBlock {
	return &mockBlock{
		jr: func(*mat.VecDense, *PositioningContext) ([][]float64, []float64) {
			return j, r
		},
	}
}

func pullToward(target ...float64) *mockBlock {
	n := len(target)
	return &mockBlock{
		jr: func(ue *mat.VecDense, _ *PositioningContext) ([][]float64, []float64) {
			r := make([]float64, n)
			for i := range n {
				r[i] = ue.AtVec(i) - target[i]
			}
			return identitySlice(n), r
		},
	}
}
