package utdoa

import (
	"math"
	"slices"
	"testing"
)

func TestZenithMeasurement_SetDim(t *testing.T) {
	var zm ZenithMeasurement
	zm.SetDim(3)
	if zm.Dim != 3 {
		t.Errorf("Dim = %d, want 3", zm.Dim)
	}
}

func TestZenithMeasurement_Calculate(t *testing.T) {
	tests := []struct {
		name string
		raw  []float64
		want []float64
	}{
		{"raw 0 is 0 deg", []float64{0}, []float64{0}},
		{"raw 450 is 45 deg", []float64{450}, []float64{math.Pi / 4}},
		{"raw 900 is 90 deg", []float64{900}, []float64{math.Pi / 2}},
		{"raw 1800 is 180 deg", []float64{1800}, []float64{math.Pi}},
		{"order preserved", []float64{1800, 0, 900}, []float64{math.Pi, 0, math.Pi / 2}},
		{"no TRPs", nil, []float64{}},
	}
	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			measSet := make([]TrpMeas, len(tc.raw))
			for i, raw := range tc.raw {
				measSet[i].zenith = raw
			}

			var zm ZenithMeasurement
			zm.Calculate(measSet)

			assertSlice(t, "Zeniths", zm.Zeniths, tc.want, 1e-9)
		})
	}
}

func TestZenithMeasurement_CovarianceMatrix(t *testing.T) {
	deg := math.Pi / 180
	measSet := []TrpMeas{
		{zenithQuality: MeasQuality{quality: 2, res: 0.5}},
		{zenithQuality: MeasQuality{quality: 1, res: 3}},
	}

	cov := ZenithMeasurement{}.CovarianceMatrix(measSet)
	if cov == nil {
		t.Fatal("covariance is nil")
	}

	assertMatrix(t, "cov", cov, [][]float64{
		{deg * deg, 0},
		{0, 9 * deg * deg},
	}, 1e-15)
}

func TestZenithMeasurement_CovarianceMatrix_ZeroVarianceReturnsNil(t *testing.T) {
	measSet := []TrpMeas{
		{zenithQuality: MeasQuality{quality: 1, res: 1}},
		{zenithQuality: MeasQuality{quality: 0, res: 1}},
	}
	if cov := (ZenithMeasurement{}).CovarianceMatrix(measSet); cov != nil {
		t.Errorf("expected nil covariance, got %v", cov)
	}
}

func TestZenithMeasurement_CovarianceUsesZenithQualityOnly(t *testing.T) {
	measSet := []TrpMeas{{
		azimuthQuality: MeasQuality{quality: 5, res: 5},
		tdoaQuality:    MeasQuality{quality: 5, res: 5},
	}}
	if cov := (ZenithMeasurement{}).CovarianceMatrix(measSet); cov != nil {
		t.Errorf("expected nil covariance when zenith quality is unset, got %v", cov)
	}
}

func TestZenithMeasurement_JacobianResidual_KnownGeometry(t *testing.T) {
	trps := newTrps([]float64{0, 0, 0})
	ue := newVec(3, 4, 12)
	dist := []float64{13}
	expectedZenith := math.Atan2(5, 12)

	t.Run("zero residual when measurement matches", func(t *testing.T) {
		zm := ZenithMeasurement{Dim: 3, Zeniths: []float64{expectedZenith}}
		ctx := allTrp(1)

		j, r := zm.JacobianResidual(ue, trps, dist, ctx)

		assertSliceSlice(t, "J", j, [][]float64{{
			12.0 * 3.0 / 845.0,
			12.0 * 4.0 / 845.0,
			-5.0 / 169.0,
		}}, 1e-12)
		assertSlice(t, "r", r, []float64{0}, 1e-12)
		if !slices.Equal(ctx.KeptTrp, []int{0}) {
			t.Errorf("KeptTrp = %v, want [0]", ctx.KeptTrp)
		}
	})

	t.Run("residual is expected minus measured", func(t *testing.T) {
		zm := ZenithMeasurement{Dim: 3, Zeniths: []float64{expectedZenith - 0.25}}
		_, r := zm.JacobianResidual(ue, trps, dist, allTrp(1))
		assertSlice(t, "r", r, []float64{0.25}, 1e-12)
	})
}

func TestZenithMeasurement_JacobianResidual_SkipsTrpDirectlyAboveUe(t *testing.T) {
	trps := newTrps([]float64{3, 4, 100}, []float64{0, 0, 0})
	ue := newVec(3, 4, 12)
	zm := ZenithMeasurement{Dim: 3, Zeniths: []float64{0, math.Atan2(5, 12)}}
	ctx := allTrp(2)

	j, r := zm.JacobianResidual(ue, trps, distancesTo(ue, trps), ctx)

	if !slices.Equal(ctx.KeptTrp, []int{1}) {
		t.Fatalf("KeptTrp = %v, want [1]", ctx.KeptTrp)
	}
	rows, cols := len(j), len(j[0])
	if rows != 1 || cols != 3 {
		t.Fatalf("J is %dx%d, want 1x3", rows, cols)
	}
	assertSlice(t, "r", r, []float64{0}, 1e-12)
}

func TestZenithMeasurement_JacobianResidual_HonoursKeptTrp(t *testing.T) {
	trps := newTrps([]float64{0, 0, 0}, []float64{50, 50, 50}, []float64{-3, -4, 0})
	ue := newVec(3, 4, 12)
	zm := ZenithMeasurement{Dim: 3, Zeniths: []float64{math.Atan2(5, 12), 0, 0.1}}
	ctx := &PositioningContext{KeptTrp: []int{0, 2}}

	j, r := zm.JacobianResidual(ue, trps, distancesTo(ue, trps), ctx)

	rows := len(j)
	if rows != 2 || len(r) != 2 {
		t.Fatalf("got %d Jacobian rows and %d residuals, want 2 and 2", rows, len(r))
	}
	assertNear(t, "r[1]", r[1], math.Atan2(10, 12)-0.1, 1e-12)
}

func TestZenithMeasurement_JacobianMatchesNumericalDerivative(t *testing.T) {
	trps := newTrps(
		[]float64{0, 0, 10},
		[]float64{100, 5, 20},
		[]float64{-30, 80, 5},
		[]float64{60, -70, 15},
	)
	ue := newVec(22, 31, 1.5)
	zm := &ZenithMeasurement{Dim: 3, Zeniths: []float64{1.0, 2.0, 0.5, 1.5}}

	analytic, _ := zm.JacobianResidual(ue, trps, distancesTo(ue, trps), allTrp(len(trps)))
	numeric := numericJacobian(t, zm, ue, trps)

	assertSliceNearMatrix(t, "J", analytic, numeric, 1e-6)
}
