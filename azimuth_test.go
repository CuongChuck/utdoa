package utdoa

import (
	"fmt"
	"math"
	"slices"
	"testing"
)

func TestAzimuthMeasurement_SetDim(t *testing.T) {
	var am AzimuthMeasurement
	am.SetDim(3)
	if am.Dim != 3 {
		t.Errorf("Dim = %d, want 3", am.Dim)
	}
	am.SetDim(2)
	if am.Dim != 2 {
		t.Errorf("Dim = %d, want 2", am.Dim)
	}
}

func TestAzimuthMeasurement_Calculate(t *testing.T) {
	tests := []struct {
		name string
		raw  []float64
		want []float64
	}{
		{"raw 0 is -180 deg", []float64{0}, []float64{-math.Pi}},
		{"raw 900 is -90 deg", []float64{900}, []float64{-math.Pi / 2}},
		{"raw 1800 is 0 deg", []float64{1800}, []float64{0}},
		{"raw 2700 is +90 deg", []float64{2700}, []float64{math.Pi / 2}},
		{"raw 3600 is +180 deg", []float64{3600}, []float64{math.Pi}},
		{"one value per TRP, order preserved", []float64{0, 1800, 2700}, []float64{-math.Pi, 0, math.Pi / 2}},
		{"no TRPs", nil, []float64{}},
	}
	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			measSet := make([]TrpMeas, len(tc.raw))
			for i, raw := range tc.raw {
				measSet[i].azimuth = raw
			}

			var am AzimuthMeasurement
			am.Calculate(measSet)

			assertSlice(t, "Azimuths", am.Azimuths, tc.want, 1e-9)
		})
	}
}

func TestAzimuthMeasurement_Calculate_ReplacesPreviousResult(t *testing.T) {
	am := AzimuthMeasurement{Azimuths: []float64{1, 2, 3, 4}}
	am.Calculate([]TrpMeas{{azimuth: 1800}})
	if len(am.Azimuths) != 1 {
		t.Errorf("len(Azimuths) = %d, want 1", len(am.Azimuths))
	}
}

func TestAzimuthMeasurement_CovarianceMatrix(t *testing.T) {
	deg := math.Pi / 180
	measSet := []TrpMeas{
		{azimuthQuality: MeasQuality{quality: 2, res: 0.5}},
		{azimuthQuality: MeasQuality{quality: 1, res: 3}},
	}

	cov := AzimuthMeasurement{}.CovarianceMatrix(measSet)
	if cov == nil {
		t.Fatal("covariance is nil")
	}

	assertMatrix(t, "cov", cov, [][]float64{
		{deg * deg, 0},
		{0, 9 * deg * deg},
	}, 1e-15)
}

func TestAzimuthMeasurement_CovarianceMatrix_ZeroVarianceReturnsNil(t *testing.T) {
	tests := []struct {
		name    string
		measSet []TrpMeas
	}{
		{"zero quality", []TrpMeas{
			{azimuthQuality: MeasQuality{quality: 1, res: 1}},
			{azimuthQuality: MeasQuality{quality: 0, res: 1}},
		}},
		{"zero resolution", []TrpMeas{
			{azimuthQuality: MeasQuality{quality: 1, res: 0}},
			{azimuthQuality: MeasQuality{quality: 1, res: 1}},
		}},
		{"unset quality", []TrpMeas{{}}},
	}
	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			if cov := (AzimuthMeasurement{}).CovarianceMatrix(tc.measSet); cov != nil {
				t.Errorf("expected nil covariance, got %v", cov)
			}
		})
	}
}

func TestAzimuthMeasurement_JacobianResidual_2D(t *testing.T) {
	trps := newTrps([]float64{0, 0})
	ue := newVec(1, 1)

	t.Run("zero residual when measurement matches", func(t *testing.T) {
		am := AzimuthMeasurement{Dim: 2, Azimuths: []float64{math.Pi / 4}}
		ctx := allTrp(1)

		j, r := am.JacobianResidual(ue, trps, nil, ctx)

		assertSliceSlice(t, "J", j, [][]float64{{-0.5, 0.5}}, 1e-12)
		assertSlice(t, "r", r, []float64{0}, 1e-12)
		if !slices.Equal(ctx.KeptTrp, []int{0}) {
			t.Errorf("KeptTrp = %v, want [0]", ctx.KeptTrp)
		}
	})

	t.Run("residual is expected minus measured", func(t *testing.T) {
		am := AzimuthMeasurement{Dim: 2, Azimuths: []float64{0}}
		_, r := am.JacobianResidual(ue, trps, nil, allTrp(1))
		assertSlice(t, "r", r, []float64{math.Pi / 4}, 1e-12)
	})

	t.Run("residual wraps across the +/-pi boundary", func(t *testing.T) {
		ue := newVec(math.Cos(math.Pi-0.1), math.Sin(math.Pi-0.1))
		am := AzimuthMeasurement{Dim: 2, Azimuths: []float64{-math.Pi + 0.1}}
		_, r := am.JacobianResidual(ue, trps, nil, allTrp(1))
		assertSlice(t, "r", r, []float64{-0.2}, 1e-12)
	})
}

func TestAzimuthMeasurement_JacobianResidual_3DAddsZeroColumn(t *testing.T) {
	trps := newTrps([]float64{0, 0, 0})
	ue := newVec(1, 1, 5)
	am := AzimuthMeasurement{Dim: 3, Azimuths: []float64{math.Pi / 4}}

	j, r := am.JacobianResidual(ue, trps, nil, allTrp(1))

	assertSliceSlice(t, "J", j, [][]float64{{-0.5, 0.5, 0}}, 1e-12)
	assertSlice(t, "r", r, []float64{0}, 1e-12)
}

func TestAzimuthMeasurement_JacobianResidual_SkipsTrpAtUePosition(t *testing.T) {
	trps := newTrps([]float64{1, 1}, []float64{0, 0})
	ue := newVec(1, 1)
	am := AzimuthMeasurement{Dim: 2, Azimuths: []float64{0, math.Pi / 4}}
	ctx := allTrp(2)

	j, r := am.JacobianResidual(ue, trps, nil, ctx)

	if !slices.Equal(ctx.KeptTrp, []int{1}) {
		t.Fatalf("KeptTrp = %v, want [1]", ctx.KeptTrp)
	}
	assertSliceSlice(t, "J", j, [][]float64{{-0.5, 0.5}}, 1e-12)
	assertSlice(t, "r", r, []float64{0}, 1e-12)
}

func TestAzimuthMeasurement_JacobianResidual_HonoursKeptTrp(t *testing.T) {
	trps := newTrps([]float64{0, 0}, []float64{10, 10}, []float64{0, 1})
	ue := newVec(1, 1)
	am := AzimuthMeasurement{Dim: 2, Azimuths: []float64{math.Pi / 4, 0, 0.5}}
	ctx := &PositioningContext{KeptTrp: []int{0, 2}}

	j, r := am.JacobianResidual(ue, trps, nil, ctx)

	assertSliceSlice(t, "J", j, [][]float64{{-0.5, 0.5}, {0, 1}}, 1e-12)
	assertSlice(t, "r", r, []float64{0, -0.5}, 1e-12)
	if !slices.Equal(ctx.KeptTrp, []int{0, 2}) {
		t.Errorf("KeptTrp = %v, want [0 2]", ctx.KeptTrp)
	}
}

func TestAzimuthMeasurement_JacobianMatchesNumericalDerivative(t *testing.T) {
	full := [][]float64{
		{0, 0, 10},
		{100, 5, 20},
		{-30, 80, 5},
		{60, -70, 15},
	}
	for _, dim := range []int{2, 3} {
		t.Run(fmt.Sprintf("%dD", dim), func(t *testing.T) {
			positions := make([][]float64, len(full))
			for i, p := range full {
				positions[i] = p[:dim]
			}
			trps := newTrps(positions...)
			ue := newVec([]float64{22, 31, 1.5}[:dim]...)

			am := &AzimuthMeasurement{Dim: dim, Azimuths: []float64{0.3, -1.2, 2.0, 0.9}}

			analytic, _ := am.JacobianResidual(ue, trps, nil, allTrp(len(trps)))
			numeric := numericJacobian(t, am, ue, trps)

			assertSliceNearMatrix(t, "J", analytic, numeric, 1e-6)
		})
	}
}
