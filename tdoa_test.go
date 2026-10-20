package utdoa

import (
	"fmt"
	"testing"
)

func TestTdoaMeasurement_Calculate(t *testing.T) {
	measSet := []TrpMeas{
		{rtoa: rtoa3gppToSec(1000, 0), factor: 0, rsrp: 50},
		{rtoa: rtoa3gppToSec(2000, 1), factor: 1, rsrp: 90},
		{rtoa: rtoa3gppToSec(3000, 2), factor: 2, rsrp: 70},
	}

	var tm TdoaMeasurement
	tm.Calculate(measSet)

	if tm.TrpRefIdx != 1 {
		t.Errorf("TrpRefIdx = %d, want 1 (highest RSRP)", tm.TrpRefIdx)
	}

	assertSlice(t, "Tdoas", tm.Tdoas, []float64{
		LIGHTSPEED * (measSet[1].rtoa - measSet[0].rtoa),
		LIGHTSPEED * (measSet[1].rtoa - measSet[2].rtoa),
	}, 1e-3)
}

func TestTdoaMeasurement_Calculate_ReferenceSelection(t *testing.T) {
	tests := []struct {
		name    string
		rsrp    []uint8
		wantRef int
	}{
		{"strongest in the middle", []uint8{10, 90, 30}, 1},
		{"strongest is last", []uint8{10, 20, 30}, 2},
		{"strongest is first", []uint8{99, 20, 30}, 0},
		{"tie keeps the first", []uint8{80, 80, 10}, 0},
		{"all zero falls back to first", []uint8{0, 0, 0}, 0},
	}
	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			measSet := make([]TrpMeas, len(tc.rsrp))
			for i, rsrp := range tc.rsrp {
				measSet[i] = TrpMeas{rtoa: rtoa3gppToSec(float64(1000*(i+1)), 0), rsrp: rsrp}
			}

			var tm TdoaMeasurement
			tm.Calculate(measSet)

			if tm.TrpRefIdx != tc.wantRef {
				t.Errorf("TrpRefIdx = %d, want %d", tm.TrpRefIdx, tc.wantRef)
			}
			if len(tm.Tdoas) != len(tc.rsrp)-1 {
				t.Errorf("len(Tdoas) = %d, want %d", len(tm.Tdoas), len(tc.rsrp)-1)
			}
		})
	}
}

func TestTdoaMeasurement_Calculate_TwoTrps(t *testing.T) {
	measSet := []TrpMeas{
		{rtoa: rtoa3gppToSec(1000, 0), rsrp: 10},
		{rtoa: rtoa3gppToSec(1500, 0), rsrp: 20},
	}
	var tm TdoaMeasurement
	tm.Calculate(measSet)

	if tm.TrpRefIdx != 1 {
		t.Errorf("TrpRefIdx = %d, want 1", tm.TrpRefIdx)
	}
	tdoa := LIGHTSPEED * (measSet[1].rtoa - measSet[0].rtoa)
	assertSlice(t, "Tdoas", tm.Tdoas, []float64{tdoa}, 1e-3)
}

func TestTdoaMeasurement_CovarianceMatrix(t *testing.T) {
	measSet := []TrpMeas{
		{tdoaQuality: MeasQuality{quality: 1, res: 1}},
		{tdoaQuality: MeasQuality{quality: 1, res: 2}},
		{tdoaQuality: MeasQuality{quality: 3, res: 1}},
	}
	tm := TdoaMeasurement{TrpRefIdx: 1}

	cov := tm.CovarianceMatrix(measSet)
	if cov == nil {
		t.Fatal("covariance is nil")
	}
	assertMatrix(t, "cov", cov, [][]float64{
		{4 + 1, 4},
		{4, 4 + 9},
	}, 1e-12)
}

func TestTdoaMeasurement_CovarianceMatrix_TwoTrps(t *testing.T) {
	measSet := []TrpMeas{
		{tdoaQuality: MeasQuality{quality: 2, res: 1}},
		{tdoaQuality: MeasQuality{quality: 1, res: 3}},
	}
	cov := (TdoaMeasurement{TrpRefIdx: 0}).CovarianceMatrix(measSet)
	if cov == nil {
		t.Fatal("covariance is nil")
	}
	assertMatrix(t, "cov", cov, [][]float64{{13}}, 1e-12)
}

func TestTdoaMeasurement_CovarianceMatrix_ZeroVarianceReturnsNil(t *testing.T) {
	good := MeasQuality{quality: 1, res: 1}
	tests := []struct {
		name    string
		measSet []TrpMeas
		ref     int
	}{
		{"reference has zero variance", []TrpMeas{{tdoaQuality: MeasQuality{}}, {tdoaQuality: good}, {tdoaQuality: good}}, 0},
		{"non-reference has zero variance", []TrpMeas{{tdoaQuality: good}, {tdoaQuality: MeasQuality{}}, {tdoaQuality: good}}, 0},
	}
	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			if cov := (TdoaMeasurement{TrpRefIdx: tc.ref}).CovarianceMatrix(tc.measSet); cov != nil {
				t.Errorf("expected nil covariance, got %v", cov)
			}
		})
	}
}

func TestTdoaMeasurement_JacobianResidual_HandComputed(t *testing.T) {
	trps := newTrps([]float64{10, 0}, []float64{0, 10}, []float64{-10, 0})
	ue := newVec(0, 0)
	tm := TdoaMeasurement{TrpRefIdx: 0, Dim: 2, Tdoas: []float64{1.5, -0.5}}
	ctx := allTrp(3)

	j, r := tm.JacobianResidual(ue, trps, []float64{10, 10, 10}, ctx)

	assertSliceSlice(t, "J", j, [][]float64{
		{-1, 1},
		{-2, 0},
	}, 1e-12)
	assertSlice(t, "r", r, []float64{-1.5, 0.5}, 1e-12)
}

func TestTdoaMeasurement_JacobianResidual_ZeroResidualForConsistentTdoa(t *testing.T) {
	trps := newTrps([]float64{0, 0}, []float64{100, 0}, []float64{0, 100})
	ue := newVec(30, 40)
	d := distancesTo(ue, trps)
	tm := TdoaMeasurement{
		TrpRefIdx: 1,
		Dim:       2,
		Tdoas:     []float64{d[1] - d[0], d[1] - d[2]},
	}

	_, r := tm.JacobianResidual(ue, trps, d, allTrp(3))

	assertSlice(t, "r", r, []float64{0, 0}, 1e-9)
}

func TestTdoaMeasurement_JacobianResidual_ZeroReferenceDistanceReturnsNil(t *testing.T) {
	trps := newTrps([]float64{0, 0}, []float64{10, 0}, []float64{0, 10})
	tm := TdoaMeasurement{TrpRefIdx: 0, Dim: 2, Tdoas: []float64{0, 0}}

	j, r := tm.JacobianResidual(newVec(0, 0), trps, []float64{0, 10, 10}, allTrp(3))

	if j != nil || r != nil {
		t.Errorf("expected (nil, nil), got (%v, %v)", j, r)
	}
}

func TestTdoaMeasurement_JacobianMatchesNumericalDerivative(t *testing.T) {
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
			set := newTrps(positions...)
			ue := newVec([]float64{22, 31, 1.5}[:dim]...)
			tm := &TdoaMeasurement{TrpRefIdx: 2, Dim: dim, Tdoas: []float64{3, -7, 12}}

			analytic, _ := tm.JacobianResidual(ue, set, distancesTo(ue, set), allTrp(len(set)))
			numeric := numericJacobian(t, tm, ue, set)

			assertSliceNearMatrix(t, "J", analytic, numeric, 1e-5)
		})
	}
}

func TestTdoaMeasurement_JacobianResidual_DroppedTrpKeepsTdoaAlignment(t *testing.T) {
	trps := newTrps(
		[]float64{10, 0},
		[]float64{0, 10},
		[]float64{-10, 0},
		[]float64{0, -10},
	)
	ue := newVec(0, 0)
	tm := TdoaMeasurement{TrpRefIdx: 0, Dim: 2, Tdoas: []float64{1, 2, 3}}
	ctx := &PositioningContext{KeptTrp: []int{0, 2, 3}}

	j, r := tm.JacobianResidual(ue, trps, []float64{10, 10, 10, 10}, ctx)

	assertSliceSlice(t, "J", j, [][]float64{
		{-2, 0},
		{-1, -1},
	}, 1e-12)
	assertSlice(t, "r", r, []float64{-2, -3}, 1e-12)
}
