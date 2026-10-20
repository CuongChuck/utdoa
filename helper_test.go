package utdoa

import (
	"math"
	"strings"
	"testing"

	"gonum.org/v1/gonum/mat"
)

func TestSquaredDist(t *testing.T) {
	tests := []struct {
		name string
		a, b []float64
		want float64
	}{
		{"2D 3-4-5 triangle", []float64{0, 0}, []float64{3, 4}, 25},
		{"3D", []float64{1, 2, 3}, []float64{4, 6, 3}, 25},
		{"identical points", []float64{7, -2, 5}, []float64{7, -2, 5}, 0},
		{"negative coordinates", []float64{-1, -1}, []float64{1, 1}, 8},
	}
	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			got := squaredDist(newVec(tc.a...), newVec(tc.b...))
			assertNear(t, "squaredDist", got, tc.want, 1e-12)

			rev := squaredDist(newVec(tc.b...), newVec(tc.a...))
			assertNear(t, "squaredDist (reversed)", rev, tc.want, 1e-12)
		})
	}
}

func TestCalculateKChiSquared_KnownValues(t *testing.T) {
	tests := []struct {
		name string
		p, d float64
		want float64
	}{
		{"90% / 2 dof", 0.90, 2, 4.605170186},
		{"90% / 3 dof", 0.90, 3, 6.251388631},
		{"95% / 1 dof", 0.95, 1, 3.841458821},
		{"50% / 2 dof", 0.50, 2, 1.386294361},
	}
	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			assertNear(t, "quantile", calculateKChiSquared(tc.p, tc.d), tc.want, 1e-6)
		})
	}
}

func TestCalculateKChiSquared_TwoDofClosedForm(t *testing.T) {
	for _, p := range []float64{0.1, 0.5, 0.9, 0.99} {
		want := -2 * math.Log(1-p)
		assertNear(t, "quantile", calculateKChiSquared(p, 2), want, 1e-6)
	}
}

func TestCalculateKChiSquared_MonotonicInProbability(t *testing.T) {
	prev := calculateKChiSquared(0.1, 3)
	for _, p := range []float64{0.5, 0.9, 0.99} {
		cur := calculateKChiSquared(p, 3)
		if cur <= prev {
			t.Errorf("quantile(%v) = %v is not greater than previous %v", p, cur, prev)
		}
		prev = cur
	}
}

func TestNormalizeAngle(t *testing.T) {
	tests := []struct {
		name  string
		angle float64
		want  float64
	}{
		{"zero", 0, 0},
		{"already normalised positive", math.Pi / 2, math.Pi / 2},
		{"already normalised negative", -math.Pi / 2, -math.Pi / 2},
		{"3pi/2 wraps to -pi/2", 3 * math.Pi / 2, -math.Pi / 2},
		{"-3pi/2 wraps to pi/2", -3 * math.Pi / 2, math.Pi / 2},
		{"full turn", 2 * math.Pi, 0},
		{"two turns plus offset", 4*math.Pi + 0.3, 0.3},
		{"just past pi", math.Pi + 0.1, -math.Pi + 0.1},
		{"just before -pi", -math.Pi - 0.1, math.Pi - 0.1},
	}
	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			assertNear(t, "normalizeAngle", normalizeAngle(tc.angle), tc.want, 1e-12)
		})
	}
}

func TestNormalizeAngle_AlwaysInRange(t *testing.T) {
	for angle := -20.0; angle <= 20.0; angle += 0.37 {
		got := normalizeAngle(angle)
		if got < -math.Pi-1e-12 || got > math.Pi+1e-12 {
			t.Errorf("normalizeAngle(%v) = %v is outside [-pi, pi]", angle, got)
		}
	}
}

func TestSubMatrix(t *testing.T) {
	src := mat.NewSymDense(3, []float64{
		1, 2, 3,
		2, 4, 5,
		3, 5, 6,
	})

	tests := []struct {
		name string
		keep []int
		want [][]float64
	}{
		{"drop middle", []int{0, 2}, [][]float64{{1, 3}, {3, 6}}},
		{"drop first", []int{1, 2}, [][]float64{{4, 5}, {5, 6}}},
		{"single element", []int{1}, [][]float64{{4}}},
		{"keep everything", []int{0, 1, 2}, [][]float64{{1, 2, 3}, {2, 4, 5}, {3, 5, 6}}},
		{"reordered indices", []int{2, 0}, [][]float64{{6, 3}, {3, 1}}},
	}
	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			got := subMatrix(src, tc.keep)
			assertMatrix(t, "subMatrix", got, tc.want, 0)
		})
	}
}

func TestSubMatrix_DoesNotModifySource(t *testing.T) {
	src := mat.NewSymDense(2, []float64{1, 2, 2, 3})
	_ = subMatrix(src, []int{1})
	assertMatrix(t, "source", src, [][]float64{{1, 2}, {2, 3}}, 0)
}

func TestRtoaConversion(t *testing.T) {
	tests := []struct {
		reported float64
		factor   float64
		value    float64
	}{
		{reported: 1.0, factor: 0.0, value: -985024.0},
		{reported: 985025.0, factor: 0.0, value: 1.0},
		{reported: 1970048.0, factor: 0.0, value: 985024.0},
		{reported: 1.0, factor: 1.0, value: -985022.0},
		{reported: 492513.0, factor: 1.0, value: 2.0},
		{reported: 985024.0, factor: 1.0, value: 985024.0},
		{reported: 1.0, factor: 2.0, value: -985020.0},
		{reported: 246257.0, factor: 2.0, value: 4.0},
		{reported: 492512.0, factor: 2.0, value: 985024.0},
		{reported: 1.0, factor: 3.0, value: -985024.0},
		{reported: 123129.0, factor: 3.0, value: 8.0},
		{reported: 246256.0, factor: 3.0, value: 985024.0},
		{reported: 1.0, factor: 4.0, value: -985024.0},
		{reported: 61565.0, factor: 4.0, value: 16.0},
		{reported: 123128.0, factor: 4.0, value: 985024.0},
		{reported: 1.0, factor: 5.0, value: -985024.0},
		{reported: 30783.0, factor: 5.0, value: 32.0},
		{reported: 61564.0, factor: 5.0, value: 985024.0},
	}
	for _, tt := range tests {
		t.Run("RTOA", func(t *testing.T) {
			assertNear(t, "", rtoa3gppToSec(tt.reported, tt.factor), tt.value*BASE_TIME_UNIT, 1e-3)
		})
	}
}

func TestBuildNormalEquation_NoCovariance(t *testing.T) {
	j := [][]float64{
		{1, 0},
		{0, 1},
		{1, 1},
	}
	r := []float64{1, 2, 3}

	h, g, _, err := buildNormalEquation(
		newVec(0, 0), nil, nil, []int{0, 1, 2},
		[]MeasurementBlock{fixedBlock(j, r)},
		[]*mat.SymDense{nil},
		2,
	)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}

	assertMatrix(t, "H", h, [][]float64{{2, 1}, {1, 2}}, 1e-12)
	assertVec(t, "g", g, []float64{4, 5}, 1e-12)
}

func TestBuildNormalEquation_WithDiagonalCovariance(t *testing.T) {
	j := [][]float64{
		{1, 0},
		{0, 1},
		{1, 1},
	}
	r := []float64{1, 2, 3}
	cov := mat.NewSymDense(3, []float64{
		2, 0, 0,
		0, 4, 0,
		0, 0, 1,
	})

	h, g, _, err := buildNormalEquation(
		newVec(0, 0), nil, nil, []int{0, 1, 2},
		[]MeasurementBlock{fixedBlock(j, r)},
		[]*mat.SymDense{cov},
		2,
	)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}

	assertMatrix(t, "H", h, [][]float64{{1.5, 1}, {1, 1.25}}, 1e-12)
	assertVec(t, "g", g, []float64{3.5, 3.5}, 1e-12)
}

func TestBuildNormalEquation_AccumulatesAcrossBlocks(t *testing.T) {
	j := [][]float64{{1, 0}, {0, 1}, {1, 1}}
	r := []float64{1, 2, 3}

	h, g, _, err := buildNormalEquation(
		newVec(0, 0), nil, nil, []int{0, 1, 2},
		[]MeasurementBlock{fixedBlock(j, r), fixedBlock(j, r)},
		[]*mat.SymDense{nil, nil},
		2,
	)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	assertMatrix(t, "H", h, [][]float64{{4, 2}, {2, 4}}, 1e-12)
	assertVec(t, "g", g, []float64{8, 10}, 1e-12)
}

func TestBuildNormalEquation_CovarianceFollowsDroppedTrp(t *testing.T) {
	block := &mockBlock{
		jr: func(_ *mat.VecDense, ctx *PositioningContext) ([][]float64, []float64) {
			ctx.KeptTrp = []int{0, 2}
			return identitySlice(2), []float64{1, 2}
		},
	}
	cov := mat.NewSymDense(3, []float64{
		2, 0, 0,
		0, 4, 0,
		0, 0, 1,
	})

	h, g, _, err := buildNormalEquation(
		newVec(0, 0), nil, nil, []int{0, 1, 2},
		[]MeasurementBlock{block},
		[]*mat.SymDense{cov},
		2,
	)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}

	assertMatrix(t, "H", h, [][]float64{{0.5, 0}, {0, 1}}, 1e-12)
	assertVec(t, "g", g, []float64{0.5, 2}, 1e-12)
}

func TestBuildNormalEquation_Errors(t *testing.T) {
	okJ := [][]float64{
		{1, 0},
		{0, 1},
		{1, 1},
	}
	okR := []float64{1, 2, 3}

	tests := []struct {
		name    string
		blocks  []MeasurementBlock
		covs    []*mat.SymDense
		wantMsg string
	}{
		{
			name: "block / covariance count mismatch",
			blocks: []MeasurementBlock{
				fixedBlock(okJ, okR),
			},
			covs:    []*mat.SymDense{nil, nil},
			wantMsg: "does not match",
		},
		{
			name: "Jacobian has wrong number of columns",
			blocks: []MeasurementBlock{
				fixedBlock(
					[][]float64{
						{1, 0, 0},
						{0, 1, 0},
						{1, 1, 0},
					},
					okR,
				),
			},
			covs:    []*mat.SymDense{nil},
			wantMsg: "invalid Jacobian dimension",
		},
		{
			name: "residual length differs from Jacobian rows",
			blocks: []MeasurementBlock{
				fixedBlock(
					okJ,
					[]float64{1, 2},
				),
			},
			covs:    []*mat.SymDense{nil},
			wantMsg: "invalid residual dimension",
		},
		{
			name: "covariance has wrong dimension",
			blocks: []MeasurementBlock{
				fixedBlock(okJ, okR),
			},
			covs: []*mat.SymDense{
				mat.NewSymDense(2, []float64{
					1, 0,
					0, 1,
				}),
			},
			wantMsg: "rank deficient",
		},
		{
			name: "covariance is not positive definite",
			blocks: []MeasurementBlock{
				fixedBlock(okJ, okR),
			},
			covs: []*mat.SymDense{
				mat.NewSymDense(3, []float64{
					1, 0, 0,
					0, 0, 0,
					0, 0, 1,
				}),
			},
			wantMsg: "rank deficient",
		},
	}

	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			h, g, _, err := buildNormalEquation(
				newVec(0, 0),
				nil,
				nil,
				[]int{0, 1, 2},
				tc.blocks,
				tc.covs,
				2,
			)

			if err == nil {
				t.Fatalf(
					"expected an error containing %q, got nil (H=%v, g=%v)",
					tc.wantMsg,
					h,
					g,
				)
			}

			if !strings.Contains(err.Error(), tc.wantMsg) {
				t.Errorf(
					"error %q does not contain %q",
					err.Error(),
					tc.wantMsg,
				)
			}

			if h != nil || g != nil {
				t.Errorf(
					"expected nil results alongside the error, got H=%v g=%v",
					h,
					g,
				)
			}
		})
	}
}
