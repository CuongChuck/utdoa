package utdoa

import (
	"strings"
	"testing"
)

func TestInitUePos_IsCentroidOfTrps(t *testing.T) {
	t.Run("2D", func(t *testing.T) {
		pos, err := initUePos(newTrps([]float64{0, 0}, []float64{6, 0}, []float64{0, 3}), 2)
		if err != nil {
			t.Fatalf("unexpected error: %v", err)
		}
		assertVec(t, "pos", pos, []float64{2, 1}, 1e-12)
	})

	t.Run("3D", func(t *testing.T) {
		pos, err := initUePos(newTrps([]float64{0, 0, 0}, []float64{2, 4, 6}), 3)
		if err != nil {
			t.Fatalf("unexpected error: %v", err)
		}
		assertVec(t, "pos", pos, []float64{1, 2, 3}, 1e-12)
	})
}

func TestInitUePos_NeedsAtLeastTwoTrps(t *testing.T) {
	for _, n := range []int{0, 1} {
		positions := make([][]float64, n)
		for i := range positions {
			positions[i] = []float64{1, 1}
		}
		pos, err := initUePos(newTrps(positions...), 2)
		if err == nil {
			t.Errorf("%d TRP(s): expected an error, got %v", n, pos)
		}
	}
}

func TestDefaultUtdoaConfig(t *testing.T) {
	cfg := defaultUtdoaConfig()

	if cfg.Dimension != 2 {
		t.Errorf("Dimension = %d, want 2", cfg.Dimension)
	}
	if len(cfg.Measurements) != 1 {
		t.Fatalf("len(Measurements) = %d, want 1", len(cfg.Measurements))
	}
	if _, ok := cfg.Measurements[0].(*TdoaMeasurement); !ok {
		t.Errorf("default measurement is %T, want *TdoaMeasurement", cfg.Measurements[0])
	}
}

func TestDefaultUtdoaConfig_ReturnsIndependentInstances(t *testing.T) {
	a, b := defaultUtdoaConfig(), defaultUtdoaConfig()
	if a.Measurements[0] == b.Measurements[0] {
		t.Error("two default configs share the same measurement block")
	}
}

func TestWithDimension(t *testing.T) {
	for _, dim := range []int{2, 3, 7} {
		cfg := defaultUtdoaConfig()
		WithDimension(dim)(&cfg)
		if cfg.Dimension != dim {
			t.Errorf("Dimension = %d, want %d", cfg.Dimension, dim)
		}
	}
}

func TestWith3D(t *testing.T) {
	cfg := defaultUtdoaConfig()
	With3D()(&cfg)
	if cfg.Dimension != 3 {
		t.Errorf("Dimension = %d, want 3", cfg.Dimension)
	}
}

func TestWithAzimuthAndZenith_AppendInOrder(t *testing.T) {
	cfg := defaultUtdoaConfig()
	WithAzimuth()(&cfg)
	WithZenith()(&cfg)

	if len(cfg.Measurements) != 3 {
		t.Fatalf("len(Measurements) = %d, want 3", len(cfg.Measurements))
	}
	if _, ok := cfg.Measurements[0].(*TdoaMeasurement); !ok {
		t.Errorf("Measurements[0] is %T, want *TdoaMeasurement", cfg.Measurements[0])
	}
	if _, ok := cfg.Measurements[1].(*AzimuthMeasurement); !ok {
		t.Errorf("Measurements[1] is %T, want *AzimuthMeasurement", cfg.Measurements[1])
	}
	if _, ok := cfg.Measurements[2].(*ZenithMeasurement); !ok {
		t.Errorf("Measurements[2] is %T, want *ZenithMeasurement", cfg.Measurements[2])
	}
}

func TestOptions_LastDimensionWins(t *testing.T) {
	cfg := defaultUtdoaConfig()
	for _, opt := range []UtdoaOption{With3D(), WithDimension(2)} {
		opt(&cfg)
	}
	if cfg.Dimension != 2 {
		t.Errorf("Dimension = %d, want 2", cfg.Dimension)
	}
}

func TestUtdoaPositioning_InvalidArguments(t *testing.T) {
	validSet := func() []TrpMeas {
		return synthMeas([]float64{30, 40}, square2D, rsrp2D, 1)
	}

	tests := []struct {
		name    string
		measSet []TrpMeas
		opts    []UtdoaOption
		wantMsg string
	}{
		{"dimension 0", validSet(), []UtdoaOption{WithDimension(0)}, "dimension must be 2 or 3"},
		{"dimension 1", validSet(), []UtdoaOption{WithDimension(1)}, "dimension must be 2 or 3"},
		{"dimension 4", validSet(), []UtdoaOption{WithDimension(4)}, "dimension must be 2 or 3"},
		{"zenith in 2D", validSet(), []UtdoaOption{WithZenith()}, "zenith"},
		{"zenith in 2D regardless of option order", validSet(), []UtdoaOption{WithZenith(), WithDimension(2)}, "zenith"},
		{"no measurements", nil, nil, "not enough measurements"},
		{"single TRP", validSet()[:1], nil, "not enough measurements"},
	}
	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			shape, err := UtdoaPositioning(tc.measSet, tc.opts...)
			if err == nil {
				t.Fatalf("expected an error, got %#v", shape)
			}
			if shape != nil {
				t.Errorf("shape should be nil on error, got %#v", shape)
			}
			if !strings.Contains(err.Error(), tc.wantMsg) {
				t.Errorf("error %q does not contain %q", err, tc.wantMsg)
			}
		})
	}
}

func TestUtdoaPositioning_2D_CommonClockOffsetCancels(t *testing.T) {
	truth := []float64{30, 40}

	base := synthMeas(truth, square2D, rsrp2D, 1)
	shifted := synthMeas(truth, square2D, rsrp2D, 1)
	for i := range shifted {
		shifted[i].rtoa += 5000
	}

	a, errA := UtdoaPositioning(base)
	b, errB := UtdoaPositioning(shifted)
	if errA != nil || errB != nil {
		t.Fatalf("unexpected errors: %v / %v", errA, errB)
	}

	pa := a.(*UncertaintyEllipse2D).Point
	pb := b.(*UncertaintyEllipse2D).Point
	assertVec(t, "Point", pb, []float64{pa.AtVec(0), pa.AtVec(1)}, 1e-6)
}

func TestUtdoaPositioning_2D_UncertaintyGrowsWithMeasurementNoise(t *testing.T) {
	truth := []float64{30, 40}
	sharp, err := UtdoaPositioning(synthMeas(truth, square2D, rsrp2D, 1))
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	blurry, err := UtdoaPositioning(synthMeas(truth, square2D, rsrp2D, 10))
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}

	s := sharp.(*UncertaintyEllipse2D)
	b := blurry.(*UncertaintyEllipse2D)

	assertNear(t, "major axis ratio", b.MajorAxis/s.MajorAxis, 10, 1e-2)
	assertNear(t, "minor axis ratio", b.MinorAxis/s.MinorAxis, 10, 1e-2)
}
