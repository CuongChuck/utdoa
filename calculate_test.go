package utdoa

import (
	"encoding/csv"
	"encoding/json"
	"fmt"
	"io"
	"log"
	"math"
	"os"
	"slices"
	"strconv"
	"strings"
	"testing"

	"gonum.org/v1/gonum/mat"
)

func TestGaussNewton_PositionInputs(t *testing.T) {
	file, err := os.Open("failed_testcase.csv")
	if err != nil {
		log.Fatalf("Failed to open testcase file: %v", err)
	}
	defer file.Close()

	reader := csv.NewReader(file)
	_, err = reader.Read()
	if err != nil {
		log.Fatalf("Error reading first row: %v", err)
	}
	for {
		options := []UtdoaOption{}
		var (
			targetUeGcs    []float64
			trpGcs         [][]float64
			tdoaQuality    [][]float64
			azimuthQuality [][]float64
			zenithQuality  [][]float64
		)
		row, err := reader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			log.Fatalf("Error reading row: %v", err)
		}

		testIdCell := row[0]
		dimCell := row[1]
		targetUeGcsCell := row[2]
		trpGcsCell := row[3]
		tdoaQualCell := row[4]
		azimuthQualCell := row[5]
		zenithQualCell := row[6]
		factorCell := row[7]
		uePosErrCell := row[8]
		testId, err := strconv.Atoi(testIdCell)
		if err != nil {
			log.Fatalf("Failed to parse test id %v: %v", testIdCell, err)
		}
		dimension, err := strconv.Atoi(dimCell)
		if err != nil {
			log.Fatalf("Failed to parse location dimension %v: %v", dimCell, err)
		}
		options = append(options, WithDimension(dimension))
		err = json.Unmarshal([]byte(targetUeGcsCell), &targetUeGcs)
		if err != nil {
			log.Fatalf("Failed to parse UE coordinates %v: %v", targetUeGcsCell, err)
		}
		err = json.Unmarshal([]byte(trpGcsCell), &trpGcs)
		if err != nil {
			log.Fatalf("Failed to parse TRP coordinates %v: %v", trpGcsCell, err)
		}
		err = json.Unmarshal([]byte(tdoaQualCell), &tdoaQuality)
		if err != nil {
			log.Fatalf("Failed to parse TDOA measurement quality %v: %v", tdoaQualCell, err)
		}
		if azimuthQualCell != "" {
			err = json.Unmarshal([]byte(azimuthQualCell), &azimuthQuality)
			if err != nil {
				log.Fatalf(
					"Failed to parse azimuth measurement quality %v: %v",
					azimuthQualCell,
					err,
				)
			}
			options = append(options, WithAzimuth())
		}
		if zenithQualCell != "" {
			err = json.Unmarshal([]byte(zenithQualCell), &zenithQuality)
			if err != nil {
				log.Fatalf(
					"Failed to parse zenith measurement quality %v: %v",
					zenithQualCell,
					err,
				)
			}
			options = append(options, WithZenith())
		}
		factor, err := strconv.ParseFloat(factorCell, 64)
		if err != nil {
			log.Fatalf("Failed to parse time granularity factor %v: %v", factorCell, err)
		}
		uePosErr, err := strconv.ParseFloat(uePosErrCell, 64)
		if err != nil {
			log.Fatalf("Failed to parse UE position error %v: %v", factorCell, err)
		}
		fmt.Printf("\nTestcase %v:\n", testId)
		targetUePcs := gcsToPcsSlice(targetUeGcs, trpGcs[0])
		numTrp := len(trpGcs)
		measSet := make([]TrpMeas, numTrp)
		for idx := range measSet {
			trpPosPcs := gcsToPcsSlice(trpGcs[idx], trpGcs[0])
			measSet[idx].id = idx
			measSet[idx].isHighAcc = idx%2 == 0
			measSet[idx].pos3gpp = gcsTo3gpp(trpGcs[idx], measSet[idx].isHighAcc)
			measSet[idx].factor = factor
			measSet[idx].rtoa3gpp = gcsToRtoa(trpPosPcs, targetUePcs, factor)
			measSet[idx].tdoaQuality = MeasQuality{
				quality: tdoaQuality[idx][0],
				res:     tdoaQuality[idx][1],
			}
			if azimuthQualCell != "" {
				measSet[idx].azimuth = pcsToAzimuth3gpp(targetUePcs, trpPosPcs)
				measSet[idx].azimuthQuality = MeasQuality{
					quality: azimuthQuality[idx][0],
					res:     azimuthQuality[idx][1],
				}
			}
			if zenithQualCell != "" && dimension == 3 {
				measSet[idx].zenith = pcsToZenith3gpp(targetUePcs, trpPosPcs)
				measSet[idx].zenithQuality = MeasQuality{
					quality: zenithQuality[idx][0],
					res:     zenithQuality[idx][1],
				}
			}
			fmt.Printf("TRP %v - test TDOA: %v\n", idx, distSlice([]float64{0.0, 0.0}, targetUePcs)-distSlice(targetUePcs, trpPosPcs))
		}
		result, err := UtdoaPositioning(measSet, options...)
		if err != nil {
			t.Errorf("\nTestcase %v failed: Error in positioning: %v", testId, err)
			continue
		}
		pointPcs, _, _ := result.GetFigure()
		pointCoord := pcsToGcsSlice(pointPcs, trpGcs[0])
		targetUeToTrp := make([]float64, numTrp)
		for idx := range targetUeToTrp {
			targetUeToTrp[idx] = distSlice(targetUePcs, measSet[idx].pos.RawVector().Data)
		}
		calculatedUeToTrp := make([]float64, numTrp)
		for idx := range calculatedUeToTrp {
			calculatedUeToTrp[idx] = distSlice(pointPcs, measSet[idx].pos.RawVector().Data)
		}
		posDiff := distSlice(gcsToPcsSlice(targetUeGcs, trpGcs[0]), pointPcs)
		if posDiff > uePosErr {
			t.Errorf(
				"\nTestcase %v failed:"+
					"\nGCS: got %v, want %v, diff %v"+
					"\nPCS (meters): got %v, want %v"+
					"\ntarget UE to TRPs %v\ncalculated UE to TRPs %v",
				testId,
				pointCoord, targetUeGcs, posDiff,
				pointPcs, targetUePcs,
				targetUeToTrp, calculatedUeToTrp,
			)
		} else {
			fmt.Printf(
				"\nTestcase %v passed:"+
					"\nGCS: got %v, want %v, diff %v"+
					"\nPCS (meters): got %v, want %v"+
					"\ntarget UE to TRPs %v\ncalculated UE to TRPs %v\n",
				testId,
				pointCoord, targetUeGcs, posDiff,
				pointPcs, targetUePcs,
				targetUeToTrp, calculatedUeToTrp)
		}
	}
}

func TestGaussNewton_ConvergesOnLinearProblem(t *testing.T) {
	trps := newTrps([]float64{0, 0}, []float64{10, 0})
	ue0 := newVec(5, 5)

	res, err := gaussNewton(
		ue0, trps,
		[]MeasurementBlock{pullToward(3, -2)},
		[]*mat.SymDense{nil},
		2,
	)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}

	if !res.IsConverged {
		t.Error("IsConverged = false, want true")
	}
	assertVec(t, "UePos", res.UePos, []float64{3, -2}, 1e-6)
	assertMatrix(t, "CovarMat", res.CovarMat, [][]float64{{1, 0}, {0, 1}}, 1e-12)
}

func TestGaussNewton_WeightingDoesNotBiasExactFit(t *testing.T) {
	trps := newTrps([]float64{0, 0}, []float64{10, 0})
	cov := mat.NewSymDense(2, []float64{4, 0, 0, 9})

	res, err := gaussNewton(
		newVec(5, 5), trps,
		[]MeasurementBlock{pullToward(3, -2)},
		[]*mat.SymDense{cov},
		2,
	)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}

	if !res.IsConverged {
		t.Error("IsConverged = false, want true")
	}
	assertVec(t, "UePos", res.UePos, []float64{3, -2}, 1e-6)
	assertMatrix(t, "CovarMat", res.CovarMat, [][]float64{{4, 0}, {0, 9.0}}, 1e-12)
}

func TestGaussNewton_CombinesBlocksIntoLeastSquaresSolution(t *testing.T) {
	trps := newTrps([]float64{0, 0}, []float64{10, 0})

	res, err := gaussNewton(
		newVec(5, 5), trps,
		[]MeasurementBlock{pullToward(0, 0), pullToward(4, 2)},
		[]*mat.SymDense{nil, nil},
		2,
	)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}

	if !res.IsConverged {
		t.Error("IsConverged = false, want true")
	}
	assertVec(t, "UePos", res.UePos, []float64{2, 1}, 1e-6)
	assertMatrix(t, "CovarMat", res.CovarMat, [][]float64{{1.0 / 2.0, 0}, {0, 1.0 / 2.0}}, 1e-12)
}

func TestGaussNewton_StartingAtSolutionConvergesImmediately(t *testing.T) {
	trps := newTrps([]float64{0, 0}, []float64{10, 0})

	res, err := gaussNewton(
		newVec(3, -2), trps,
		[]MeasurementBlock{pullToward(3, -2)},
		[]*mat.SymDense{nil},
		2,
	)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if !res.IsConverged {
		t.Error("IsConverged = false, want true")
	}
	assertVec(t, "UePos", res.UePos, []float64{3, -2}, 1e-9)
}

func TestGaussNewton_StopsAtMaxLoopWhenNotConverging(t *testing.T) {
	trps := newTrps([]float64{0, 0}, []float64{10, 0})
	block := fixedBlock(identitySlice(2), []float64{1, 1})

	res, err := gaussNewton(
		newVec(5, 5), trps,
		[]MeasurementBlock{block},
		[]*mat.SymDense{nil},
		2,
	)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}

	if res.IsConverged {
		t.Error("IsConverged = true, want false")
	}
	if block.calls != MAX_LOOP {
		t.Errorf("Jacobian evaluated %d times, want %d", block.calls, MAX_LOOP)
	}
	if res.CovarMat == nil {
		t.Error("CovarMat is nil, want the last normal matrix")
	}
}

func TestGaussNewton_DoesNotModifyInitialPosition(t *testing.T) {
	trps := newTrps([]float64{0, 0}, []float64{10, 0})
	ue0 := newVec(5, 5)

	res, err := gaussNewton(
		ue0, trps,
		[]MeasurementBlock{fixedBlock(identitySlice(2), []float64{1, 1})},
		[]*mat.SymDense{nil},
		2,
	)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}

	assertVec(t, "initial position", ue0, []float64{5, 5}, 0)
	if res.UePos == ue0 {
		t.Error("result aliases the caller's initial position vector")
	}
}

func TestGaussNewton_ExcludesTrpCoincidentWithUe(t *testing.T) {
	trps := newTrps([]float64{0, 0}, []float64{10, 0})
	block := fixedBlock(identitySlice(2), []float64{1, 1})

	_, err := gaussNewton(
		newVec(0, 0), trps,
		[]MeasurementBlock{block},
		[]*mat.SymDense{nil},
		2,
	)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}

	if len(block.seenKept) == 0 {
		t.Fatal("block was never evaluated")
	}
	if got := block.seenKept[0]; !slices.Equal(got, []int{1}) {
		t.Errorf("first iteration kept TRPs %v, want [1]", got)
	}
}

func TestGaussNewton_Errors(t *testing.T) {
	trps := newTrps([]float64{0, 0}, []float64{10, 0})

	t.Run("block and covariance counts differ", func(t *testing.T) {
		res, err := gaussNewton(
			newVec(5, 5), trps,
			[]MeasurementBlock{pullToward(1, 1)},
			[]*mat.SymDense{nil, nil},
			2,
		)
		if err == nil {
			t.Fatalf("expected an error, got %+v", res)
		}
		if res != nil {
			t.Error("result should be nil when an error is returned")
		}
	})

	t.Run("singular normal matrix", func(t *testing.T) {
		zeroJ := fixedBlock([][]float64{{0, 0}, {0, 0}}, []float64{0, 0})
		res, err := gaussNewton(
			newVec(5, 5), trps,
			[]MeasurementBlock{zeroJ},
			[]*mat.SymDense{nil},
			2,
		)
		if err == nil {
			t.Fatalf("expected an error, got %+v", res)
		}
		if !strings.Contains(err.Error(), "rank deficient") {
			t.Errorf("error %q does not mention rank deficiency", err)
		}
	})
}

func TestUncertaintyCovariance(t *testing.T) {
	dim := 2
	targetUePos := []float64{30, 40}
	measSet := synthMeas(targetUePos, square2D, rsrp2D, 5.0)
	tdoaMeas := &TdoaMeasurement{}
	tdoaMeas.SetDim(dim)
	tdoaMeas.Calculate(measSet)
	calculateTrpPos(measSet, 2)
	measurements := []MeasurementBlock{tdoaMeas}
	covarianceMats := []*mat.SymDense{tdoaMeas.CovarianceMatrix(measSet)}
	initialUePos, err := initUePos(measSet, 2)
	if err != nil {
		t.Errorf("Testing calculating initial UE position failed: %v", err)
	}

	gnResult, err := gaussNewton(
		initialUePos,
		measSet,
		measurements,
		covarianceMats,
		dim,
	)
	if err != nil {
		t.Fatalf("gaussNewton failed: %v", err)
	}

	numTrp := len(measSet)
	distTrps := make([]float64, numTrp)
	keep := make([]int, 0, numTrp)
	for i := range numTrp {
		sqrDist := squaredDist(gnResult.UePos, measSet[i].pos)
		if sqrDist > MIN_DIST2 {
			distTrps[i] = math.Sqrt(sqrDist)
			keep = append(keep, i)
		}
	}

	coefficient, _, _, err := buildNormalEquation(
		gnResult.UePos,
		measSet,
		distTrps,
		keep,
		measurements,
		covarianceMats,
		dim,
	)
	if err != nil {
		t.Fatalf("buildNormalEquation failed: %v", err)
	}

	var chol mat.Cholesky
	if ok := chol.Factorize(coefficient); !ok {
		t.Fatalf("coefficient matrix is rank deficient")
	}
	expectedCov := mat.NewSymDense(dim, nil)
	if err := chol.InverseTo(expectedCov); err != nil {
		t.Fatalf("failed to invert coefficient matrix: %v", err)
	}

	assertMatricesNear(t, "Covariance Matrix", gnResult.CovarMat, expectedCov, 1e-6)
}

func TestUncertaintyEllipse2D_Figure(t *testing.T) {
	tests := []struct {
		name       string
		cov        []float64
		eigenMajor float64
		eigenMinor float64
		angle      float64
	}{
		{"axis aligned, x dominant", []float64{4, 0, 0, 1}, 4, 1, 0},
		{"axis aligned, y dominant", []float64{1, 0, 0, 4}, 4, 1, math.Pi / 2},
		{"rotated by 45 degrees", []float64{2, 1, 1, 2}, 3, 1, math.Pi / 4},
		{"isotropic", []float64{1, 0, 0, 1}, 1, 1, 0},
		{"negative eigenvalue is clamped to zero", []float64{1, 2, 2, 1}, 3, 0, math.Pi / 4},
	}
	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			var e UncertaintyEllipse2D
			err := e.Figure(mat.NewSymDense(2, tc.cov), nil)
			if err != nil {
				t.Fatalf("unexpected error: %v", err)
			}

			assertNear(t, "MajorAxis", e.MajorAxis, math.Sqrt(chi2P90Dof2*tc.eigenMajor), 1e-6)
			assertNear(t, "MinorAxis", e.MinorAxis, math.Sqrt(chi2P90Dof2*tc.eigenMinor), 1e-6)
			assertNear(t, "Angle", e.Angle, tc.angle, 1e-12)
			if e.MajorAxis < e.MinorAxis {
				t.Errorf("MajorAxis %v < MinorAxis %v", e.MajorAxis, e.MinorAxis)
			}
		})
	}
}

func TestUncertaintyEllipse2D_Figure_KeepsPoint(t *testing.T) {
	e := UncertaintyEllipse2D{Point: newVec(1, 2)}
	if err := e.Figure(mat.NewSymDense(2, []float64{4, 0, 0, 1}), nil); err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	assertVec(t, "Point", e.Point, []float64{1, 2}, 0)
}

func TestUncertaintyEllipse2D_Figure_ScalesWithSqrtOfCovariance(t *testing.T) {
	var small, big UncertaintyEllipse2D
	_ = small.Figure(mat.NewSymDense(2, []float64{2, 0.5, 0.5, 1}), nil)
	_ = big.Figure(mat.NewSymDense(2, []float64{8, 2, 2, 4}), nil)

	assertNear(t, "major ratio", big.MajorAxis/small.MajorAxis, 2, 1e-9)
	assertNear(t, "minor ratio", big.MinorAxis/small.MinorAxis, 2, 1e-9)
	assertNear(t, "angle", big.Angle, small.Angle, 1e-12)
}

func TestUncertaintyEllipsoid3D_Figure(t *testing.T) {
	tests := []struct {
		name  string
		cov   []float64
		eigen [3]float64
		angle float64
	}{
		{
			"diagonal, x largest",
			[]float64{9, 0, 0, 0, 4, 0, 0, 0, 1},
			[3]float64{9, 4, 1}, 0,
		},
		{
			"diagonal, z largest",
			[]float64{1, 0, 0, 0, 4, 0, 0, 0, 9},
			[3]float64{9, 4, 1}, math.Pi / 2,
		},
		{
			"correlated x/y with eigenvalues 1, 2, 3",
			[]float64{2, 1, 0, 1, 2, 0, 0, 0, 2},
			[3]float64{3, 2, 1}, math.Pi / 4,
		},
	}
	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			var e UncertaintyEllipsoid3D
			if err := e.Figure(mat.NewSymDense(3, tc.cov), nil); err != nil {
				t.Fatalf("unexpected error: %v", err)
			}

			assertNear(t, "MajorAxis", e.MajorAxis, math.Sqrt(chi2P90Dof3*tc.eigen[0]), 1e-6)
			assertNear(t, "MidAxis", e.MidAxis, math.Sqrt(chi2P90Dof3*tc.eigen[1]), 1e-6)
			assertNear(t, "MinorAxis", e.MinorAxis, math.Sqrt(chi2P90Dof3*tc.eigen[2]), 1e-6)
			assertNear(t, "Angle", e.Angle, tc.angle, 1e-12)

			if !(e.MajorAxis >= e.MidAxis && e.MidAxis >= e.MinorAxis) {
				t.Errorf("axes not ordered: major=%v mid=%v minor=%v",
					e.MajorAxis, e.MidAxis, e.MinorAxis)
			}
		})
	}
}

func TestUncertaintyEllipsoid3D_Figure_NegativeEigenvalueClampedToZero(t *testing.T) {
	cov := mat.NewSymDense(3, []float64{
		1, 2, 0,
		2, 1, 0,
		0, 0, 2,
	})
	var e UncertaintyEllipsoid3D
	if err := e.Figure(cov, nil); err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	assertNear(t, "MinorAxis", e.MinorAxis, 0, 1e-12)
	assertNear(t, "MajorAxis", e.MajorAxis, math.Sqrt(chi2P90Dof3*3), 1e-6)
}
