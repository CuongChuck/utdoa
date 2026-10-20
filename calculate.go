package utdoa

import (
	"errors"
	"fmt"
	"math"

	"gonum.org/v1/gonum/mat"
)

/*
 * Description: Gauss-Newton algorithm/loop to compute UE estimate position.
 * Parameter : uePos
 * Type : *mat.VecDense
 * Parameter : measSet
 * Type : []TrpMeas
 * Parameter : measurements
 * Type : []MeasurementBlock
 * Parameter : covarianceMats
 * Type : []*mat.SymDense
 * Parameter : dim
 * Type : int
 */
func gaussNewton(
	uePos *mat.VecDense,
	measSet []TrpMeas,
	measurements []MeasurementBlock,
	covarianceMats []*mat.SymDense,
	dim int,
) (*GaussNewtonResult, error) {
	numTrp := len(measSet)
	ueCurrPos := mat.NewVecDense(dim, nil)
	ueCurrPos.CopyVec(uePos)
	distTrps := make([]float64, numTrp)
	var covarMat *mat.SymDense
	isConverged := false

	evalAt := func(pos *mat.VecDense) (*mat.SymDense, *mat.VecDense, float64, error) {
		var keep []int
		for i := range numTrp {
			sq := squaredDist(pos, measSet[i].pos)
			if sq > MIN_DIST2 {
				distTrps[i] = math.Sqrt(sq)
				keep = append(keep, i)
			}
		}
		return buildNormalEquation(pos, measSet, distTrps, keep,
			measurements, covarianceMats, dim)
	}

	for loop := range MAX_LOOP {
		fmt.Printf("UE:\n%v\n", mat.Formatted(ueCurrPos))
		var (
			delta mat.VecDense
			chol  mat.Cholesky
		)
		coefficient, constant, cost, err := evalAt(ueCurrPos)
		if err != nil {
			return nil, err
		}

		if isConverged || loop == MAX_LOOP-1 {
			if ok := chol.Factorize(coefficient); !ok {
				return nil, errors.New("normal matrix is rank deficient")
			}
			uncertaintyCovariance := mat.NewSymDense(coefficient.SymmetricDim(), nil)
			if err := chol.InverseTo(uncertaintyCovariance); err != nil {
				return nil, err
			}
			covarMat = uncertaintyCovariance
			fmt.Printf("Residual block:\n%v\n", mat.Formatted(constant))
			break
		}

		if ok := chol.Factorize(coefficient); !ok {
			return nil, errors.New("normal matrix is rank deficient")
		}

		if err := chol.SolveVecTo(&delta, constant); err != nil {
			return nil, err
		}

		delta.ScaleVec(-1.0, &delta)
		trial := mat.NewVecDense(dim, nil)
		step := 1.0
		accepted := false
		for ; step >= 1e-6; step /= 2 {
			trial.AddScaledVec(ueCurrPos, step, &delta)
			_, _, trialCost, err := evalAt(trial)
			if err == nil && trialCost <= cost {
				accepted = true
				break
			}
		}
		if !accepted {
			// No descent step exists: we are at a stationary point.
			isConverged = true
			continue
		}

		stepNorm := step * mat.Norm(&delta, 2)
		if stepNorm < MIN_DIST_DIFF {
			isConverged = true
		}
		ueCurrPos.CopyVec(trial)
	}

	return &GaussNewtonResult{
		ueCurrPos,
		covarMat,
		isConverged,
	}, nil
}

/*
 * Description: Computes the uncertainty ellipse around the UE position estimate.
 * Parameter : ellipse
 * Type : UncertaintyEllipse2D
 * Parameter : covariance
 * Type : *mat.SymDense
 */
func (ellipse *UncertaintyEllipse2D) Figure(
	covariance *mat.SymDense,
	origin *mat.VecDense,
) error {
	varianceX := covariance.At(0, 0)
	varianceY := covariance.At(1, 1)
	varianceXY := covariance.At(0, 1)
	diff := (varianceX - varianceY) * 0.5

	sqr := math.Hypot(diff, varianceXY)
	eigenMajor := 0.5*(varianceX+varianceY) + sqr
	eigenMinor := 0.5*(varianceX+varianceY) - sqr

	k := calculateKChiSquared(CONFIDENCE, 2.0)
	angle := 0.5 * math.Atan2(2.0*varianceXY, diff)
	ellipse.MajorAxis = math.Sqrt(k * math.Max(0.0, eigenMajor))
	ellipse.MinorAxis = math.Sqrt(k * math.Max(0.0, eigenMinor))
	ellipse.Angle = angle / DEG_TO_RAD
	ellipse.Point = pcsToGcs(ellipse.Point, origin)
	return nil
}

/*
 * Description: Computes the uncertainty ellipsoid around the UE position estimate.
 * Parameter : ellipsoid
 * Type : UncertaintyEllipsoid3D
 * Parameter : covariance
 * Type : *mat.SymDense
 */
func (ellipsoid *UncertaintyEllipsoid3D) Figure(
	covariance *mat.SymDense,
	origin *mat.VecDense,
) error {
	varianceX := covariance.At(0, 0)
	varianceY := covariance.At(1, 1)
	varianceXY := covariance.At(0, 1)

	var eigsym mat.EigenSym
	if ok := eigsym.Factorize(covariance, false); !ok {
		return errors.New("cannot calculate 3 eigenvalues")
	}
	eigens := eigsym.Values(nil)

	k := calculateKChiSquared(CONFIDENCE, 3.0)
	angle := 0.5 * math.Atan2(2.0*varianceXY, varianceX-varianceY)
	ellipsoid.MajorAxis = math.Sqrt(k * math.Max(0.0, eigens[2]))
	ellipsoid.MidAxis = math.Sqrt(k * math.Max(0.0, eigens[1]))
	ellipsoid.MinorAxis = math.Sqrt(k * math.Max(0.0, eigens[0]))
	ellipsoid.Angle = angle / DEG_TO_RAD
	ellipsoid.Point = pcsToGcs(ellipsoid.Point, origin)
	return nil
}

/*
 * Description: Return the center point, major, minor axis and angle of orientation of ellipse.
 * Parameter : ellipse
 * Type : UncertaintyEllipse2D
 */
func (ellipse *UncertaintyEllipse2D) GetFigure() (
	point []float64,
	axis []float64,
	angle float64,
) {
	point = ellipse.Point.RawVector().Data
	axis = []float64{ellipse.MajorAxis, ellipse.MinorAxis}
	angle = ellipse.Angle
	return
}

/*
 * Description: Return the center point, major, mid, minor axis and angle of orientation of ellipsoid.
 * Parameter : ellipsoid
 * Type : UncertaintyEllipsoid3D
 */
func (ellipsoid *UncertaintyEllipsoid3D) GetFigure() (
	point []float64,
	axis []float64,
	angle float64,
) {
	point = ellipsoid.Point.RawVector().Data
	axis = []float64{ellipsoid.MajorAxis, ellipsoid.MidAxis, ellipsoid.MinorAxis}
	angle = ellipsoid.Angle
	return
}
