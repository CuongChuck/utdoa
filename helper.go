package utdoa

import (
	"errors"
	"fmt"
	"math"

	"gonum.org/v1/gonum/mat"
	"gonum.org/v1/gonum/stat/distuv"
)

/*
 * Description: Calculates the squared Euclidean distance between two 3D positions.
 * Parameter : uePos
 * Type : *mat.VecDense
 * Parameter : trpPos
 * Type : *mat.VecDense
 */
func squaredDist(uePos, trpPos *mat.VecDense) float64 {
	var sum float64
	for i := range uePos.Len() {
		d := uePos.AtVec(i) - trpPos.AtVec(i)
		sum += d * d
	}
	return sum
}

/*
 * Description: Calculates the chi-squared quantile for a given probability and degrees of freedom.
 * Parameter : p
 * Type : float64
 * Parameter : d
 * Type : float64
 */
func calculateKChiSquared(p float64, d float64) float64 {
	chiDist := distuv.ChiSquared{
		K:   d,
		Src: nil,
	}
	return chiDist.Quantile(p)
}

/*
 * Description: Normalize angle residual
 * Parameter : angle
 * Type : float64
 */
func normalizeAngle(angle float64) float64 {
	return math.Atan2(math.Sin(angle), math.Cos(angle))
}

/*
 * Description: Remove rows and columns to get subMatrix
 * Parameter : matrix
 * Type : *mat.Dense
 * Parameter : keep
 * Type : []int
 */
func subMatrix(matrix *mat.SymDense, keep []int) *mat.SymDense {
	size := len(keep)
	result := mat.NewSymDense(size, nil)
	for row := range size {
		for col := 0; col <= row; col++ {
			result.SetSym(row, col, matrix.At(keep[row], keep[col]))
		}
	}
	return result
}

/*
 * Description: Convert coordinates from PCS (meters) format to GCS.
 * Parameter : pos
 * Type : *mat.VecDense
 * Parameter : origin
 * Type : *mat.VecDense
 */
func pcsToGcs(pos *mat.VecDense, origin *mat.VecDense) *mat.VecDense {
	if pos == nil || origin == nil {
		return pos
	}
	dim, _ := pos.Dims()
	result := make([]float64, dim)
	lat0Rad := origin.AtVec(0) * DEG_TO_RAD
	result[0] = origin.AtVec(0) + (pos.AtVec(0)/EARTH_RADIUS)/DEG_TO_RAD
	result[1] = origin.AtVec(1) + (pos.AtVec(1)/(EARTH_RADIUS*math.Cos(lat0Rad)))/DEG_TO_RAD
	if dim == 3 {
		result[2] = pos.AtVec(2)
	}
	return mat.NewVecDense(dim, result)
}

/*
 * Description: Convert coordinates from GCS format to PCS (meters).
 * Parameter : pos
 * Type : *mat.VecDense
 * Parameter : origin
 * Type : *mat.VecDense
 */
func gcsToPcs(pos *mat.VecDense, origin *mat.VecDense) *mat.VecDense {
	dim, _ := pos.Dims()
	result := make([]float64, dim)
	lat0Rad := origin.AtVec(0) * DEG_TO_RAD
	dLong := (pos.AtVec(1) - origin.AtVec(1)) * DEG_TO_RAD
	dLat := (pos.AtVec(0) - origin.AtVec(0)) * DEG_TO_RAD
	result[1] = dLong * EARTH_RADIUS * math.Cos(lat0Rad)
	result[0] = dLat * EARTH_RADIUS
	if dim == 3 {
		result[2] = pos.AtVec(2)
	}
	return mat.NewVecDense(dim, result)
}

/*
 * Description: Convert RTOA from 3GPP format to seconds
 * Parameter : rtoa3gpp
 * Type : float64
 * Parameter : factor
 * Type : float64
 */
func rtoa3gppToSec(rtoa3gpp, factor float64) float64 {
	return BASE_TIME_UNIT * (LOWEST_RTOA + rtoa3gpp*math.Pow(BASE_TIME, factor))
}

/*
 * Description: Build normal equation components in Gauss-Newton.
 * Parameter : uePos
 * Type : *mat.VecDense
 * Parameter : measSet
 * Type : []TrpMeas
 * Parameter : distTrps
 * Type : []float64
 * Parameter : keep
 * Type : []int
 * Parameter : measurements
 * Type : []MeasurementBlock
 * Parameter : covarianceMats
 * Type : []*mat.SymDense
 * Parameter : dim
 * Type : int
 */
func buildNormalEquation(
	uePos *mat.VecDense,
	measSet []TrpMeas,
	distTrps []float64,
	keep []int,
	measurements []MeasurementBlock,
	covarianceMats []*mat.SymDense,
	dim int,
) (*mat.SymDense, *mat.VecDense, float64, error) {
	isWeight := true
	numMeas := len(measurements)
	if numMeas != len(covarianceMats) {
		return nil, nil, 0, fmt.Errorf(
			"number of measurements (%d) does not match number of covariance matrices (%d)",
			numMeas,
			len(covarianceMats),
		)
	}
	ctx := &PositioningContext{KeptTrp: keep}
	jacobianSets := make([][][]float64, numMeas)
	residualSets := make([][]float64, numMeas)
	covariances := make([]*mat.SymDense, numMeas)
	coefficient := mat.NewSymDense(dim, nil)
	constant := mat.NewVecDense(dim, nil)
	var (
		jacobians     [][]float64
		residuals     []float64
		chol          mat.Cholesky
		weight        mat.SymDense
		jW            mat.Dense
		residualBlock mat.VecDense
		coeffBlock    mat.Dense
		cost          float64
	)
	totalRows := 0

	for idx, measurement := range measurements {
		jacobian, residual := measurement.JacobianResidual(
			uePos,
			measSet,
			distTrps,
			ctx,
		)
		rowsJ := len(jacobian)
		for row, values := range jacobian {
			if len(values) != dim {
				return nil, nil, 0, fmt.Errorf(
					"invalid Jacobian dimension at row %d: got %d columns, want %d",
					row,
					len(values),
					dim,
				)
			}
		}
		if rowsJ == 0 {
			return nil, nil, 0, fmt.Errorf("invalid Jacobian: empty Jacobian")
		}
		colsJ := len(jacobian[0])
		rowsR := len(residual)
		if colsJ != dim {
			return nil, nil, 0, fmt.Errorf(
				"invalid Jacobian dimension: %dx%d",
				rowsJ,
				colsJ,
			)
		}
		if rowsR != rowsJ {
			return nil, nil, 0, fmt.Errorf(
				"invalid residual dimension: got %d, want %d",
				rowsR,
				rowsJ,
			)
		}
		jacobianSets[idx] = jacobian
		residualSets[idx] = residual
	}

	for idx, measurement := range measurements {
		jacobian, residual := measurement.FilterJacobianResidual(jacobianSets[idx], residualSets[idx], ctx)
		rowsJ := len(jacobian)
		if rowsJ == 0 {
			return nil, nil, 0, fmt.Errorf("invalid Jacobian: empty Jacobian")
		}
		for row, values := range jacobian {
			if len(values) != dim {
				return nil, nil, 0, fmt.Errorf(
					"invalid Jacobian dimension at row %d: got %d columns, want %d",
					row,
					len(values),
					dim,
				)
			}
		}
		colsJ := len(jacobian[0])
		rowsR := len(residual)
		if colsJ != dim {
			return nil, nil, 0, fmt.Errorf(
				"invalid Jacobian dimension: %dx%d",
				rowsJ,
				colsJ,
			)
		}
		if rowsR != rowsJ {
			return nil, nil, 0, fmt.Errorf(
				"invalid residual dimension: got %d, want %d",
				rowsR,
				rowsJ,
			)
		}
		jacobians = append(jacobians, jacobian...)
		residuals = append(residuals, residual...)
		covariance := covarianceMats[idx]
		if covariance != nil {
			if len(ctx.KeptTrp) != len(keep) {
				covariances[idx] = subMatrix(covariance, ctx.KeptTrp)
			} else {
				covariances[idx] = covariance
			}
		} else {
			isWeight = false
		}
		totalRows += rowsJ
	}
	j := mat.NewDense(totalRows, dim, nil)
	r := mat.NewVecDense(totalRows, nil)
	for row := range totalRows {
		j.SetRow(row, jacobians[row])
		r.SetVec(row, residuals[row])
	}
	fmt.Printf("Jacobian:\n%v\n", mat.Formatted(j))
	fmt.Printf("Residual:\n%v\n", mat.Formatted(r))
	if !isWeight {
		coeffBlock.Mul(j.T(), j)
		residualBlock.MulVec(j.T(), r)
		cost = mat.Dot(r, r)
	} else {
		c := mat.NewSymDense(totalRows, nil)
		offset := 0
		for _, covariance := range covariances {
			d := covariance.SymmetricDim()
			for i := range d {
				for j := i; j < d; j++ {
					c.SetSym(offset+i, offset+j, covariance.At(i, j))
				}
			}
			offset += d
		}
		if ok := chol.Factorize(c); !ok {
			return nil, nil, 0, errors.New(
				"measurement covariance is rank deficient",
			)
		}
		if err := chol.InverseTo(&weight); err != nil {
			return nil, nil, 0, err
		}
		jW.Mul(j.T(), &weight)
		coeffBlock.Mul(&jW, j)
		var wr mat.VecDense
		wr.MulVec(&weight, r)
		residualBlock.MulVec(&jW, r)
		cost = mat.Dot(r, &wr)
	}
	for row := range dim {
		for col := 0; col <= row; col++ {
			coefficient.SetSym(row, col, coeffBlock.At(row, col))
		}
	}
	constant.AddVec(constant, &residualBlock)
	return coefficient, constant, cost, nil
}
