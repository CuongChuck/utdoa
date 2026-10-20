package utdoa

import (
	"fmt"
	"math"

	"gonum.org/v1/gonum/mat"
)

/*
 * Description: Set dimension for measurement
 * Parameter : am
 * Type : *AzimuthMeasurement
 * Parameter : dim
 * Type : int
 */
func (am *AzimuthMeasurement) SetDim(dim int) {
	am.Dim = dim
}

/*
 * Description: Calculates azimuth measurements from TRP measurements.
 * Parameter : am
 * Type : *AzimuthMeasurement
 * Parameter : measSet
 * Type : []TrpMeas
 */
func (am *AzimuthMeasurement) Calculate(measSet []TrpMeas) {
	numTrp := len(measSet)
	am.Azimuths = make([]float64, numTrp)
	for idx, meas := range measSet {
		degree := meas.azimuth*ANGLE_RESOLUTION - DEGREE_PER_PI
		am.Azimuths[idx] = degree * DEG_TO_RAD
		fmt.Printf("TRP %v - Azimuth: %v\n", idx, am.Azimuths[idx])
	}
}

/*
 * Description: Computes the azimuth weighting matrix from TRP measurements.
 * Parameter : am
 * Type : AzimuthMeasurement
 * Parameter : measSet
 * Type : []TrpMeas
 */
func (am AzimuthMeasurement) CovarianceMatrix(measSet []TrpMeas) *mat.SymDense {
	numTrp := len(measSet)
	variances := make([]float64, numTrp)
	for i := range numTrp {
		uncertaintyAzimuth := measSet[i].azimuthQuality.quality * measSet[i].azimuthQuality.res
		uncertaintyAzimuthRadian := uncertaintyAzimuth * DEG_TO_RAD
		variances[i] = uncertaintyAzimuthRadian * uncertaintyAzimuthRadian
		if variances[i] == 0 {
			return nil
		}
	}

	covarianceMat := mat.NewSymDense(numTrp, nil)
	for i := range numTrp {
		covarianceMat.SetSym(i, i, variances[i])
	}
	return covarianceMat
}

/*
 * Description: Computes the azimuth Jacobian matrix and Residual vector for the current UE position
 * Parameter : am
 * Type : AzimuthMeasurement
 * Parameter : uePos
 * Type : *mat.VecDense
 * Parameter : measSet
 * Type : []TrpMeas
 * Parameter : distTrps
 * Type : []float64
 * Parameter : ctx
 * Type : *PositioningContext
 */
func (am AzimuthMeasurement) JacobianResidual(
	uePos *mat.VecDense,
	measSet []TrpMeas,
	distTrps []float64,
	ctx *PositioningContext,
) ([][]float64, []float64) {
	numTrp := len(ctx.KeptTrp)
	filteredTrp := ctx.KeptTrp[:0]
	j := make([][]float64, 0, numTrp)
	r := make([]float64, 0, numTrp)
	for _, trpIdx := range ctx.KeptTrp {
		dX := uePos.AtVec(0) - measSet[trpIdx].pos.AtVec(0)
		dY := uePos.AtVec(1) - measSet[trpIdx].pos.AtVec(1)
		dX2 := dX * dX
		dY2 := dY * dY
		dXY := dX2 + dY2
		if dXY == 0 {
			continue
		}
		expected_azimuth := math.Atan2(dY, dX)
		jRow := []float64{-(dY / dXY), dX / dXY}
		if am.Dim == 3 {
			jRow = append(jRow, 0)
		}
		j = append(j, jRow)
		r = append(r, normalizeAngle(expected_azimuth-am.Azimuths[trpIdx]))
		filteredTrp = append(filteredTrp, trpIdx)
	}
	ctx.KeptTrp = filteredTrp
	return j, r
}

/*
 * Description: Filter Azimuth Jacobian matrix and Residual vector from removed TRP measurements.
 * Parameter : am
 * Type : AzimuthMeasurement
 * Parameter : jacobian
 * Type : [][]float64
 * Parameter : residual
 * Type : []float64
 * Parameter : ctx
 * Type : *PositioningContext
 */
func (am AzimuthMeasurement) FilterJacobianResidual(
	jacobian [][]float64,
	residual []float64,
	ctx *PositioningContext,
) ([][]float64, []float64) {
	numRowJ := len(jacobian)
	numTrp := len(ctx.KeptTrp)
	if numRowJ == numTrp {
		return jacobian, residual
	}
	numRowR := len(residual)
	newJacobian := make([][]float64, 0, numRowJ)
	newResidual := make([]float64, 0, numRowR)
	for _, trpIdx := range ctx.KeptTrp {
		newJacobian = append(newJacobian, jacobian[trpIdx])
		newResidual = append(newResidual, residual[trpIdx])
	}
	return newJacobian, newResidual
}
