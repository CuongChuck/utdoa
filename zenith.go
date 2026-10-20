package utdoa

import (
	"fmt"
	"math"

	"gonum.org/v1/gonum/mat"
)

/*
 * Description: Set dimension for measurement
 * Parameter : zm
 * Type : *ZenithMeasurement
 * Parameter : dim
 * Type : int
 */
func (zm *ZenithMeasurement) SetDim(dim int) {
	zm.Dim = dim
}

/*
 * Description: Calculates zenith measurements from TRP measurements.
 * Parameter : zm
 * Type : *ZenithMeasurement
 * Parameter : measSet
 * Type : []TrpMeas
 */
func (zm *ZenithMeasurement) Calculate(measSet []TrpMeas) {
	numTrp := len(measSet)
	zm.Zeniths = make([]float64, numTrp)
	for idx, meas := range measSet {
		zm.Zeniths[idx] = meas.zenith * ANGLE_RESOLUTION * DEG_TO_RAD
		fmt.Printf("TRP %v - Zenith: %v\n", idx, zm.Zeniths[idx])
	}
}

/*
 * Description: Computes the zenith weighting matrix from TRP measurements.
 * Parameter : zm
 * Type : ZenithMeasurement
 * Parameter : measSet
 * Type : []TrpMeas
 */
func (zm ZenithMeasurement) CovarianceMatrix(measSet []TrpMeas) *mat.SymDense {
	numTrp := len(measSet)
	variances := make([]float64, numTrp)
	for i := range numTrp {
		uncertaintyZenith := measSet[i].zenithQuality.quality * measSet[i].zenithQuality.res
		uncertaintyZenithRadian := uncertaintyZenith * DEG_TO_RAD
		variances[i] = uncertaintyZenithRadian * uncertaintyZenithRadian
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
 * Description: Computes the Zenith Jacobian matrix and Residual vector for the current UE position
 * Parameter : zm
 * Type : ZenithMeasurement
 * Parameter : uePos
 * Type : *mat.VecDense
 * Parameter : measSet
 * Type : []TrpMeas
 * Parameter : distTrps
 * Type : []float64
 * Parameter : ctx
 * Type : *PositioningContext
 */
func (zm ZenithMeasurement) JacobianResidual(
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
		dZ := uePos.AtVec(2) - measSet[trpIdx].pos.AtVec(2)
		dX2 := dX * dX
		dY2 := dY * dY
		rho := dX2 + dY2
		if rho == 0 {
			continue
		}
		sqrtDeltaXY := math.Sqrt(rho)
		dist2 := distTrps[trpIdx] * distTrps[trpIdx]
		denom := dist2 * sqrtDeltaXY
		expected_zenith := math.Atan2(sqrtDeltaXY, dZ)
		j = append(j, []float64{((dZ * dX) / denom), ((dZ * dY) / denom), -(sqrtDeltaXY / dist2)})
		r = append(r, normalizeAngle(expected_zenith-zm.Zeniths[trpIdx]))
		filteredTrp = append(filteredTrp, trpIdx)
	}
	ctx.KeptTrp = filteredTrp
	return j, r
}

/*
 * Description: Filter Zenith Jacobian matrix and Residual vector from removed TRP measurements.
 * Parameter : zm
 * Type : ZenithMeasurement
 * Parameter : jacobian
 * Type : [][]float64
 * Parameter : residual
 * Type : []float64
 * Parameter : ctx
 * Type : *PositioningContext
 */
func (zm ZenithMeasurement) FilterJacobianResidual(
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
