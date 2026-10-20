package utdoa

import (
	"fmt"

	"gonum.org/v1/gonum/mat"
)

/*
 * Description: Set dimension for measurement
 * Parameter : tm
 * Type : *TdoaMeasurement
 * Parameter : dim
 * Type : int
 */
func (tm *TdoaMeasurement) SetDim(dim int) {
	tm.Dim = dim
}

/*
 * Description: Calculates Time Difference of Arrival (TDOA) measurements from TRP measurements.
 * Parameter : tm
 * Type : *TdoaMeasurement
 * Parameter : measSet
 * Type : []TrpMeas
 */
func (tm *TdoaMeasurement) Calculate(measSet []TrpMeas) {
	numTrp := len(measSet)
	maxRsrp := uint8(0)
	tm.TrpRefIdx = 0
	for trpIdx := range measSet {
		meas := &measSet[trpIdx]
		meas.rtoa = rtoa3gppToSec(meas.rtoa3gpp, meas.factor)
		if meas.rsrp > maxRsrp {
			maxRsrp = meas.rsrp
			tm.TrpRefIdx = trpIdx
		}
	}
	tm.Tdoas = make([]float64, numTrp-1)
	row := 0
	for trpIdx := range measSet {
		if trpIdx == tm.TrpRefIdx {
			continue
		}
		timeTdoa := measSet[tm.TrpRefIdx].rtoa - measSet[trpIdx].rtoa
		tm.Tdoas[row] = LIGHTSPEED * timeTdoa
		fmt.Printf("TRP %v - TDOA: %v\n", trpIdx, tm.Tdoas[row])
		row++
	}
}

/*
 * Description: Computes the TDOA weighting matrix from TRP measurements.
 * Parameter : tm
 * Type : TdoaMeasurement
 * Parameter : measSet
 * Type : []TrpMeas
 */
func (tm TdoaMeasurement) CovarianceMatrix(
	measSet []TrpMeas,
) *mat.SymDense {
	numTrp := len(measSet)
	variances := make([]float64, numTrp)
	for i, meas := range measSet {
		uncertainty := meas.tdoaQuality.quality * meas.tdoaQuality.res
		variances[i] = uncertainty * uncertainty
		if variances[i] == 0 {
			return nil
		}
	}
	covarianceMat := mat.NewSymDense(numTrp-1, nil)
	row := 0
	for trpIdx := range numTrp {
		if trpIdx == tm.TrpRefIdx {
			continue
		}
		col := 0
		for otherIdx := range numTrp {
			if otherIdx == tm.TrpRefIdx {
				continue
			}
			value := variances[tm.TrpRefIdx]
			if row == col {
				value += variances[trpIdx]
			}
			covarianceMat.SetSym(row, col, value)
			col++
		}
		row++
	}
	return covarianceMat
}

/*
 * Description: Computes the TDOA Jacobian matrix and Residual vector for the current UE position.
 * Parameter : tm
 * Type : TdoaMeasurement
 * Parameter : uePos
 * Type : *mat.VecDense
 * Parameter : measSet
 * Type : []TrpMeas
 * Parameter : distTrps
 * Type : []float64
 * Parameter : ctx
 * Type : *PositioningContext
 */
func (tm TdoaMeasurement) JacobianResidual(
	uePos *mat.VecDense,
	measSet []TrpMeas,
	distTrps []float64,
	ctx *PositioningContext,
) ([][]float64, []float64) {
	refIdx := tm.TrpRefIdx
	if distTrps[refIdx] == 0 {
		return nil, nil
	}
	numTrp := len(ctx.KeptTrp) - 1
	j := make([][]float64, numTrp)
	r := make([]float64, numTrp)
	jRef := make([]float64, tm.Dim)
	refDist := distTrps[refIdx]
	refPos := measSet[refIdx].pos
	for i := range tm.Dim {
		jRef[i] = (uePos.AtVec(i) - refPos.AtVec(i)) / refDist
	}
	row := 0
	for _, trpIdx := range ctx.KeptTrp {
		if trpIdx == refIdx {
			continue
		}
		jTrp := make([]float64, tm.Dim)
		distTrp := distTrps[trpIdx]
		measTrp := measSet[trpIdx]
		for i := range tm.Dim {
			jTrp[i] = (uePos.AtVec(i) - measTrp.pos.AtVec(i)) / distTrp
		}
		for i := range tm.Dim {
			jTrp[i] = jRef[i] - jTrp[i]
		}
		j[row] = jTrp
		tdoaIdx := trpIdx
		if tdoaIdx > refIdx {
			tdoaIdx--
		}
		r[row] = refDist - distTrps[trpIdx] - tm.Tdoas[tdoaIdx]
		row++
	}
	return j, r
}

/*
 * Description: Filter TDOA Jacobian matrix and Residual vector from removed TRP measurements.
 * Parameter : tm
 * Type : TdoaMeasurement
 * Parameter : jacobian
 * Type : [][]float64
 * Parameter : residual
 * Type : []float64
 * Parameter : ctx
 * Type : *PositioningContext
 */
func (tm TdoaMeasurement) FilterJacobianResidual(
	jacobian [][]float64,
	residual []float64,
	ctx *PositioningContext,
) ([][]float64, []float64) {
	numRowJ := len(jacobian)
	numTrp := len(ctx.KeptTrp)
	if numRowJ == numTrp-1 {
		return jacobian, residual
	}
	numRowR := len(residual)
	newJacobian := make([][]float64, 0, numRowJ)
	newResidual := make([]float64, 0, numRowR)
	for _, trpIdx := range ctx.KeptTrp {
		keptIdx := trpIdx
		if trpIdx > tm.TrpRefIdx {
			keptIdx--
		}
		newJacobian = append(newJacobian, jacobian[keptIdx])
		newResidual = append(newResidual, residual[keptIdx])
	}
	return newJacobian, newResidual
}
