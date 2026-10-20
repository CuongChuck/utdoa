package utdoa

import (
	"fmt"
	"math"

	"gonum.org/v1/gonum/mat"
)

/*
 * Description: Initializes UE position as the centroid of all TRP positions.
 * Parameter : measSet
 * Type : []trpMeas
 * Parameter : dim
 * Type : int
 */
func initUePos(measSet []TrpMeas, dim int) (*mat.VecDense, error) {
	numTrp := len(measSet)
	pos := mat.NewVecDense(dim, nil)
	for _, meas := range measSet {
		pos.AddVec(pos, meas.pos)
	}
	pos.ScaleVec(1.0/float64(numTrp), pos)
	return pos, nil
}

/*
 * Description: Convert TRP coordinates into degrees from measurements.
 * Parameter : measSet
 * Type : []trpMeas
 * Parameter : dim
 * Type : int
 */
func calculateTrpPos(measSet []TrpMeas, dim int) {
	for idx := range measSet {
		measSet[idx].posGcs = mat.NewVecDense(dim, nil)
		if measSet[idx].isHighAcc {
			measSet[idx].posGcs.SetVec(1,
				measSet[idx].pos3gpp.AtVec(1)/math.Pow(2.0, HIGH_FACTOR)*LONGITUDE_CONST)
			measSet[idx].posGcs.SetVec(0,
				measSet[idx].pos3gpp.AtVec(0)/math.Pow(2.0, HIGH_FACTOR)*LATITUDE_CONST)
			if dim == 3 {
				measSet[idx].posGcs.SetVec(2, measSet[idx].pos3gpp.AtVec(2)*math.Pow(2.0, -7.0))
			}
		} else {
			measSet[idx].posGcs.SetVec(1,
				measSet[idx].pos3gpp.AtVec(1)/math.Pow(2.0, DIRECT_FACTOR)*LONGITUDE_CONST)
			measSet[idx].posGcs.SetVec(0,
				measSet[idx].pos3gpp.AtVec(0)/math.Pow(2.0, DIRECT_FACTOR)*LATITUDE_CONST)
			if dim == 3 {
				measSet[idx].posGcs.SetVec(2, measSet[idx].pos3gpp.AtVec(2))
			}
		}
	}
	measSet[0].pos = mat.NewVecDense(dim, nil)
	if dim == 3 {
		measSet[0].pos.SetVec(2, measSet[0].posGcs.AtVec(2))
	}
	fmt.Printf("TRP 0 Cartesian coordinates: %v\n", measSet[0].pos.RawVector().Data)
	for idx := 1; idx < len(measSet); idx++ {
		measSet[idx].pos = mat.NewVecDense(dim, nil)
		measSet[idx].pos.CopyVec(gcsToPcs(measSet[idx].posGcs, measSet[0].posGcs))
		fmt.Printf("TRP %v Cartesian coordinates: %v\n", idx, measSet[idx].pos.RawVector().Data)
	}
}
