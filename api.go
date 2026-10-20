package utdoa

import (
	"math"

	"gonum.org/v1/gonum/mat"
)

/*
 * Description: Create a TRP measurement from 3GPP encoded values.
 * Parameter : id
 * Type : int
 * Parameter : pos3gpp
 * Type : *mat.VecDense
 * Parameter : rtoa3gpp
 * Type : float64
 * Parameter : factor
 * Type : float64
 * Parameter : rsrp
 * Type : uint8
 * Parameter : highAccuracy
 * Type : bool
 */
func NewTrpMeas(
	id int,
	pos3gpp *mat.VecDense,
	rtoa3gpp float64,
	factor float64,
	rsrp uint8,
	highAccuracy bool,
) TrpMeas {
	return TrpMeas{
		id:        id,
		pos3gpp:   pos3gpp,
		rtoa3gpp:  rtoa3gpp,
		factor:    factor,
		rsrp:      rsrp,
		isHighAcc: highAccuracy,
	}
}

/*
 * Description: Set the azimuth measurement of a TRP in 3GPP encoded units.
 * Parameter : meas
 * Type : *TrpMeas
 * Parameter : azimuth
 * Type : float64
 */
func (meas *TrpMeas) SetAzimuth(azimuth float64) {
	meas.azimuth = azimuth
}

/*
 * Description: Set the zenith measurement of a TRP in 3GPP encoded units.
 * Parameter : meas
 * Type : *TrpMeas
 * Parameter : zenith
 * Type : float64
 */
func (meas *TrpMeas) SetZenith(zenith float64) {
	meas.zenith = zenith
}

/*
 * Description: Set the TDOA measurement quality of a TRP.
 * Parameter : meas
 * Type : *TrpMeas
 * Parameter : quality
 * Type : float64
 * Parameter : res
 * Type : float64
 */
func (meas *TrpMeas) SetTdoaQuality(quality float64, res float64) {
	meas.tdoaQuality = MeasQuality{quality: quality, res: res}
}

/*
 * Description: Set the azimuth measurement quality of a TRP.
 * Parameter : meas
 * Type : *TrpMeas
 * Parameter : quality
 * Type : float64
 * Parameter : res
 * Type : float64
 */
func (meas *TrpMeas) SetAzimuthQuality(quality float64, res float64) {
	meas.azimuthQuality = MeasQuality{quality: quality, res: res}
}

/*
 * Description: Set the zenith measurement quality of a TRP.
 * Parameter : meas
 * Type : *TrpMeas
 * Parameter : quality
 * Type : float64
 * Parameter : res
 * Type : float64
 */
func (meas *TrpMeas) SetZenithQuality(quality float64, res float64) {
	meas.zenithQuality = MeasQuality{quality: quality, res: res}
}

/*
 * Description: Set the position accuracy flag of a TRP.
 * Parameter : meas
 * Type : *TrpMeas
 * Parameter : highAccuracy
 * Type : bool
 */
func (meas *TrpMeas) SetHighAccuracy(highAccuracy bool) {
	meas.isHighAcc = highAccuracy
}

/*
 * Description: Return the identifier of a TRP measurement.
 * Type : *TrpMeas
 */
func (meas *TrpMeas) ID() int {
	return meas.id
}

/*
 * Description: Create a measurement quality descriptor.
 * Parameter : quality
 * Type : float64
 * Parameter : res
 * Type : float64
 */
func NewMeasQuality(quality float64, res float64) MeasQuality {
	return MeasQuality{quality: quality, res: res}
}

/*
 * Description: Encode a geographic position into the 3GPP TRP position format.
 * Parameter : latDeg
 * Type : float64
 * Parameter : lonDeg
 * Type : float64
 * Parameter : heightM
 * Type : float64
 * Parameter : highAccuracy
 * Type : bool
 */
func EncodeTrpPosition(
	latDeg float64,
	lonDeg float64,
	heightM float64,
	highAccuracy bool,
) *mat.VecDense {
	factor := DIRECT_FACTOR
	if highAccuracy {
		factor = HIGH_FACTOR
	}
	pos := mat.NewVecDense(3, []float64{
		latDeg / LATITUDE_CONST * math.Pow(2.0, factor),
		lonDeg / LONGITUDE_CONST * math.Pow(2.0, factor),
		heightM,
	})
	if highAccuracy {
		pos.SetVec(2, heightM*math.Pow(2.0, 7.0))
	}
	return pos
}

/*
 * Description: Encode a propagation distance into a 3GPP RTOA value.
 * Parameter : distanceM
 * Type : float64
 * Parameter : factor
 * Type : float64
 */
func EncodeRtoa(distanceM float64, factor float64) float64 {
	return (distanceM/LIGHTSPEED/BASE_TIME_UNIT - LOWEST_RTOA) /
		math.Pow(BASE_TIME, factor)
}

/*
 * Description: Encode an azimuth angle in degrees into the 3GPP azimuth code.
 * Parameter : azimuthDeg
 * Type : float64
 */
func EncodeAzimuth(azimuthDeg float64) float64 {
	return (azimuthDeg + DEGREE_PER_PI) / ANGLE_RESOLUTION
}

/*
 * Description: Encode a zenith angle in degrees into the 3GPP zenith code.
 * Parameter : zenithDeg
 * Type : float64
 */
func EncodeZenith(zenithDeg float64) float64 {
	return zenithDeg / ANGLE_RESOLUTION
}
