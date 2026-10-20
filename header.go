package utdoa

import (
	"math"

	"gonum.org/v1/gonum/mat"
)

type MeasQuality struct {
	quality float64
	res     float64
}

type TrpMeas struct {
	id             int
	rtoa3gpp       float64
	rtoa           float64
	factor         float64
	rsrp           uint8
	pos            *mat.VecDense
	posGcs         *mat.VecDense
	pos3gpp        *mat.VecDense
	isHighAcc      bool
	azimuth        float64
	zenith         float64
	tdoaQuality    MeasQuality
	azimuthQuality MeasQuality
	zenithQuality  MeasQuality
}

type UncertaintyEllipse2D struct {
	Point     *mat.VecDense
	MajorAxis float64
	MinorAxis float64
	Angle     float64
}

type UncertaintyEllipsoid3D struct {
	UncertaintyEllipse2D
	MidAxis float64
}

type UncertaintyShape interface {
	Figure(covariance *mat.SymDense, originGcs *mat.VecDense) error
	GetFigure() (point []float64, axis []float64, angle float64)
}

type MeasurementBlock interface {
	Calculate(measSet []TrpMeas)

	JacobianResidual(
		uePos *mat.VecDense,
		measSet []TrpMeas,
		distTrps []float64,
		ctx *PositioningContext,
	) ([][]float64, []float64)

	FilterJacobianResidual(
		jacobian [][]float64,
		residual []float64,
		ctx *PositioningContext,
	) ([][]float64, []float64)

	CovarianceMatrix(measSet []TrpMeas) *mat.SymDense

	SetDim(dim int)
}

type PositioningContext struct {
	KeptTrp []int
}

type AzimuthMeasurement struct {
	Azimuths []float64
	Dim      int
}

type ZenithMeasurement struct {
	Zeniths []float64
	Dim     int
}

type TdoaMeasurement struct {
	Tdoas     []float64
	TrpRefIdx int
	Dim       int
}

type UtdoaConfig struct {
	Dimension    int
	Measurements []MeasurementBlock
}

type GaussNewtonResult struct {
	UePos       *mat.VecDense
	CovarMat    *mat.SymDense
	IsConverged bool
}

type UtdoaOption func(*UtdoaConfig)

const (
	BASE_TIME_UNIT   float64 = 1.0 / (480000.0 * 4096.0)
	LOWEST_RTOA              = -985024.0
	BASE_TIME                = 2.0
	LIGHTSPEED               = 299792458.0
	CONFIDENCE               = 0.9
	ANGLE_RESOLUTION         = 0.1
	DEGREE_PER_PI            = 180.0
	LATITUDE_CONST           = 90.0
	LONGITUDE_CONST          = 180.0
	DIRECT_FACTOR            = 23.0
	HIGH_FACTOR              = 31.0
	EARTH_RADIUS             = 6371000.0
	DEG_TO_RAD               = math.Pi / 180.0

	MAX_LOOP      = 50
	MIN_DIST_DIFF = 1e-3
	MIN_DIST      = 1e-2
	MIN_DIST2     = MIN_DIST * MIN_DIST
)
