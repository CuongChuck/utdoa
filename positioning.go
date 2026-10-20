package utdoa

import (
	"errors"

	"gonum.org/v1/gonum/mat"
)

/*
 * Description: Caller assign dimension of position vector
 * Parameter : dim
 * Type : int
 */
func WithDimension(dim int) UtdoaOption {
	return func(cfg *UtdoaConfig) {
		cfg.Dimension = dim
	}
}

/*
 * Description: Caller assign 3D dimension of position vector
 */
func With3D() UtdoaOption {
	return WithDimension(3)
}

/*
 * Description: Caller add azimuth measurements to the algorithm
 */
func WithAzimuth() UtdoaOption {
	return func(cfg *UtdoaConfig) {
		cfg.Measurements = append(
			cfg.Measurements,
			&AzimuthMeasurement{},
		)
	}
}

/*
 * Description: Caller add zenith measurements to the algorithm
 */
func WithZenith() UtdoaOption {
	return func(cfg *UtdoaConfig) {
		cfg.Measurements = append(
			cfg.Measurements,
			&ZenithMeasurement{},
		)
	}
}

/*
 * Description: Default configuration for UTDOA positioning
 */
func defaultUtdoaConfig() UtdoaConfig {
	return UtdoaConfig{
		Dimension:    2,
		Measurements: []MeasurementBlock{&TdoaMeasurement{}},
	}
}

/*
 * Description: Main function handling the overall UTDOA positioning process.
 * Parameter : measSet
 * Type : []trpMeas
 * Parameter : options
 * Type : ...UtdoaOption
 */
func UtdoaPositioning(
	measSet []TrpMeas,
	options ...UtdoaOption,
) (UncertaintyShape, error) {
	cfg := defaultUtdoaConfig()
	for _, option := range options {
		option(&cfg)
	}
	if cfg.Dimension != 2 && cfg.Dimension != 3 {
		return nil, errors.New("position location dimension must be 2 or 3")
	}
	if cfg.Dimension == 2 {
		for _, m := range cfg.Measurements {
			if _, ok := m.(*ZenithMeasurement); ok {
				return nil, errors.New(
					"zenith measurement requires 3D positioning",
				)
			}
		}
	}
	if len(measSet) < 2 {
		return nil, errors.New("not enough measurements")
	}

	numMeasurement := len(cfg.Measurements)
	calculateTrpPos(measSet, cfg.Dimension)
	initialUePos, err := initUePos(measSet, cfg.Dimension)
	if err != nil {
		return nil, err
	}
	covarianceMats := make([]*mat.SymDense, numMeasurement)
	for idx, measurement := range cfg.Measurements {
		measurement.SetDim(cfg.Dimension)
		measurement.Calculate(measSet)
		covarianceMats[idx] = measurement.CovarianceMatrix(measSet)
	}
	gNResult, err := gaussNewton(
		initialUePos,
		measSet,
		cfg.Measurements,
		covarianceMats,
		cfg.Dimension,
	)
	if err != nil {
		return nil, err
	}
	var result UncertaintyShape
	switch cfg.Dimension {
	case 2:
		result = &UncertaintyEllipse2D{Point: gNResult.UePos}
	case 3:
		result = &UncertaintyEllipsoid3D{
			UncertaintyEllipse2D: UncertaintyEllipse2D{Point: gNResult.UePos},
		}
	}
	if err := result.Figure(gNResult.CovarMat, measSet[0].posGcs); err != nil {
		return nil, err
	}
	return result, nil
}
