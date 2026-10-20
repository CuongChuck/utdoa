# utdoa

`utdoa` is a Go library for estimating the position of a user equipment (UE) from
uplink radio measurements collected by a set of transmission/reception points
(TRPs).

The core solver is a **Gauss-Newton weighted least squares** estimator that fuses

- **TDOA** – time difference of arrival between a reference TRP and every other TRP
  (the classic UTDOA observable),
- **Azimuth** – the direction of arrival reported by a TRP (2D and 3D),
- **Zenith** – the elevation/zenith angle reported by a TRP (3D only),

and returns both the estimated UE position and the **uncertainty ellipse**
(2D) or **uncertainty ellipsoid** (3D) that surrounds it at a 90 % confidence
level.

The library speaks the 3GPP data encodings used by 5G NR positioning
(`RTOA`, TRP position codes, azimuth/zenith angle codes), so measurements taken
from a real gNB or from a test harness can be fed in directly.

---

## Table of contents

- [Features](#features)
- [Requirements](#requirements)
- [Installation](#installation)
- [Quick start (2D, TDOA only)](#quick-start-2d-tdoa-only)
- [3D with azimuth and zenith](#3d-with-azimuth-and-zenith)
- [API reference](#api-reference)
- [Data encoding reference](#data-encoding-reference)
- [Coordinate systems](#coordinate-systems)
- [Algorithm overview](#algorithm-overview)
- [Error handling](#error-handling)
- [Testing](#testing)
- [Project layout](#project-layout)
- [Known limitations](#known-limitations)

---

## Features

- TDOA positioning with automatic reference-TRP selection (the TRP with the
  strongest RSRP becomes the reference).
- Optional azimuth (angle of arrival) and zenith angle observations, mixed
  freely with TDOA in the same least-squares problem.
- 2D (latitude/longitude) and 3D (latitude/longitude/height) positioning.
- Per-measurement weighting through a quality/resolution descriptor; TDOA
  correlations with the reference TRP are modelled in the covariance matrix.
- Weighted Gauss-Newton with a backtracking line search, so a non-linear or
  badly conditioned problem degrades gracefully instead of diverging.
- Analytic uncertainty figure (ellipse / ellipsoid) derived from the inverse
  normal matrix and a chi-squared quantile.
- Small, dependency-light API built on top of
  [`gonum`](https://github.com/gonum/gonum) for the linear algebra.

## Requirements

| Requirement | Version |
| ----------- | ------- |
| Go          | >= 1.26 (see `go.mod`) |
| `gonum.org/v1/gonum` | v0.17.0 |

## Installation

The module is currently declared as `module utdoa` in `go.mod`, i.e. the import
path is simply `utdoa`. For a local checkout:

```bash
go mod init myapp
go mod edit -require=utdoa@v0.0.0
go mod edit -replace=utdoa=/absolute/path/to/utdoa
go mod tidy
```

If you publish the library under a real module path, change the `module` line
in `go.mod` and import it accordingly.

## Quick start (2D, TDOA only)

```go
package main

import (
	"fmt"
	"math"

	"utdoa"
)

const earthRadiusM = 6371000.0

type trpSite struct {
	name string
	lat  float64
	lon  float64
	rsrp uint8
}

// distanceM approximates the distance between two geographic points in meters.
func distanceM(aLat, aLon, bLat, bLon float64) float64 {
	lat0 := aLat * math.Pi / 180
	dLat := (bLat - aLat) * math.Pi / 180
	dLon := (bLon - aLon) * math.Pi / 180
	return math.Hypot(dLon*math.Cos(lat0), dLat) * earthRadiusM
}

func main() {
	// Ground truth position of the UE.
	ueLat, ueLon := 10.784950, 106.675865

	trps := []trpSite{
		{"TRP-0", 10.785962, 106.675179, 50},
		{"TRP-1", 10.784850, 106.673406, 90},
		{"TRP-2", 10.784781, 106.676000, 70},
		{"TRP-3", 10.788500, 106.679000, 60},
	}

	measSet := make([]utdoa.TrpMeas, 0, len(trps))
	for i, t := range trps {
		distance := distanceM(t.lat, t.lon, ueLat, ueLon)
		meas := utdoa.NewTrpMeas(
			i,                                        // TRP identifier
			utdoa.EncodeTrpPosition(t.lat, t.lon, 0, false), // 3GPP TRP position
			utdoa.EncodeRtoa(distance, 0),            // 3GPP RTOA value
			0,                                        // time granularity factor
			t.rsrp,                                   // RSRP used to pick the reference TRP
			false,                                    // position accuracy flag
		)
		meas.SetTdoaQuality(1, 0.1) // uncertainty = quality * res
		measSet = append(measSet, meas)
	}

	result, err := utdoa.UtdoaPositioning(measSet)
	if err != nil {
		fmt.Println("positioning failed:", err)
		return
	}

	point, axis, angle := result.GetFigure()
	fmt.Printf("UE position : lat=%.6f lon=%.6f\n", point[0], point[1])
	fmt.Printf("truth      : lat=%.6f lon=%.6f\n", ueLat, ueLon)
	fmt.Printf("uncertainty: major=%.2f m minor=%.2f m angle=%.2f deg\n",
		axis[0], axis[1], angle)
}
```

Output:

```text
UE position : lat=10.784950 lon=106.675865
truth      : lat=10.784950 lon=106.675865
uncertainty: major=0.17 m minor=0.14 m angle=38.00 deg
```

The same program is available as a runnable example:

```bash
go run ./examples/basic
```

## 3D with azimuth and zenith

Add `utdoa.With3D()`, `utdoa.WithAzimuth()` and `utdoa.WithZenith()` and provide
the angle observations of every TRP. Note that zenith measurements are only
meaningful in 3D — combining them with a 2D solve is rejected with an error.

```go
north, east := localOffsetM(trpLat, trpLon, ueLat-trpLat, ueLon-trpLon)
azimuthDeg := math.Atan2(east, north) * 180 / math.Pi        // from north, clockwise
zenithDeg := math.Atan2(math.Hypot(north, east), ueHeight-trpHeight) * 180 / math.Pi

meas.SetAzimuth(utdoa.EncodeAzimuth(azimuthDeg))
meas.SetZenith(utdoa.EncodeZenith(zenithDeg))

result, err := utdoa.UtdoaPositioning(
	measSet,
	utdoa.With3D(),
	utdoa.WithAzimuth(),
	utdoa.WithZenith(),
)
```

A complete version lives in `examples/angles`:

```bash
go run ./examples/angles
```

## API reference

### Entry point

```go
func UtdoaPositioning(measSet []TrpMeas, options ...UtdoaOption) (UncertaintyShape, error)
```

Runs the whole pipeline: coordinate conversion, measurement preparation,
Gauss-Newton solving, uncertainty computation. It returns an error when the
configuration or the data is invalid.

### Options

| Option | Description |
| ------ | ----------- |
| `WithDimension(dim int)` | Set the position dimension explicitly. Only `2` and `3` are valid. |
| `With3D()` | Shortcut for `WithDimension(3)`. |
| `WithAzimuth()` | Add an azimuth (AoA) measurement block to the problem. |
| `WithZenith()` | Add a zenith angle measurement block (requires 3D). |

The default configuration is 2D with a single TDOA block. Options are applied in
order, so the last one wins.

### Input data

```go
type TrpMeas struct { /* unexported fields */ }
```

A `TrpMeas` holds everything the solver knows about one TRP. Because the fields
are internal, construct it with the helpers below.

| Function | Description |
| -------- | ----------- |
| `NewTrpMeas(id int, pos3gpp *mat.VecDense, rtoa3gpp, factor float64, rsrp uint8, highAccuracy bool) TrpMeas` | Creates a TRP measurement from 3GPP encoded values. `id` is the TRP identifier, `pos3gpp` the encoded TRP position, `rtoa3gpp` the encoded round-trip time, `factor` the time granularity, `rsrp` the received power used to select the reference TRP and `highAccuracy` selects the high-accuracy position encoding. |
| `(*TrpMeas).SetAzimuth(azimuth float64)` | Sets the azimuth code (3GPP units). |
| `(*TrpMeas).SetZenith(zenith float64)` | Sets the zenith code (3GPP units). |
| `(*TrpMeas).SetTdoaQuality(quality, res float64)` | Sets the TDOA uncertainty. |
| `(*TrpMeas).SetAzimuthQuality(quality, res float64)` | Sets the azimuth uncertainty. |
| `(*TrpMeas).SetZenithQuality(quality, res float64)` | Sets the zenith uncertainty. |
| `(*TrpMeas).SetHighAccuracy(highAccuracy bool)` | Switches the position encoding. |
| `(*TrpMeas).ID() int` | Returns the TRP identifier. |

Encoding helpers:

| Function | Description |
| -------- | ----------- |
| `EncodeTrpPosition(latDeg, lonDeg, heightM float64, highAccuracy bool) *mat.VecDense` | Encodes a geographic position into the 3GPP TRP position vector. |
| `EncodeRtoa(distanceM float64, factor float64) float64` | Encodes a propagation distance into a 3GPP RTOA value. |
| `EncodeAzimuth(azimuthDeg float64) float64` | Encodes an azimuth angle in degrees into the 3GPP azimuth code. |
| `EncodeZenith(zenithDeg float64) float64` | Encodes a zenith angle in degrees into the 3GPP zenith code. |
| `NewMeasQuality(quality, res float64) MeasQuality` | Creates a quality descriptor. |

Measurement quality is expressed as `uncertainty = quality × res`:

- TDOA: the uncertainty is a time in seconds (the variance is its square).
- Azimuth / zenith: the uncertainty is an angle in degrees; it is converted to
  radians internally.

If any TRP of a measurement block reports a zero uncertainty, that block is
solved **unweighted** (plain least squares).

### Result

```go
type UncertaintyShape interface {
	Figure(covariance *mat.SymDense, originGcs *mat.VecDense) error
	GetFigure() (point []float64, axis []float64, angle float64)
}
```

| Type | Fields | Meaning |
| ---- | ------ | ------- |
| `UncertaintyEllipse2D` | `Point`, `MajorAxis`, `MinorAxis`, `Angle` | 2D result. `Point` is `[lat, lon]` in degrees, the axes are metres, `Angle` is the orientation of the major axis in degrees. |
| `UncertaintyEllipsoid3D` | `Point`, `MajorAxis`, `MidAxis`, `MinorAxis`, `Angle` | 3D result. `Point` is `[lat, lon, height]` (degrees, degrees, metres). |

`GetFigure()` returns the same information as a tuple:

- `point` – the estimated UE position in GCS (degrees, degrees, and metres in 3D),
- `axis` – `[major, minor]` in 2D or `[major, mid, minor]` in 3D, in metres,
- `angle` – orientation of the major axis in degrees.

Type-assert the result if you need the concrete struct:

```go
if ellipse, ok := result.(*utdoa.UncertaintyEllipse2D); ok {
	fmt.Println(ellipse.MajorAxis, ellipse.MinorAxis, ellipse.Angle)
}
```

### Custom measurement types

`MeasurementBlock` is the extension point for observables that are not built
in:

```go
type MeasurementBlock interface {
	Calculate(measSet []TrpMeas)
	JacobianResidual(uePos *mat.VecDense, measSet []TrpMeas, distTrps []float64,
		ctx *PositioningContext) ([][]float64, []float64)
	FilterJacobianResidual(jacobian [][]float64, residual []float64,
		ctx *PositioningContext) ([][]float64, []float64)
	CovarianceMatrix(measSet []TrpMeas) *mat.SymDense
	SetDim(dim int)
}
```

`Calculate` derives the raw observations, `JacobianResidual` returns the
Jacobian and the residual vector at the current UE position,
`FilterJacobianResidual` drops the rows of TRPs that were removed during the
iteration, and `CovarianceMatrix` returns the measurement covariance used for
weighting (return `nil` for an unweighted block).

## Data encoding reference

### TRP position

`EncodeTrpPosition` produces a `[x, y, z]` vector that is decoded as:

| Flag | Latitude | Longitude | Height |
| ---- | -------- | --------- | ------ |
| direct (`highAccuracy = false`) | `x / 2^23 × 90` | `y / 2^23 × 180` | `z` (metres) |
| high accuracy (`highAccuracy = true`) | `x / 2^31 × 90` | `y / 2^31 × 180` | `z / 2^7` (metres) |

In 2D the height component is ignored.

### RTOA

The round-trip time value is decoded to seconds with

```text
rtoa_seconds = BASE_TIME_UNIT × (LOWEST_RTOA + rtoa3gpp × 2^factor)
BASE_TIME_UNIT = 1 / (480000 × 4096)      // 1/15.36 MHz sub-frame granularity
LOWEST_RTOA   = -985024
```

`factor` is the time granularity exponent (0 → 1 unit, 1 → 2 units, …).
`EncodeRtoa(distanceM, factor)` is the inverse of that formula and is the
recommended way to build the value from a known distance.

### Angles

| Angle | 3GPP code | Decoded value |
| ----- | ---------- | ------------- |
| Azimuth | `code` | `code × 0.1 − 180` degrees, measured clockwise from north |
| Zenith  | `code` | `code × 0.1` degrees above the horizon |

### Reference TRP

For the TDOA block the TRP with the **highest RSRP** is used as the reference.
Every other TRP produces one TDOA observation relative to it, so a common
clock offset between the TRPs cancels out.

## Coordinate systems

- **GCS** – geographic coordinates: latitude and longitude in degrees, height
  in metres.
- **PCS** – a local Cartesian frame in metres where the origin is the first
  element of `measSet` (TRP 0) and the axes point north / east (and up in 3D).

Internally the solver works in PCS; the reported position and the uncertainty
figure are converted back to GCS relative to TRP 0.

## Algorithm overview

1. **Preparation** – TRP positions are decoded from 3GPP into GCS and then into
   PCS around TRP 0.
2. **Initial guess** – the UE position starts at the centroid of all TRP
   positions.
3. **Measurement blocks** – every block converts its raw data into observations
   (`TdoaMeasurement`, `AzimuthMeasurement`, `ZenithMeasurement`) and builds a
   covariance matrix for weighting.
4. **Gauss-Newton loop** (max. 50 iterations):
   - build the normal equations `JᵀWJ` and `JᵀWr` from all blocks, dropping
     TRPs that are closer than 1 cm to the current estimate,
   - solve the system with a Cholesky factorisation,
   - apply a backtracking line search (step halved from `1.0` down to `1e-6`)
     that only accepts a step when the cost does not increase,
   - stop when the step norm drops below `1e-3` or no descent step exists.
5. **Uncertainty** – the covariance matrix is the inverse of the final normal
   matrix; the ellipse/ellipsoid is obtained from its eigenvalues scaled by the
   chi-squared quantile for a 90 % confidence level (4.605 for 2 dof, 6.251 for
   3 dof).

## Error handling

`UtdoaPositioning` returns an error for the following situations:

| Error | Cause |
| ----- | ----- |
| `position location dimension must be 2 or 3` | `WithDimension` was called with an unsupported value. |
| `zenith measurement requires 3D positioning` | `WithZenith()` combined with a 2D solve. |
| `not enough measurements` | Fewer than two TRPs were supplied. |
| `normal matrix is rank deficient` | The geometry is degenerate (e.g. all TRPs collinear). |
| `measurement covariance is rank deficient` | The measurement covariance cannot be inverted. |

## Testing

```bash
go test ./...
```

The suite contains unit tests for every helper plus a CSV-driven regression
test that replays recorded scenarios:

| File | Purpose |
| ---- | ------- |
| `testcase.csv` | Template with the expected column layout. |
| `all_testcase.csv` | Full scenario set (`test_id`, `dim`, `ue_pos`, `trp_pos`, `tdoa_quality`, `azimuth_quality`, `zenith_quality`, `factor`, `error`). |
| `failed_testcase.csv` | Scenarios that did not meet the accuracy target; used by `TestGaussNewton_PositionInputs`. |

The `examples/` directory contains two runnable programs (`basic`, `angles`)
that double as end-to-end checks.

## Project layout

| File | Content |
| ---- | ------- |
| `api.go` | Exported constructors, setters and encoding helpers. |
| `positioning.go` | Public entry point `UtdoaPositioning` and the functional options. |
| `header.go` | Shared types and constants. |
| `init.go` | TRP coordinate conversion and the initial UE position. |
| `calculate.go` | Gauss-Newton solver and the uncertainty ellipse/ellipsoid. |
| `helper.go` | Coordinate conversions, chi-squared quantile, normal equation assembly. |
| `tdoa.go` | TDOA measurement block. |
| `azimuth.go` | Azimuth measurement block. |
| `zenith.go` | Zenith measurement block. |
| `testutil.go` | Test helpers and synthetic measurement generators. |

## Known limitations

- The solver writes progress information (TRP coordinates, Jacobian, residual,
  convergence trace) to `stdout` with `fmt.Printf`. Redirect or filter the
  output if that is unwanted in production.
- The TDOA covariance assumes that all TDOAs are perfectly correlated with the
  reference TRP and independent of each other.
- Zenith measurements are ignored in 2D solves (and rejected by the API).
- The module path in `go.mod` is the placeholder `utdoa`; change it before
  publishing the library.
- A few of the repository's own tests are currently failing (the CSV-driven
  regression test and the angle-unit assertions in `calculate_test.go`); they
  are pre-existing and unrelated to using the library.
