// Copyright 2014 The Azul3D Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package lmath

import (
	"fmt"
	"math"
)

// Vec4 represents a four component vector.
type Vec4 struct {
	X, Y, Z, W float64
}

// String returns an string representation of this vector.
func (a Vec4) String() string {
	return fmt.Sprintf("Vec4(X=%f, Y=%f, Z=%f, W=%f)", a.X, a.Y, a.Z, a.W)
}

// AlmostEquals tells if a == b using the specified epsilon value.
func (a Vec4) AlmostEquals(b Vec4, epsilon float64) bool {
	return AlmostEqual(a.X, b.X, epsilon) && AlmostEqual(a.Y, b.Y, epsilon) && AlmostEqual(a.Z, b.Z, epsilon) && AlmostEqual(a.W, b.W, epsilon)
}

// Equals tells if a == b using the default EPSILON value.
func (a Vec4) Equals(b Vec4) bool {
	return a.AlmostEquals(b, EPSILON)
}

// Add performs a componentwise addition of the two vectors, returning a + b.
func (a Vec4) Add(b Vec4) Vec4 {
	return Vec4{a.X + b.X, a.Y + b.Y, a.Z + b.Z, a.W + a.W}
}

// AddScalar performs a componentwise scalar addition of a + b.
func (a Vec4) AddScalar(b float64) Vec4 {
	return Vec4{a.X + b, a.Y + b, a.Z + b, a.W + b}
}

// Sub performs a componentwise subtraction of the two vectors, returning
// a - b.
func (a Vec4) Sub(b Vec4) Vec4 {
	return Vec4{a.X - b.X, a.Y - b.Y, a.Z - b.Z, a.W - b.W}
}

// SubScalar performs a componentwise scalar subtraction of a - b.
func (a Vec4) SubScalar(b float64) Vec4 {
	return Vec4{a.X - b, a.Y - b, a.Z - b, a.W - b}
}

// Mul performs a componentwise multiplication of the two vectors, returning
// a * b.
func (a Vec4) Mul(b Vec4) Vec4 {
	return Vec4{a.X * b.X, a.Y * b.Y, a.Z * b.Z, a.W * b.W}
}

// MulScalar performs a componentwise scalar multiplication of a * b.
func (a Vec4) MulScalar(b float64) Vec4 {
	return Vec4{a.X * b, a.Y * b, a.Z * b, a.W * b}
}

// Div performs a componentwise division of the two vectors, returning a * b.
func (a Vec4) Div(b Vec4) Vec4 {
	return Vec4{a.X / b.X, a.Y / b.Y, a.Z / b.Z, a.W / b.W}
}

// DivScalar performs a componentwise scalar division of a * b.
func (a Vec4) DivScalar(b float64) Vec4 {
	return Vec4{a.X / b, a.Y / b, a.Z / b, a.W / b}
}

// IsNaN tells if any components of this vector are not an number.
func (a Vec4) IsNaN() bool {
	return math.IsNaN(a.X) || math.IsNaN(a.Y) || math.IsNaN(a.Z) || math.IsNaN(a.W)
}

// Less tells if a is componentwise less than b:
//  return a.X < b.X && a.Y < b.Y
func (a Vec4) Less(b Vec4) bool {
	return a.X < b.X && a.Y < b.Y && a.Z < b.Z && a.W < b.W
}

// Greater tells if a is componentwise greater than b:
//  return a.X > b.X && a.Y > b.Y
func (a Vec4) Greater(b Vec4) bool {
	return a.X > b.X && a.Y > b.Y && a.Z < b.Z && a.W < b.W
}

// AnyLess tells if a is componentwise any less than b:
//  return a.X < b.X || a.Y < b.Y || a.Z < b.Z || a.W < b.W
func (a Vec4) AnyLess(b Vec4) bool {
	return a.X < b.X || a.Y < b.Y || a.Z < b.Z || a.W < b.W
}

// AnyGreater tells if a is componentwise any greater than b:
//  return a.X > b.X || a.Y > b.Y || a.Z > b.Z || a.W > b.W
func (a Vec4) AnyGreater(b Vec4) bool {
	return a.X > b.X || a.Y > b.Y || a.Z > b.Z || a.W > b.W
}

// Clamp clamps each value in the vector to the range of [min, max] and returns
// it.
func (a Vec4) Clamp(min, max float64) Vec4 {
	return Vec4{
		Clamp(a.X, min, max),
		Clamp(a.Y, min, max),
		Clamp(a.Z, min, max),
		Clamp(a.W, min, max),
	}
}

// Radians converts each value in the vector from degrees to radians and
// returns it.
func (a Vec4) Radians() Vec4 {
	return Vec4{
		Radians(a.X),
		Radians(a.Y),
		Radians(a.Z),
		Radians(a.W),
	}
}

// Degrees converts each value in the vector from radians to degrees and
// returns it.
func (a Vec4) Degrees() Vec4 {
	return Vec4{
		Degrees(a.X),
		Degrees(a.Y),
		Degrees(a.Z),
		Degrees(a.W),
	}
}

// Rounded rounds each value in the vector to the nearest whole number and
// returns it.
func (a Vec4) Rounded() Vec4 {
	return Vec4{
		Rounded(a.X),
		Rounded(a.Y),
		Rounded(a.Z),
		Rounded(a.W),
	}
}

// Dot returns the dot product of a and b.
func (a Vec4) Dot(b Vec4) float64 {
	return a.X*b.X + a.Y*b.Y + a.Z*b.Z + a.W*b.W
}

// LengthSq returns the magnitude squared of this vector, useful for comparing
// distances.
func (a Vec4) LengthSq() float64 {
	return a.X*a.X + a.Y*a.Y + a.Z*a.Z + a.W*a.W
}

// Length returns the magnitude of this vector. To avoid a sqrt call when
// strictly comparing distances, LengthSq can be used instead.
func (a Vec4) Length() float64 {
	return math.Sqrt(a.X*a.X + a.Y*a.Y + a.Z*a.Z + a.W*a.W)
}

// Normalized returns the normalized (i.e. length/magnitude == 1) vector of a.
// If the vector's length is zero (and division by zero would occur) then
// [Vec4Zero, false] is returned.
func (a Vec4) Normalized() (v Vec4, ok bool) {
	length := math.Sqrt(a.X*a.X + a.Y*a.Y + a.Z*a.Z + a.W*a.W)
	if Equal(length, 0) {
		return Vec4Zero, false
	}
	return Vec4{
		a.X / length,
		a.Y / length,
		a.Z / length,
		a.W / length,
	}, true
}

// Proj returns a vector representing the projection of vector a onto b.
func (a Vec4) Proj(b Vec4) Vec4 {
	return b.MulScalar(a.Dot(b) / b.LengthSq())
}

// Min returns a vector representing the smallest components of both the
// vectors.
func (a Vec4) Min(b Vec4) Vec4 {
	var r Vec4
	if a.X < b.X {
		r.X = a.X
	} else {
		r.X = b.X
	}
	if a.Y < b.Y {
		r.Y = a.Y
	} else {
		r.Y = b.Y
	}
	if a.Z < b.Z {
		r.Z = a.Z
	} else {
		r.Z = b.Z
	}
	if a.W < b.W {
		r.W = a.W
	} else {
		r.W = b.W
	}
	return r
}

// Max returns a vector representing the largest components of both the
// vectors.
func (a Vec4) Max(b Vec4) Vec4 {
	var r Vec4
	if a.X > b.X {
		r.X = a.X
	} else {
		r.X = b.X
	}
	if a.Y > b.Y {
		r.Y = a.Y
	} else {
		r.Y = b.Y
	}
	if a.Z > b.Z {
		r.Z = a.Z
	} else {
		r.Z = b.Z
	}
	if a.W > b.W {
		r.W = a.W
	} else {
		r.W = b.W
	}
	return r
}

// Lerp returns a vector representing the linear interpolation between the
// vectors a and b. The t parameter is the amount to interpolate (0.0 - 1.0)
// between the vectors.
func (a Vec4) Lerp(b Vec4, t float64) Vec4 {
	return a.Mul(b.MulScalar(t))
}

// Transform transforms this vector by the matrix (vector * matrix) and returns
// the result.
// This is an fully general operation.
func (a Vec4) Transform(b Mat4) Vec4 {
	return Vec4{
		a.X*b[0][0] + a.Y*b[1][0] + a.Z*b[2][0] + a.W*b[3][0],
		a.X*b[0][1] + a.Y*b[1][1] + a.Z*b[2][1] + a.W*b[3][1],
		a.X*b[0][2] + a.Y*b[1][2] + a.Z*b[2][2] + a.W*b[3][2],
		a.X*b[0][3] + a.Y*b[1][3] + a.Z*b[2][3] + a.W*b[3][3],
	}
}

// Quat converts this vector to a quaternion. It's short-hand for:
//  Quat{a.X, a.Y, a.Z, a.W}
func (a Vec4) Quat() Quat {
	return Quat{a.X, a.Y, a.Z, a.W}
}

// Vec3 converts this four-component vector to a three-component one. It's
// short-hand for:
//  Vec3{a.X, a.Y, a.Z}
func (a Vec4) Vec3() Vec3 {
	return Vec3{a.X, a.Y, a.Z}
}

var (
	Vec4One   = Vec4{1, 1, 1, 1}
	Vec4XUnit = Vec4{1, 0, 0, 0}
	Vec4YUnit = Vec4{0, 1, 0, 0}
	Vec4ZUnit = Vec4{0, 0, 1, 0}
	Vec4WUnit = Vec4{0, 0, 0, 1}
	Vec4Zero  = Vec4{0, 0, 0, 0}
)
