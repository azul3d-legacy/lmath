// Copyright 2014 The Azul3D Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package lmath

import (
	"fmt"
	"math"
)

// Vec2 represents a 2D vector or point.
type Vec2 struct {
	X, Y float64
}

// String returns an string representation of this vector.
func (a Vec2) String() string {
	return fmt.Sprintf("Vec2(X=%f, Y=%f)", a.X, a.Y)
}

// AlmostEquals tells if a == b using the specified epsilon value.
func (a Vec2) AlmostEquals(b Vec2, epsilon float64) bool {
	return AlmostEqual(a.X, b.X, epsilon) && AlmostEqual(a.Y, b.Y, epsilon)
}

// Equals tells if a == b using the default EPSILON value.
func (a Vec2) Equals(b Vec2) bool {
	return a.AlmostEquals(b, EPSILON)
}

// Add performs a componentwise addition of the two vectors, returning a + b.
func (a Vec2) Add(b Vec2) Vec2 {
	return Vec2{a.X + b.X, a.Y + b.Y}
}

// AddScalar performs a componentwise scalar addition of a + b.
func (a Vec2) AddScalar(b float64) Vec2 {
	return Vec2{a.X + b, a.Y + b}
}

// Sub performs a componentwise subtraction of the two vectors, returning
// a - b.
func (a Vec2) Sub(b Vec2) Vec2 {
	return Vec2{a.X - b.X, a.Y - b.Y}
}

// SubScalar performs a componentwise scalar subtraction of a - b.
func (a Vec2) SubScalar(b float64) Vec2 {
	return Vec2{a.X - b, a.Y - b}
}

// Mul performs a componentwise multiplication of the two vectors, returning
// a * b.
func (a Vec2) Mul(b Vec2) Vec2 {
	return Vec2{a.X * b.X, a.Y * b.Y}
}

// MulScalar performs a componentwise scalar multiplication of a * b.
func (a Vec2) MulScalar(b float64) Vec2 {
	return Vec2{a.X * b, a.Y * b}
}

// Div performs a componentwise division of the two vectors, returning a * b.
func (a Vec2) Div(b Vec2) Vec2 {
	return Vec2{a.X / b.X, a.Y / b.Y}
}

// DivScalar performs a componentwise scalar division of a * b.
func (a Vec2) DivScalar(b float64) Vec2 {
	return Vec2{a.X / b, a.Y / b}
}

// IsNaN tells if any components of this vector are not an number.
func (a Vec2) IsNaN() bool {
	return math.IsNaN(a.X) || math.IsNaN(a.Y)
}

// Less tells if a is componentwise less than b:
//  return a.X < b.X && a.Y < b.Y
func (a Vec2) Less(b Vec2) bool {
	return a.X < b.X && a.Y < b.Y
}

// Greater tells if a is componentwise greater than b:
//  return a.X > b.X && a.Y > b.Y
func (a Vec2) Greater(b Vec2) bool {
	return a.X > b.X && a.Y > b.Y
}

// AnyLess tells if a is componentwise any less than b:
//  return a.X < b.X || a.Y < b.Y
func (a Vec2) AnyLess(b Vec2) bool {
	return a.X < b.X || a.Y < b.Y
}

// AnyGreater tells if a is componentwise any greater than b:
//  return a.X > b.X || a.Y > b.Y
func (a Vec2) AnyGreater(b Vec2) bool {
	return a.X > b.X || a.Y > b.Y
}

// Clamp clamps each value in the vector to the range of [min, max] and returns
// it.
func (a Vec2) Clamp(min, max float64) Vec2 {
	return Vec2{
		Clamp(a.X, min, max),
		Clamp(a.Y, min, max),
	}
}

// Radians converts each value in the vector from degrees to radians and
// returns it.
func (a Vec2) Radians() Vec2 {
	return Vec2{
		Radians(a.X),
		Radians(a.Y),
	}
}

// Degrees converts each value in the vector from radians to degrees and
// returns it.
func (a Vec2) Degrees() Vec2 {
	return Vec2{
		Degrees(a.X),
		Degrees(a.Y),
	}
}

// Rounded rounds each value in the vector to the nearest whole number and
// returns it.
func (a Vec2) Rounded() Vec2 {
	return Vec2{
		Rounded(a.X),
		Rounded(a.Y),
	}
}

// Dot returns the dot product of a and b.
func (a Vec2) Dot(b Vec2) float64 {
	return a.X*b.X + a.Y*b.Y
}

// Inverse returns the inverse (negated) vector -a.
func (a Vec2) Inverse() Vec2 {
	return Vec2{-a.X, -a.Y}
}

// LengthSq returns the magnitude squared of this vector, useful for comparing
// distances.
func (a Vec2) LengthSq() float64 {
	return a.X*a.X + a.Y*a.Y
}

// Length returns the magnitude of this vector. To avoid a sqrt call when
// strictly comparing distances, LengthSq can be used instead.
func (a Vec2) Length() float64 {
	return math.Sqrt(a.X*a.X + a.Y*a.Y)
}

// Normalized returns the normalized (i.e. length/magnitude == 1) vector of a.
// If the vector's length is zero (and division by zero would occur) then
// [Vec2Zero, false] is returned.
func (a Vec2) Normalized() (v Vec2, ok bool) {
	length := math.Sqrt(a.X*a.X + a.Y*a.Y)
	if Equal(length, 0) {
		return Vec2Zero, false
	}
	return Vec2{
		a.X / length,
		a.Y / length,
	}, true
}

// Proj returns a vector representing the projection of vector a onto b.
func (a Vec2) Proj(b Vec2) Vec2 {
	return b.MulScalar(a.Dot(b) / b.LengthSq())
}

// Min returns a vector representing the smallest components of both the
// vectors.
func (a Vec2) Min(b Vec2) Vec2 {
	var r Vec2
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
	return r
}

// Max returns a vector representing the largest components of both the
// vectors.
func (a Vec2) Max(b Vec2) Vec2 {
	var r Vec2
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
	return r
}

// Lerp returns a vector representing the linear interpolation between the
// vectors a and b. The t parameter is the amount to interpolate (0.0 - 1.0)
// between the vectors.
func (a Vec2) Lerp(b Vec2, t float64) Vec2 {
	return a.Mul(b.MulScalar(t))
}

// Angle returns the angle in radians between the two vectors.
func (a Vec2) Angle(b Vec2) float64 {
	return math.Atan2(b.Y-a.Y, b.X-a.X)
}

// TransformVec2 transforms a 2-component point vector by the matrix (without
// translation component) and returns the result.
// This function assumes that the matrix is an affine transformation.
func (a Vec2) TransformVec2(b Mat3) Vec2 {
	return Vec2{
		a.X*b[0][0] + a.Y*b[1][0],
		a.X*b[0][1] + a.Y*b[1][1],
	}
}

// TransformPointVec2 transforms a 2-component point vector by the matrix (with
// translation component) and returns the result.
// This function assumes that the matrix is an affine transformation.
func (a Vec2) TransformPointVec2(b Mat3) Vec2 {
	return Vec2{
		a.X*b[0][0] + a.Y*b[1][0] + b[2][0],
		a.X*b[0][1] + a.Y*b[1][1] + b[2][1],
	}
}

var (
	Vec2One   = Vec2{1, 1}
	Vec2XUnit = Vec2{1, 0}
	Vec2YUnit = Vec2{0, 1}
	Vec2Zero  = Vec2{0, 0}
)
