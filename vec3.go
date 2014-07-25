// Copyright 2014 The Azul3D Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package lmath

import (
	"fmt"
	"math"
)

// Vec3 represents a 3D vector or point.
type Vec3 struct {
	X, Y, Z float64
}

// String returns an string representation of this vector.
func (a Vec3) String() string {
	return fmt.Sprintf("Vec3(X=%f, Y=%f, Z=%f)", a.X, a.Y, a.Z)
}

// AlmostEquals tells if a == b using the specified epsilon value.
func (a Vec3) AlmostEquals(b Vec3, epsilon float64) bool {
	return AlmostEqual(a.X, b.X, epsilon) && AlmostEqual(a.Y, b.Y, epsilon) && AlmostEqual(a.Z, b.Z, epsilon)
}

// Equals tells if a == b using the default EPSILON value.
func (a Vec3) Equals(b Vec3) bool {
	return a.AlmostEquals(b, EPSILON)
}

// Add performs a componentwise addition of the two vectors, returning a + b.
func (a Vec3) Add(b Vec3) Vec3 {
	return Vec3{a.X + b.X, a.Y + b.Y, a.Z + b.Z}
}

// AddScalar performs a componentwise scalar addition of a + b.
func (a Vec3) AddScalar(b float64) Vec3 {
	return Vec3{a.X + b, a.Y + b, a.Z + b}
}

// Sub performs a componentwise subtraction of the two vectors, returning
// a - b.
func (a Vec3) Sub(b Vec3) Vec3 {
	return Vec3{a.X - b.X, a.Y - b.Y, a.Z - b.Z}
}

// SubScalar performs a componentwise scalar subtraction of a - b.
func (a Vec3) SubScalar(b float64) Vec3 {
	return Vec3{a.X - b, a.Y - b, a.Z - b}
}

// Mul performs a componentwise multiplication of the two vectors, returning
// a * b.
func (a Vec3) Mul(b Vec3) Vec3 {
	return Vec3{a.X * b.X, a.Y * b.Y, a.Z * b.Z}
}

// MulScalar performs a componentwise scalar multiplication of a * b.
func (a Vec3) MulScalar(b float64) Vec3 {
	return Vec3{a.X * b, a.Y * b, a.Z * b}
}

// Div performs a componentwise division of the two vectors, returning a * b.
func (a Vec3) Div(b Vec3) Vec3 {
	return Vec3{a.X / b.X, a.Y / b.Y, a.Z / b.Z}
}

// DivScalar performs a componentwise scalar division of a * b.
func (a Vec3) DivScalar(b float64) Vec3 {
	return Vec3{a.X / b, a.Y / b, a.Z / b}
}

// IsNaN tells if any components of this vector are not an number.
func (a Vec3) IsNaN() bool {
	return math.IsNaN(a.X) || math.IsNaN(a.Y) || math.IsNaN(a.Z)
}

// Less tells if a is componentwise less than b:
//  return a.X < b.X && a.Y < b.Y && a.Z < b.Z
func (a Vec3) Less(b Vec3) bool {
	return a.X < b.X && a.Y < b.Y && a.Z < b.Z
}

// Greater tells if a is componentwise greater than b:
//  return a.X > b.X && a.Y > b.Y && a.Z > b.Z
func (a Vec3) Greater(b Vec3) bool {
	return a.X > b.X && a.Y > b.Y && a.Z > b.Z
}

// AnyLess tells if a is componentwise any less than b:
//  return a.X < b.X || a.Y < b.Y || a.Z < b.Z
func (a Vec3) AnyLess(b Vec3) bool {
	return a.X < b.X || a.Y < b.Y || a.Z < b.Z
}

// AnyGreater tells if a is componentwise any greater than b:
//  return a.X > b.X || a.Y > b.Y || a.Z > b.Z
func (a Vec3) AnyGreater(b Vec3) bool {
	return a.X > b.X || a.Y > b.Y || a.Z > b.Z
}

// Clamp clamps each value in the vector to the range of [min, max] and returns
// it.
func (a Vec3) Clamp(min, max float64) Vec3 {
	return Vec3{
		Clamp(a.X, min, max),
		Clamp(a.Y, min, max),
		Clamp(a.Z, min, max),
	}
}

// Radians converts each value in the vector from degrees to radians and
// returns it.
func (a Vec3) Radians() Vec3 {
	return Vec3{
		Radians(a.X),
		Radians(a.Y),
		Radians(a.Z),
	}
}

// Degrees converts each value in the vector from radians to degrees and
// returns it.
func (a Vec3) Degrees() Vec3 {
	return Vec3{
		Degrees(a.X),
		Degrees(a.Y),
		Degrees(a.Z),
	}
}

// Rounded rounds each value in the vector to the nearest whole number and
// returns it.
func (a Vec3) Rounded() Vec3 {
	return Vec3{
		Rounded(a.X),
		Rounded(a.Y),
		Rounded(a.Z),
	}
}

// Dot returns the dot product of a and b.
func (a Vec3) Dot(b Vec3) float64 {
	return a.X*b.X + a.Y*b.Y + a.Z*b.Z
}

// Inverse returns the inverse (negated) vector -a.
func (a Vec3) Inverse() Vec3 {
	return Vec3{-a.X, -a.Y, -a.Z}
}

// LengthSq returns the magnitude squared of this vector, useful for comparing
// distances.
func (a Vec3) LengthSq() float64 {
	return a.X*a.X + a.Y*a.Y + a.Z*a.Z
}

// Length returns the magnitude of this vector. To avoid a sqrt call when
// strictly comparing distances, LengthSq can be used instead.
func (a Vec3) Length() float64 {
	return math.Sqrt(a.X*a.X + a.Y*a.Y + a.Z*a.Z)
}

// Normalized returns the normalized (i.e. length/magnitude == 1) vector of a.
// If the vector's length is zero (and division by zero would occur) then
// [Vec3Zero, false] is returned.
func (a Vec3) Normalized() (v Vec3, ok bool) {
	length := math.Sqrt(a.X*a.X + a.Y*a.Y + a.Z*a.Z)
	if Equal(length, 0) {
		return Vec3Zero, false
	}
	return Vec3{
		a.X / length,
		a.Y / length,
		a.Z / length,
	}, true
}

// Proj returns a vector representing the projection of vector a onto b.
func (a Vec3) Proj(b Vec3) Vec3 {
	return b.MulScalar(a.Dot(b) / b.LengthSq())
}

// Min returns a vector representing the smallest components of both the
// vectors.
func (a Vec3) Min(b Vec3) Vec3 {
	var r Vec3
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
	return r
}

// Max returns a vector representing the largest components of both the
// vectors.
func (a Vec3) Max(b Vec3) Vec3 {
	var r Vec3
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
	return r
}

// Lerp returns a vector representing the linear interpolation between the
// vectors a and b. The t parameter is the amount to interpolate (0.0 - 1.0)
// between the vectors.
func (a Vec3) Lerp(b Vec3, t float64) Vec3 {
	return a.Mul(b.MulScalar(t))
}

// Cross returns the cross product of the two vectors.
func (a Vec3) Cross(b Vec3) Vec3 {
	return Vec3{
		a.Y*b.Z - b.Y*a.Z,
		b.X*a.Z - a.X*b.Z,
		a.X*b.Y - b.X*a.Y,
	}
}

// Angle returns the unsigned angle between the vectors a and b, in radians.
func (a Vec3) Angle(b Vec3) float64 {
	a, _ = a.Normalized()
	b, _ = b.Normalized()
	var n float64
	if a.Dot(b) < 0 {
		n = a.Add(b).Length() / 2
		return math.Pi - 2.0*math.Asin(math.Min(n, 1))
	}
	n = a.Sub(b).Length() / 2
	return 2.0 * math.Asin(math.Min(n, 1))
}

// SignedAngle returns the signed angle between the vectors a and b, in
// radians.
// The returned angle is positive if the rotation from a to b is clockwise when
// looking in the direction of the reference vector.
func (a Vec3) SignedAngle(b, reference Vec3) float64 {
	a, _ = a.Normalized()
	b, _ = b.Normalized()
	angle := a.Angle(b)
	if a.Cross(b).Dot(reference) < 0.0 {
		angle = -angle
	}
	return angle
}

// TransformMat3 transforms this point vector by the matrix (vector * matrix),
// and returns the result.
// Can operate on orthonormal transformation matrices.
func (a Vec3) TransformMat3(b Mat3) Vec3 {
	return Vec3{
		a.X*b[0][0] + a.Y*b[1][0] + a.Z*b[2][0],
		a.X*b[0][1] + a.Y*b[1][1] + a.Z*b[2][1],
		a.X*b[0][2] + a.Y*b[1][2] + a.Z*b[2][2],
	}
}

// TransformGeneralMat3 transforms this vector by the matrix (vector * matrix)
// without translation component, and returns the result, as a fully general
// operation.
func (a Vec3) TransformGeneralMat3(b Mat3) Vec3 {
	i, _ := b.Inverse()
	return a.TransformMat3(i)
}

// TransformMat4 transforms this point vector by the affine transformation
// matrix (vector * matrix) and returns the result.
// The matrix parameter must be an affine transformation matrix.
func (a Vec3) TransformMat4(b Mat4) Vec3 {
	return Vec3{
		a.X*b[0][0] + a.Y*b[1][0] + a.Z*b[2][0] + b[3][0],
		a.X*b[0][1] + a.Y*b[1][1] + a.Z*b[2][1] + b[3][1],
		a.X*b[0][2] + a.Y*b[1][2] + a.Z*b[2][2] + b[3][2],
	}
}

// TransformVecMat4 transforms this vector (without translation component) by
// the orthonormal matrix and returns the result.
func (a Vec3) TransformVecMat4(b Mat4) Vec3 {
	return Vec3{
		a.X*b[0][0] + a.Y*b[1][0] + a.Z*b[2][0],
		a.X*b[0][1] + a.Y*b[1][1] + a.Z*b[2][1],
		a.X*b[0][2] + a.Y*b[1][2] + a.Z*b[2][2],
	}
}

// TransformGeneralMat4 transforms this vector by the matrix (vector * matrix)
// without translation component, and returns the result, as a fully general
// operation.
func (a Vec3) TransformGeneralMat4(b Mat4) Vec3 {
	i, _ := b.UpperMat3().InverseTransposed()
	return a.TransformMat3(i)
}

// HprToXyz converts Hew, Pitch and Roll rotation to X, Y, and Z axis rotation.
func (v Vec3) HprToXyz() Vec3 {
	return Vec3{v.Y, v.Z, v.X}
}

// XyzToHpr converts X, Y, and Z axis rotation to Hew, Pitch, and Roll
// rotation.
func (v Vec3) XyzToHpr() Vec3 {
	return Vec3{v.Z, v.X, v.Y}
}

var (
	Vec3One   = Vec3{1, 1, 1}
	Vec3XUnit = Vec3{1, 0, 0}
	Vec3YUnit = Vec3{0, 1, 0}
	Vec3ZUnit = Vec3{0, 0, 1}
	Vec3Zero  = Vec3{0, 0, 0}
)
