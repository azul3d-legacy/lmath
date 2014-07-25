// Copyright 2014 The Azul3D Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package lmath

import (
	"fmt"
	"math"
)

// Mat3 represents an 3x3 matrix, indices are in m[row][column] order.
type Mat3 [3][3]float64

// String returns an string representation of this matrix.
func (a Mat3) String() string {
	return fmt.Sprintf("Mat3(\n    %f, %f, %f\n    %f, %f, %f\n    %f, %f, %f\n)", a[0][0], a[0][1], a[0][2], a[1][0], a[1][1], a[1][2], a[2][0], a[2][1], a[2][2])
}

// SetRow sets the values of the specified matrix row to the values in the
// specified three-component vector.
// The row parameter must be 0, 1, or 2 or else a panic will occur.
func (a Mat3) SetRow(row int, values Vec3) Mat3 {
	a[row][0] = values.X
	a[row][1] = values.Y
	a[row][2] = values.Z
	return a
}

// Row returns the values in the specified matrix row as an three-component
// vector.
// The row parameter must be 0, 1, or 2 or else a panic will occur.
func (a Mat3) Row(row int) Vec3 {
	return Vec3{a[row][0], a[row][1], a[row][2]}
}

// SetRowVec2 sets the values in the specified matrix row to the values in the
// specified three-component vector, leaving the third element in the row
// untouched.
// The row parameter must be 0, 1, or 2 or else a panic will occur.
func (a Mat3) SetRowVec2(row int, values Vec2) Mat3 {
	a[row][0] = values.X
	a[row][1] = values.Y
	return a
}

// RowVec2 returns the values in the specified matrix row as an two-component
// vector, the third element of the row is ignored.
// The row parameter must be 0, 1, or 2 or else a panic will occur.
func (a Mat3) RowVec2(row int) Vec2 {
	return Vec2{a[row][0], a[row][1]}
}

// SetCol sets the values in the specified matrix column to the values in the
// specified three-component vector.
// The column parameter must be 0, 1, or 2 or else a panic will occur.
func (a Mat3) SetCol(column int, values Vec3) Mat3 {
	a[0][column] = values.X
	a[1][column] = values.Y
	a[2][column] = values.Z
	return a
}

// Col returns the values in the specified matrix column as an three-component
// vector.
// The column parameter must be 0, 1, or 2 or else a panic will occur.
func (a Mat3) Col(column int) Vec3 {
	return Vec3{a[0][column], a[1][column], a[2][column]}
}

// SetColVec2 sets the values in the specified matrix column to the values in
// the specified two-component vector, leaving the third element in the column
// untouched.
// The column parameter must be 0, 1, or 2 or else an panic will occur.
func (a Mat3) SetColVec2(column int, values Vec2) Mat3 {
	a[0][column] = values.X
	a[1][column] = values.Y
	return a
}

// ColVec2 returns the values in the specified matrix column as an
// two-component vector, ignoring the third element of the column.
// The column parameter must be 0, 1, or 2 or else an panic will occur.
func (a Mat3) ColVec2(column int) Vec2 {
	return Vec2{a[0][column], a[1][column]}
}

// AlmostEquals tells whether a is memberwise equal to b using the specified
// epsilon value.
func (a Mat3) AlmostEquals(b Mat3, epsilon float64) bool {
	return AlmostEqual(a[0][0], b[0][0], epsilon) &&
		AlmostEqual(a[0][1], b[0][1], epsilon) &&
		AlmostEqual(a[0][2], b[0][2], epsilon) &&

		AlmostEqual(a[1][0], b[1][0], epsilon) &&
		AlmostEqual(a[1][1], b[1][1], epsilon) &&
		AlmostEqual(a[1][2], b[1][2], epsilon) &&

		AlmostEqual(a[2][0], b[2][0], epsilon) &&
		AlmostEqual(a[2][1], b[2][1], epsilon) &&
		AlmostEqual(a[2][2], b[2][2], epsilon)
}

// Equals tells whether a is memberwise equal to b using the default EPSILON
// value.
func (a Mat3) Equals(b Mat3) bool {
	return a.AlmostEquals(b, EPSILON)
}

// AddScalar performs memberwise addition a + s and returns the result.
func (a Mat3) AddScalar(s float64) Mat3 {
	return Matrix3(
		a[0][0]+s, a[0][1]+s, a[0][2]+s,
		a[1][0]+s, a[1][1]+s, a[1][2]+s,
		a[2][0]+s, a[2][1]+s, a[2][2]+s,
	)
}

// SubScalar performs memberwise subtraction a - s and returns the result.
func (a Mat3) SubScalar(s float64) Mat3 {
	return Matrix3(
		a[0][0]-s, a[0][1]-s, a[0][2]-s,
		a[1][0]-s, a[1][1]-s, a[1][2]-s,
		a[2][0]-s, a[2][1]-s, a[2][2]-s,
	)
}

// MulScalar performs memberwise multiplication a * s and returns the result.
func (a Mat3) MulScalar(s float64) Mat3 {
	return Matrix3(
		a[0][0]*s, a[0][1]*s, a[0][2]*s,
		a[1][0]*s, a[1][1]*s, a[1][2]*s,
		a[2][0]*s, a[2][1]*s, a[2][2]*s,
	)
}

// DivScalar performs memberwise division a / s and returns the result.
func (a Mat3) DivScalar(s float64) Mat3 {
	return Matrix3(
		a[0][0]/s, a[0][1]/s, a[0][2]/s,
		a[1][0]/s, a[1][1]/s, a[1][2]/s,
		a[2][0]/s, a[2][1]/s, a[2][2]/s,
	)
}

// Add performs memberwise addition a + b and returns the result.
func (a Mat3) Add(b Mat3) Mat3 {
	return Matrix3(
		a[0][0]+b[0][0], a[0][1]+b[0][1], a[0][2]+b[0][2],
		a[1][0]+b[1][0], a[1][1]+b[1][1], a[1][2]+b[1][2],
		a[2][0]+b[2][0], a[2][1]+b[2][1], a[2][2]+b[2][2],
	)
}

// Sub performs in-place subtraction a - b and returns the result.
func (a Mat3) Sub(b Mat3) Mat3 {
	return Matrix3(
		a[0][0]-b[0][0], a[0][1]-b[0][1], a[0][2]-b[0][2],
		a[1][0]-b[1][0], a[1][1]-b[1][1], a[1][2]-b[1][2],
		a[2][0]-b[2][0], a[2][1]-b[2][1], a[2][2]-b[2][2],
	)
}

// Mul returns the result of a * b
func (a Mat3) Mul(b Mat3) Mat3 {
	var out Mat3
	out[0][0] = a[0][0]*b[0][0] + a[0][1]*b[1][0] + a[0][2]*b[2][0]
	out[0][1] = a[0][0]*b[0][1] + a[0][1]*b[1][1] + a[0][2]*b[2][1]
	out[0][2] = a[0][0]*b[0][2] + a[0][1]*b[1][2] + a[0][2]*b[2][2]
	out[1][0] = a[1][0]*b[0][0] + a[1][1]*b[1][0] + a[1][2]*b[2][0]
	out[1][1] = a[1][0]*b[0][1] + a[1][1]*b[1][1] + a[1][2]*b[2][1]
	out[1][2] = a[1][0]*b[0][2] + a[1][1]*b[1][2] + a[1][2]*b[2][2]
	out[2][0] = a[2][0]*b[0][0] + a[2][1]*b[1][0] + a[2][2]*b[2][0]
	out[2][1] = a[2][0]*b[0][1] + a[2][1]*b[1][1] + a[2][2]*b[2][1]
	out[2][2] = a[2][0]*b[0][2] + a[2][1]*b[1][2] + a[2][2]*b[2][2]
	return out
}

// Transposed returns the transposed version of matrix a.
func (a Mat3) Transposed() Mat3 {
	a01 := a[0][1]
	a[0][1] = a[1][0]
	a[1][0] = a01

	a02 := a[0][2]
	a[0][2] = a[2][0]
	a[2][0] = a02

	a12 := a[1][2]
	a[1][2] = a[2][1]
	a[2][1] = a12
	return a
}

// Determinant calculates and returns the determinant of the matrix.
func (a Mat3) Determinant() float64 {
	det2 := func(e00, e01, e10, e11 float64) float64 {
		return e00*e11 - e10*e01
	}

	return a[0][0]*det2(a[1][1], a[1][2], a[2][1], a[2][2]) - a[0][1]*det2(a[1][0], a[1][2], a[2][0], a[2][2]) + a[0][2]*det2(a[1][0], a[1][1], a[2][0], a[2][1])
}

// MulQuat multiplies the matrix by the specified quaternion rotation and
// returns the result,
func (a Mat3) MulQuat(b Quat) Mat3 {
	return a.Mul(b.ExtractToMat3())
}

// Inverse returns the inverse of the matrix a. This is a fully general
// operation and has no requirements about the transformation represented by
// this matrix.
// Returned is result=Mat3Identity, ok=false if the matrix was singular and
// could not be inverted.
func (a Mat3) Inverse() (result Mat3, ok bool) {
	det := a.Determinant()

	// We throw the value out only if it's smaller than our "small"
	// tolerence squared.  This helps reduce overly-sensitive
	// rejections.
	if Equal(det, EPSILON*EPSILON) {
		return Mat3Identity, false
	}

	det = 1.0 / det

	det2 := func(e00, e01, e10, e11 float64) float64 {
		return e00*e11 - e10*e01
	}

	result[0][0] = det * det2(a[1][1], a[1][2], a[2][1], a[2][2])
	result[1][0] = -det * det2(a[1][0], a[1][2], a[2][0], a[2][2])
	result[2][0] = det * det2(a[1][0], a[1][1], a[2][0], a[2][1])

	result[0][1] = -det * det2(a[0][1], a[0][2], a[2][1], a[2][2])
	result[1][1] = det * det2(a[0][0], a[0][2], a[2][0], a[2][2])
	result[2][1] = -det * det2(a[0][0], a[0][1], a[2][0], a[2][1])

	result[0][2] = det * det2(a[0][1], a[0][2], a[1][1], a[1][2])
	result[1][2] = -det * det2(a[0][0], a[0][2], a[1][0], a[1][2])
	result[2][2] = det * det2(a[0][0], a[0][1], a[1][0], a[1][1])
	return result, true
}

// InverseTransposed simultaneously inverts and transposes the matrix and
// returns it.
// Returned is result=Mat3Identity, ok=false if the matrix was singular and
// could not be inverted.
func (a Mat3) InverseTransposed() (result Mat3, ok bool) {
	det := a.Determinant()

	// We throw the value out only if it's smaller than our "small"
	// tolerence squared.  This helps reduce overly-sensitive
	// rejections.
	if Equal(det, EPSILON*EPSILON) {
		return Mat3Identity, false
	}

	det = 1.0 / det

	det2 := func(e00, e01, e10, e11 float64) float64 {
		return e00*e11 - e10*e01
	}

	result[0][0] = det * det2(a[1][1], a[1][2], a[2][1], a[2][2])
	result[0][1] = -det * det2(a[1][0], a[1][2], a[2][0], a[2][2])
	result[0][2] = det * det2(a[1][0], a[1][1], a[2][0], a[2][1])

	result[1][0] = -det * det2(a[0][1], a[0][2], a[2][1], a[2][2])
	result[1][1] = det * det2(a[0][0], a[0][2], a[2][0], a[2][2])
	result[1][2] = -det * det2(a[0][0], a[0][1], a[2][0], a[2][1])

	result[2][0] = det * det2(a[0][1], a[0][2], a[1][1], a[1][2])
	result[2][1] = -det * det2(a[0][0], a[0][2], a[1][0], a[1][2])
	result[2][2] = det * det2(a[0][0], a[0][1], a[1][0], a[1][1])
	return result, true
}

// unwindYUpRotation extracts the rotation about the x, y, and z axes from the
// given hpr & scale matrix. Adjusts the matrix to eliminate the rotation.
// This function assumes the matrix is stored in a right-handed Y-up coordinate
// system.
func (a Mat3) unwindYUpRotation(cs CoordSys) (updated Mat3, hpr Vec3) {
	// Extract the axes from the matrix.
	x := a.Row(0)
	y := a.Row(1)
	z := a.Row(2)

	// Project Z into the XZ plane.
	xz, _ := Vec2{z.X, z.Z}.Normalized()

	// Compute the rotation about the +Y (up) axis.  This is yaw / heading.
	heading := math.Atan2(xz.X, xz.Y)

	// Unwind the heading, and continue.
	rotY := Mat3FromAxisAngle(Vec3{0, 1, 0}, -heading, CoordSysYUpRight)

	x = x.TransformMat3(rotY)
	y = y.TransformMat3(rotY)
	z = z.TransformMat3(rotY)

	// Project the rotated Z into the YZ plane.
	yz, _ := Vec2{z.Y, z.Z}.Normalized()

	// Compute the rotation about the +X (right) axis.  This is pitch.
	pitch := -math.Atan2(yz.X, yz.Y)

	// Unwind the pitch.
	rotX := Mat3FromAxisAngle(Vec3{1, 0, 0}, -pitch, CoordSysYUpRight)

	x = x.TransformMat3(rotX)
	y = y.TransformMat3(rotX)
	z = z.TransformMat3(rotX)

	// Project the rotated X onto the XY plane.
	xy, _ := Vec2{x.X, x.Y}.Normalized()

	// Compute the rotation about the +Z (back) axis.  This is roll.
	roll := -math.Atan2(xy.Y, xy.X)

	// Unwind the roll from the axes, and continue.
	rotZ := Mat3FromAxisAngle(Vec3{0, 0, 1}, roll, CoordSysZUpRight)

	x = x.TransformMat3(rotZ)
	y = y.TransformMat3(rotZ)
	z = z.TransformMat3(rotZ)

	// Reset the matrix to reflect the unwinding.
	a = a.SetRow(0, x)
	a = a.SetRow(1, y)
	a = a.SetRow(2, z)

	// Return the three rotation components.
	return a, Vec3{heading, pitch, roll}
}

// unwindZUpRotation extracts the rotation about the x, y, and z axes from the
// given hpr & scale matrix. Adjusts the matrix to eliminate the rotation.
// This function assumes the matrix is stored in a right-handed Z-up coordinate
// system.
func (a Mat3) unwindZUpRotation(cs CoordSys) (updated Mat3, hpr Vec3) {
	// Extract the axes from the matrix.
	x := a.Row(0)
	y := a.Row(1)
	z := a.Row(2)

	// Project Y into the XY plane.
	xy, _ := Vec2{y.X, y.Y}.Normalized()

	// Compute the rotation about the +Z (up) axis.  This is yaw / heading.
	heading := -math.Atan2(xy.X, xy.Y)

	// Unwind the heading, and continue.
	rotZ := Mat3FromAxisAngle(Vec3{0, 0, 1}, -heading, CoordSysZUpRight)

	x = x.TransformMat3(rotZ)
	y = y.TransformMat3(rotZ)
	z = z.TransformMat3(rotZ)

	// Project the rotated Y into the YZ plane.
	yz, _ := Vec2{y.Y, y.Z}.Normalized()

	// Compute the rotation about the +X (right) axis.  This is pitch.
	pitch := math.Atan2(yz.Y, yz.X)

	// Unwind the pitch.
	rotX := Mat3FromAxisAngle(Vec3{1, 0, 0}, -pitch, CoordSysZUpRight)

	x = x.TransformMat3(rotX)
	y = y.TransformMat3(rotX)
	z = z.TransformMat3(rotX)

	// Project X into the XZ plane.
	xz, _ := Vec2{x.X, x.Z}.Normalized()

	// Compute the rotation about the -Y (back) axis.  This is roll.
	roll := -math.Atan2(xz.Y, xz.X)

	// Unwind the roll from the axes, and continue.
	rotY := Mat3FromAxisAngle(Vec3{0, 1, 0}, -roll, CoordSysZUpRight)

	x = x.TransformMat3(rotY)
	y = y.TransformMat3(rotY)
	z = z.TransformMat3(rotY)

	// Reset the matrix to reflect the unwinding.
	a = a.SetRow(0, x)
	a = a.SetRow(1, y)
	a = a.SetRow(2, z)

	// Return the three rotation components.
	return a, Vec3{heading, pitch, roll}
}

// Decompose extracts out the scaling, shearing, and hew/pitch/roll components
// from the composed rotation matrix.
// The coordinate system must be one of: CoordSysZUpRight, CoordSysZUpLeft, or
// CoordSysYUpLeft.
func (a Mat3) Decompose(cs CoordSys) (scale, shear, hpr Vec3) {
	// Extract the rotation and scale, according to the coordinate
	// system of choice.
	switch cs {
	case CoordSysZUpRight:
		a, hpr = a.unwindZUpRotation(cs)

	case CoordSysZUpLeft:
		a[0][2] = -a[0][2]
		a[1][2] = -a[1][2]
		a[2][0] = -a[2][0]
		a[2][1] = -a[2][1]
		a, hpr = a.unwindZUpRotation(cs)
		hpr.X = -hpr.X
		hpr.Z = -hpr.Z

	case CoordSysYUpLeft:
		a[0][2] = -a[0][2]
		a[1][2] = -a[1][2]
		a[2][0] = -a[2][0]
		a[2][1] = -a[2][1]
		a, hpr = a.unwindYUpRotation(cs)

	default:
		panic(fmt.Sprintf("Decompose(): Unexpected coordinate system %d", cs))
	}

	scale = Vec3{a[0][0], a[1][1], a[2][2]}

	// Normalize the scale out of the shear components, and return the
	// shear.
	if !Equal(scale.X, 0) {
		a[0][1] /= scale.X
		a[0][2] /= scale.X
	}
	if !Equal(scale.Y, 0) {
		a[1][0] /= scale.Y
		a[1][2] /= scale.Y
	}
	if !Equal(scale.Z, 0) {
		a[2][0] /= scale.Z
		a[2][1] /= scale.Z
	}

	shear = Vec3{
		a[0][1] + a[1][0],
		a[2][0] + a[0][2],
		a[2][1] + a[1][2],
	}
	return
}

// IsNan tells if any components of this matrix are not an number.
func (a Mat3) IsNaN() bool {
	return math.IsNaN(a[0][0]) || math.IsNaN(a[0][1]) || math.IsNaN(a[0][2]) ||
		math.IsNaN(a[1][0]) || math.IsNaN(a[1][1]) || math.IsNaN(a[1][2]) ||
		math.IsNaN(a[2][0]) || math.IsNaN(a[2][1]) || math.IsNaN(a[2][2])
}

// Matrix3 returns an new Mat3 given the specified matrix components.
func Matrix3(m00, m01, m02, m10, m11, m12, m20, m21, m22 float64) Mat3 {
	return Mat3{
		{m00, m01, m02},
		{m10, m11, m12},
		{m20, m21, m22},
	}
}

// Mat3FromAxisAngle returns a rotation matrix that will rotate by the given
// angle in radians counterclockwise about the indicated axis.
func Mat3FromAxisAngle(axis Vec3, angle float64, cs CoordSys) Mat3 {
	axis, _ = axis.Normalized()

	// In a left-handed coordinate system, counterclockwise is the
	// other direction.
	if cs.LeftHanded() {
		angle = -angle
	}

	s := math.Sin(angle)
	c := math.Cos(angle)
	t := 1.0 - c

	t0 := t * axis.X
	t1 := t * axis.Y
	t2 := t * axis.Z
	s0 := s * axis.X
	s1 := s * axis.Y
	s2 := s * axis.Z

	return Matrix3(
		t0*axis.X+c,
		t0*axis.Y+s2,
		t0*axis.Z-s1,

		t1*axis.X-s2,
		t1*axis.Y+c,
		t1*axis.Z+s0,

		t2*axis.X+s1,
		t2*axis.Y-s0,
		t2*axis.Z+c,
	)
}

// Mat3Compose composes a transformation matrix that applies the given scaling,
// shearing, and hew/pitch/roll euler rotation values for the given coordinate
// system.
func Mat3Compose(scale, shear, hpr Vec3, cs CoordSys) Mat3 {
	a := Mat3FromScaleShear(scale, shear, cs)

	if !Equal(hpr.Z, 0) {
		r := Mat3FromAxisAngle(cs.Forward(), hpr.Z, cs)
		a = a.Mul(r)
	}
	if !Equal(hpr.Y, 0) {
		r := Mat3FromAxisAngle(cs.Right(), hpr.Y, cs)
		a = a.Mul(r)
	}
	if !Equal(hpr.X, 0) {
		r := Mat3FromAxisAngle(cs.Up(), hpr.X, cs)
		a = a.Mul(r)
	}
	return a
}

// Mat3FromScaleShear returns a matrix that will apply the given scaling and
// shearing values along their respective axis in the specified coordinate
// system.
// A panic will occur if the coordinate system is invalid.
func Mat3FromScaleShear(scale, shear Vec3, cs CoordSys) Mat3 {
	// We have to match the placement of the shear components in the
	// matrix to the way we extract out the rotation in Decompose(). Therefore,
	// the shear is sensitive to the coordinate system.
	switch cs {
	case CoordSysZUpRight:
		return Matrix3(
			scale.X, shear.X*scale.X, 0,
			0, scale.Y, 0,
			shear.Y*scale.Z, shear.Z*scale.Z, scale.Z,
		)

	case CoordSysZUpLeft:
		return Matrix3(
			scale.X, shear.X*scale.X, 0,
			0, scale.Y, 0,
			-shear.Y*scale.Z, -shear.Z*scale.Z, scale.Z,
		)

	case CoordSysYUpRight:
		return Matrix3(
			scale.X, 0, shear.Y*scale.X,
			shear.X*scale.Y, scale.Y, shear.Z*scale.Y,
			0, 0, scale.Z,
		)

	case CoordSysYUpLeft:
		return Matrix3(
			scale.X, 0, -shear.Y*scale.X,
			shear.X*scale.Y, scale.Y, -shear.Z*scale.Y,
			0, 0, scale.Z,
		)

	default:
		panic(fmt.Sprintf("Mat3FromScaleShear(): Invalid coordinate system %d", cs))
	}
}

// Mat3FromTranslation returns a matrix that will apply the given translation
// vector.
func Mat3FromTranslation(translation Vec2) Mat3 {
	return Matrix3(
		1, 0, 0,
		0, 1, 0,
		translation.X, translation.Y, 1,
	)
}

var (
	Mat3Identity = Matrix3(
		1, 0, 0,
		0, 1, 0,
		0, 0, 1,
	)

	Mat3YToZUp = Matrix3(
		1, 0, 0,
		0, 0, 1,
		0, -1, 0,
	)

	Mat3ZToYUp = Matrix3(
		1, 0, 0,
		0, 0, -1,
		0, 1, 0,
	)

	Mat3FlipY = Matrix3(
		1, 0, 0,
		0, -1, 0,
		0, 0, 1,
	)

	Mat3FlipZ = Matrix3(
		1, 0, 0,
		0, 1, 0,
		0, 0, -1,
	)

	Mat3LZToRY = Mat3FlipY.Mul(Mat3ZToYUp)
	Mat3LYToRZ = Mat3FlipZ.Mul(Mat3YToZUp)
)
